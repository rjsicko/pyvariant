import os.path
import sys
from ftplib import FTP
from pathlib import Path
from string import punctuation
from typing import List, Union

import pandas as pd
from appdirs import user_data_dir
from gtfparse import read_gtf
from pyfaidx import Fasta

from .constants import CACHE_DIR_ENV, NAME
from .files import bgzip, is_bgzipped, read_fasta
from .normalize import normalize_df
from .utils import strip_version

# Dictionary used to replace punctuation in a string
PUNCTUATION_TO_UNDERSCORE = str.maketrans(punctuation + " ", "_" * len(punctuation + " "))

# Ensembl FTP URL
ENSEMBL_FTP_SERVER = "ftp.ensembl.org"

# FASTA directory example: pub/release-100/fasta/homo_sapiens
FASTA_SUBDIR_TEMPLATE = "pub/release-{release}/fasta/{species}/{type}"

# FASTA file name example: Homo_sapiens.GRCh38.100.dna.toplevel.fa.gz
# NOTE: The file name format changed after Ensembl release 75
FASTA_FILENAME_TEMPLATE_OLD = {
    "cdna": "{species}.{reference}.{release}.{type}.all.fa.gz",
    "dna": "{species}.{reference}.{release}.{type}.toplevel.fa.gz",
    "ncrna": "{species}.{reference}.{release}.{type}.fa.gz",
    "pep": "{species}.{reference}.{release}.{type}.all.fa.gz",
}
FASTA_FILENAME_TEMPLATE_NEW = {
    "cdna": "{species}.{reference}.{type}.all.fa.gz",
    "dna": "{species}.{reference}.{type}.toplevel.fa.gz",
    "ncrna": "{species}.{reference}.{type}.fa.gz",
    "pep": "{species}.{reference}.{type}.all.fa.gz",
}

# GTF annotation directory example: pub/release-100/gtf/homo_sapiens
GTF_SUBDIR_TEMPLATE = "pub/release-{release}/gtf/{species}"

# GTF annotation file example: Homo_sapiens.GRCh38.100.gtf.gz
GTF_FILENAME_TEMPLATE = "{species}.{reference}.{release}.gtf.gz"


class EnsemblCache:
    """Class for managing Ensembl files."""

    def __init__(self, species: str, release: int, cache_dir: str = ""):
        self.species = normalize_species(species)
        self.release = normalize_release(release)
        self.reference = reference_by_release(self.release)

        if cache_dir:
            self.cache_dir = os.path.abspath(cache_dir)
        else:
            self.cache_dir = get_cache_dir()

    def install(
        self,
        clean: bool = True,
        recache: bool = False,
        redownload: bool = False,
        restrict_genes: List[str] = [],
    ):
        """Download missing data, process, and cache."""

        # Create the cache directory structure
        self.make_release_cache_dir()

        # Download each FASTA file
        for fasta, index, downloadf in [
            (
                self.local_cdna_fasta_filepath,
                self.local_cdna_index_filepath,
                self.download_cdna_fasta,
            ),
            (self.local_dna_fasta_filepath, self.local_dna_index_filepath, self.download_dna_fasta),
            (self.local_pep_fasta_filepath, self.local_pep_index_filepath, self.download_pep_fasta),
            (
                self.local_ncrna_fasta_filepath,
                self.local_ncrna_index_filepath,
                self.download_ncrna_fasta,
            ),
        ]:
            if not os.path.exists(fasta) or redownload:
                downloadf()

        # Re-compress (if required) and index each FASTA file
        for fasta, index, indexf in [
            (self.local_cdna_fasta_filepath, self.local_cdna_index_filepath, self.index_cdna_fasta),
            (self.local_dna_fasta_filepath, self.local_dna_index_filepath, self.index_dna_fasta),
            (self.local_pep_fasta_filepath, self.local_pep_index_filepath, self.index_pep_fasta),
            (
                self.local_ncrna_fasta_filepath,
                self.local_ncrna_index_filepath,
                self.index_ncrna_fasta,
            ),
        ]:
            # NOTE: pyfaidx only supports compressed FASTA in BGZF format. Ensembl FASTA comes in
            # GZ format, so we need to compress the data with bgzip.
            if not is_bgzipped(fasta):
                bgzip(fasta)
            if not os.path.exists(index) or redownload:
                indexf()

        # Download the GTF file
        gtf_missing = not os.path.exists(self.local_gtf_filepath)
        cache_missing = not os.path.exists(self.local_gtf_cache_filepath)
        if (gtf_missing and cache_missing) or redownload:
            self.download_gtf()

        # Process and cache the GTF file
        if cache_missing or recache:
            df = read_gtf(self.local_gtf_filepath)
            if restrict_genes:
                df = df[df["gene_name"].isin(restrict_genes)]

            df = normalize_df(df)
            self.cache_df(df)
            if clean:
                self.delete_gtf()

    # ---------------------------------------------------------------------------------------------
    # Release cache directory
    # ---------------------------------------------------------------------------------------------
    @property
    def release_cache_dir(self) -> str:
        """Local path to the directory containing the release data files."""
        return os.path.join(self.cache_dir, self.species, self.reference, str(self.release))

    def make_release_cache_dir(self):
        """Create the release cache directory, if it doesn't already exist."""
        Path(self.release_cache_dir).mkdir(exist_ok=True, parents=True)

    # ---------------------------------------------------------------------------------------------
    # Remote FASTA file paths
    # ---------------------------------------------------------------------------------------------
    @property
    def remote_cdna_fasta_subdir(self) -> str:
        """Path to the directory on the FTP server that contains the cDNA FASTA file."""
        return self._remote_fasta_subdir("cdna")

    @property
    def remote_dna_fasta_subdir(self) -> str:
        """Path to the directory on the FTP server that contains the DNA FASTA file."""
        return self._remote_fasta_subdir("dna")

    @property
    def remote_ncrna_fasta_subdir(self) -> str:
        """Path to the directory on the FTP server that contains the ncRNA FASTA file."""
        return self._remote_fasta_subdir("ncrna")

    @property
    def remote_pep_fasta_subdir(self) -> str:
        """Path to the directory on the FTP server that contains the peptide FASTA file."""
        return self._remote_fasta_subdir("pep")

    def _remote_fasta_subdir(self, fasta_type: str) -> str:
        return FASTA_SUBDIR_TEMPLATE.format(
            release=self.release, species=self.species, type=fasta_type
        )

    @property
    def remote_cdna_fasta_filename(self) -> str:
        """Name of the cDNA FASTA file on the FTP server."""
        return self._remote_fasta_filename("cdna")

    @property
    def remote_dna_fasta_filename(self) -> str:
        """Name of the DNA FASTA file on the FTP server."""
        return self._remote_fasta_filename("dna")

    @property
    def remote_ncrna_fasta_filename(self) -> str:
        """Name of the ncRNA FASTA file on the FTP server."""
        return self._remote_fasta_filename("ncrna")

    @property
    def remote_pep_fasta_filename(self) -> str:
        """Name of the peptide FASTA file on the FTP server."""
        return self._remote_fasta_filename("pep")

    def _remote_fasta_filename(self, fasta_type: str) -> str:
        if self.release <= 75:
            template = FASTA_FILENAME_TEMPLATE_OLD[fasta_type]
            return template.format(
                species=self.species.capitalize(),
                reference=self.reference,
                release=self.release,
                type=fasta_type,
            )
        else:
            template = FASTA_FILENAME_TEMPLATE_NEW[fasta_type]
            return template.format(
                species=self.species.capitalize(), reference=self.reference, type=fasta_type
            )

    # ---------------------------------------------------------------------------------------------
    # Local FASTA file paths
    # ---------------------------------------------------------------------------------------------
    @property
    def local_cdna_fasta_filepath(self) -> str:
        """Local name of the cDNA FASTA file."""
        return self._local_fasta_filepath(self.remote_cdna_fasta_filename)

    @property
    def local_dna_fasta_filepath(self) -> str:
        """Local name of the DNA FASTA file."""
        return self._local_fasta_filepath(self.remote_dna_fasta_filename)

    @property
    def local_ncrna_fasta_filepath(self) -> str:
        """Local name of the ncRNA FASTA file."""
        return self._local_fasta_filepath(self.remote_ncrna_fasta_filename)

    @property
    def local_pep_fasta_filepath(self) -> str:
        """Local name of the peptide FASTA file."""
        return self._local_fasta_filepath(self.remote_pep_fasta_filename)

    def _local_fasta_filepath(self, remote_fasta_filename: str) -> str:
        return os.path.join(self.release_cache_dir, remote_fasta_filename)

    # ---------------------------------------------------------------------------------------------
    # Download FASTA files
    # ---------------------------------------------------------------------------------------------
    def download_cdna_fasta(self):
        """Download the cDNA FASTA file from the FTP server."""
        self._ftp_download(
            ENSEMBL_FTP_SERVER,
            self.remote_cdna_fasta_subdir,
            self.remote_cdna_fasta_filename,
            self.local_cdna_fasta_filepath,
        )

    def download_dna_fasta(self):
        """Download the DNA FASTA file from the FTP server."""
        self._ftp_download(
            ENSEMBL_FTP_SERVER,
            self.remote_dna_fasta_subdir,
            self.remote_dna_fasta_filename,
            self.local_dna_fasta_filepath,
        )

    def download_ncrna_fasta(self):
        """Download the ncRNA FASTA file from the FTP server."""
        self._ftp_download(
            ENSEMBL_FTP_SERVER,
            self.remote_ncrna_fasta_subdir,
            self.remote_ncrna_fasta_filename,
            self.local_ncrna_fasta_filepath,
        )

    def download_pep_fasta(self):
        """Download the peptide FASTA file from the FTP server."""
        self._ftp_download(
            ENSEMBL_FTP_SERVER,
            self.remote_pep_fasta_subdir,
            self.remote_pep_fasta_filename,
            self.local_pep_fasta_filepath,
        )

    # ---------------------------------------------------------------------------------------------
    # Local FASTA index paths
    # ---------------------------------------------------------------------------------------------
    @property
    def local_cdna_index_filepath(self) -> str:
        """Local name of the cDNA FASTA file."""
        return self._fasta_index_path(self.local_cdna_fasta_filepath)

    @property
    def local_dna_index_filepath(self) -> str:
        """Local name of the DNA FASTA file."""
        return self._fasta_index_path(self.local_dna_fasta_filepath)

    @property
    def local_ncrna_index_filepath(self) -> str:
        """Local name of the ncRNA FASTA file."""
        return self._fasta_index_path(self.local_ncrna_fasta_filepath)

    @property
    def local_pep_index_filepath(self) -> str:
        """Local name of the peptide FASTA file."""
        return self._fasta_index_path(self.local_pep_fasta_filepath)

    def _fasta_index_path(self, local_fasta_filename: str) -> str:
        """Return the path to the FASTA index file."""
        return local_fasta_filename + ".fai"

    # ---------------------------------------------------------------------------------------------
    # Index FASTA files
    # ---------------------------------------------------------------------------------------------
    def index_cdna_fasta(self):
        """(Re)build the index file for the cDNA Fasta."""
        return self._index_fasta(self.local_cdna_fasta_filepath)

    def index_dna_fasta(self):
        """(Re)build the index file for the DNA Fasta."""
        return self._index_fasta(self.local_dna_fasta_filepath)

    def index_ncrna_fasta(self):
        """(Re)build the index file for the cDNA Fasta."""
        return self._index_fasta(self.local_ncrna_fasta_filepath)

    def index_pep_fasta(self):
        """(Re)build the index file for the peptide Fasta."""
        return self._index_fasta(self.local_pep_fasta_filepath)

    def _index_fasta(self, local_fasta_filename: str):
        print(f"Indexing {local_fasta_filename}...", file=sys.stderr)
        _ = Fasta(
            local_fasta_filename,
            key_function=strip_version,
            as_raw=True,
            sequence_always_upper=True,
            build_index=True,
            rebuild=True,
        )

    # ---------------------------------------------------------------------------------------------
    # Load FASTA files
    # ---------------------------------------------------------------------------------------------
    def load_cdna_fasta(self) -> Fasta:
        """Load and return the cDNA a `pyfaidx.Fasta` object."""
        return read_fasta(self.local_cdna_fasta_filepath)

    def load_dna_fasta(self) -> Fasta:
        """Load and return the DNA a `pyfaidx.Fasta` object."""
        return read_fasta(self.local_dna_fasta_filepath)

    def load_ncrna_fasta(self) -> Fasta:
        """Load and return the ncDNA a `pyfaidx.Fasta` object."""
        return read_fasta(self.local_ncrna_fasta_filepath)

    def load_pep_fasta(self) -> Fasta:
        """Load and return the peptide a `pyfaidx.Fasta` object."""
        return read_fasta(self.local_pep_fasta_filepath)

    # ---------------------------------------------------------------------------------------------
    # GTF files
    # ---------------------------------------------------------------------------------------------
    @property
    def remote_gtf_subdir(self) -> str:
        """Path to the directory on the FTP server that contains the GTF file."""
        return GTF_SUBDIR_TEMPLATE.format(release=self.release, species=self.species)

    @property
    def remote_gtf_filename(self) -> str:
        """Name of the GTF file on the FTP server."""
        return GTF_FILENAME_TEMPLATE.format(
            species=self.species.capitalize(), reference=self.reference, release=self.release
        )

    @property
    def local_gtf_filepath(self) -> str:
        """Local name of the GTF file."""
        return os.path.join(self.release_cache_dir, self.remote_gtf_filename)

    @property
    def local_gtf_cache_filepath(self) -> str:
        """Local name of the cache file."""
        return os.path.join(
            self.release_cache_dir, self.remote_gtf_filename.replace(".gtf.gz", ".pickle")
        )

    def download_gtf(self):
        """Download the GTF file from the FTP server."""
        self._ftp_download(
            ENSEMBL_FTP_SERVER,
            self.remote_gtf_subdir,
            self.remote_gtf_filename,
            self.local_gtf_filepath,
        )

    def delete_gtf(self):
        """Delete the GTF file from the local filesystem."""
        print(f"Removing {self.local_gtf_filepath}", file=sys.stderr)
        os.remove(self.local_gtf_filepath)

    # ---------------------------------------------------------------------------------------------
    # Database cache
    # ---------------------------------------------------------------------------------------------
    def cache_df(self, df: pd.DataFrame) -> pd.DataFrame:
        """Convert the Ensembl GTF to a pandas DataFrame and cache it."""
        print(
            f"Converting {self.local_gtf_filepath} to {self.local_gtf_cache_filepath}",
            file=sys.stderr,
        )
        df.to_pickle(self.local_gtf_cache_filepath)

    def load_df(self) -> pd.DataFrame:
        """Return the cached Ensembl GTF data."""
        return pd.read_pickle(self.local_gtf_cache_filepath)

    # ---------------------------------------------------------------------------------------------
    # Generic functions
    # ---------------------------------------------------------------------------------------------
    def _ftp_download(self, server: str, subdir: str, remote_file: str, local_file: str):
        """Download a file from an FTP server."""
        try:
            ftp = FTP(server)
            print(f"Connecting to {server}", file=sys.stderr)
            ftp.login()
            url = "ftp://" + server + "/" + subdir + "/" + remote_file
            print(f"Downloading {url} to {local_file}", file=sys.stderr)
            ftp.cwd(subdir)
            self.make_release_cache_dir()
            with open(local_file, "wb") as fh:
                ftp.retrbinary(f"RETR {remote_file}", fh.write)
            print("Download successful", file=sys.stderr)
            ftp.quit()
        except Exception as exc:
            raise RuntimeError(f"Download failed: {exc}")


def get_cache_dir() -> str:
    f"""Get the cache root directory. If the environmental variable '{CACHE_DIR_ENV}' is set, this
    package will use that as the directory. Otherwise, this package will use the default appdata
    dir for the user (platform dependant).
    """
    try:
        return os.environ[CACHE_DIR_ENV]
    except KeyError:
        return os.path.join(user_data_dir(), NAME)


def normalize_release(release: Union[float, int, str]) -> int:
    """Normalize a release number."""
    return int(release)


def normalize_species(species: str) -> str:
    """Normalize a species name."""
    return species.lower().translate(PUNCTUATION_TO_UNDERSCORE)


def reference_by_release(release: int) -> str:
    """Given the Ensembl release number, return the reference name.

    Examples:
        >>> reference_by_release(69)
        'GRCh37'
        >>> reference_by_release(100)
        'GRCh38'
    """
    if release == 54:
        return "GRCh36"
    elif 54 < release <= 75:
        return "GRCh37"
    elif 75 < release:
        return "GRCh38"
    else:
        raise ValueError(f"Unknown reference for release '{release}'")
