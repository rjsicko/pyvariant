from ftplib import FTP
from logzero import logger
import os.path
from pathlib import Path

from gtfparse import read_gtf
import pandas

# default cache directory
DEFAULT_CACHE_DIR = "."

# Ensembl FTP URL
ENSEMBL_FTP_SERVER = "ftp.ensembl.org"

# FASTA directory example: pub/release-100/fasta/homo_sapiens
FASTA_SUBDIR_TEMPLATE = "pub/release-{release}/fasta/{species}/{type}"

# FASTA example: Homo_sapiens.GRCh38.100.dna.toplevel.fa.gz
# the file name format changed after Ensembl release 75
FASTA_FILENAME_TEMPLATE_DNA_OLD = "{species}.{reference}.{release}.dna.toplevel.fa.gz"
FASTA_FILENAME_TEMPLATE_DNA_NEW = "{species}.{reference}.dna.toplevel.fa.gz"

# GTF annotation directory example: pub/release-100/gtf/homo_sapiens
GTF_SUBDIR_TEMPLATE = "pub/release-{release}/gtf/{species}"

# GTF annotation file example: Homo_sapiens.GRCh38.100.gtf.gz
GTF_FILENAME_TEMPLATE = "{species}.{reference}.{release}.gtf.gz"


class DataCache:
    """Class for managing the data files required by this package."""

    def __init__(
        self, species: str, reference: str, release: int, cache_dir: str = DEFAULT_CACHE_DIR
    ):
        self.species = species
        self.reference = reference
        self.release = release
        self.cache_dir = os.path.abspath(cache_dir)

    @property
    def release_cache_dir(self) -> str:
        """Local path to the directory containing the release data files."""
        return os.path.join(self.cache_dir, self.species, self.reference, str(self.release))

    def make_release_cache_dir(self):
        """Create the release cache directory, if it doesn't already exist."""
        Path(self.release_cache_dir).mkdir(exist_ok=True, parents=True)

    def download_all(self):
        """Download required data and cache it."""
        self.make_release_cache_dir()

        if not os.path.exists(self.local_fasta_filepath):
            self.download_fasta()

        if not os.path.exists(self.local_gtf_cache_filepath):
            if not os.path.exists(self.local_gtf_filepath):
                self.download_gtf()
            self.cache_gtf()

    # ---------------------------------------------------------------------------------------------
    # FASTA file
    # ---------------------------------------------------------------------------------------------
    @property
    def remote_fasta_subdir(self) -> str:
        """Path to the directory on the FTP server that contains the FASTA file."""
        return FASTA_SUBDIR_TEMPLATE.format(release=self.release, species=self.species)

    @property
    def remote_fasta_filename(self) -> str:
        """Name of the FASTA file on the FTP server."""
        if self.release <= 75:
            return FASTA_FILENAME_TEMPLATE_DNA_OLD.format(
                species=self.species.capitalize(), reference=self.reference, release=self.release
            )
        else:
            return FASTA_FILENAME_TEMPLATE_DNA_NEW.format(
                species=self.species.capitalize(), reference=self.reference
            )

    @property
    def local_fasta_filepath(self) -> str:
        """Local name of the FASTA file."""
        return os.path.join(self.release_cache_dir, self.remote_fasta_filename)

    def download_fasta(self):
        """Download the FASTA file from the FTP server."""
        self._ftp_download(
            ENSEMBL_FTP_SERVER,
            self.remote_fasta_subdir,
            self.remote_fasta_filename,
            self.local_fasta_filepath,
        )

    # ---------------------------------------------------------------------------------------------
    # GTF file
    # ---------------------------------------------------------------------------------------------
    @property
    def remote_gtf_subdir(self) -> str:
        """Path to the directory on the FTP server that contains the GTF file."""
        return GTF_SUBDIR_TEMPLATE.format(release=self.release, species=self.species)

    @property
    def remote_gtf_filename(self) -> str:
        """Name of the GTF file on the FTP server."""
        return GTF_FILENAME_TEMPLATE.format(
            species=self.species.capitalize(),
            reference=self.reference,
            release=self.release,
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

    def cache_gtf(self):
        """Convert the Ensembl GTF to a pandas DataFrame and cache it."""
        logger.debug(f"Converting {self.local_gtf_filepath} to {self.local_gtf_cache_filepath}")
        read_gtf(self.local_gtf_filepath).to_pickle(self.local_gtf_cache_filepath)

    def load_gtf(self) -> pandas.DataFrame:
        """Return the cached Ensembl GTF data."""
        return pandas.read_pickle(self.local_gtf_cache_filepath)

    # ---------------------------------------------------------------------------------------------
    # generic functions
    # ---------------------------------------------------------------------------------------------
    def _ftp_download(self, server: str, subdir: str, remote_file: str, local_file: str):
        """Download a file from an FTP server."""
        try:
            ftp = FTP(server)
            logger.debug(f"Connecting to {server}")
            ftp.login()
            url = "ftp://" + server + subdir + "/" + remote_file
            logger.debug(f"Downloading {url} to {local_file}")
            ftp.cwd(subdir)
            self.make_release_cache_dir()
            with open(local_file, "wb") as fh:
                ftp.retrbinary(f"RETR {remote_file}", fh.write)
            logger.debug("Download successful")
            ftp.quit()
        except Exception as exc:
            logger.exception(f"Download failed: {exc}")


if __name__ == "__main__":
    DataCache(
        "homo_sapiens", "GRCh38", 100, "/home/matt/Downloads/ensembl_map_data/"
    ).download_all()  # DEBUG
    DataCache(
        "homo_sapiens", "GRCh37", 69, "/home/matt/Downloads/ensembl_map_data/"
    ).download_all()  # DEBUG
