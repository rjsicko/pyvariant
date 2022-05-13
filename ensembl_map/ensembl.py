import os
from typing import Optional

from logzero import logger
from pyensembl import MAX_ENSEMBL_RELEASE, EnsemblRelease
from pyensembl.download_cache import DownloadCache
from pyensembl.ensembl_url_templates import ENSEMBL_FTP_SERVER, FASTA_SUBDIR_TEMPLATE
from pyfaidx import Fasta

from .utils import bgzip, is_bgzipped, reverse_complement

# the file name format changed after Ensembl release 75
OLD_FASTA_FILENAME_TEMPLATE_DNA = "{species}.{reference}.{release}.dna.toplevel.fa.gz"
NEW_FASTA_FILENAME_TEMPLATE_DNA = "{species}.{reference}.dna.toplevel.fa.gz"


# --------------------------------------
#  EnsemblRelease functions
# --------------------------------------
def load_ensembl(
    release: int = MAX_ENSEMBL_RELEASE,
    species: str = "homo_sapiens",
    cache_dir: Optional[str] = None,
    download_if_missing: bool = False,
) -> EnsemblRelease:
    """Load the specified Ensembl release and species.

    Examples:
        >>> load_ensembl(69)
        EnsemblRelease(release=69, species='homo_sapiens')

        >>> load_ensembl(100)
        EnsemblRelease(release=100, species='homo_sapiens')
    """
    logger.debug(f"Loading Ensembl release {release} (species: {species})")

    # pyensembl cache directory is read from an enviornmental variable
    if cache_dir:
        os.environ["PYENSEMBL_CACHE_DIR"] = cache_dir

    # set the Ensembl release
    ensembl = EnsemblRelease(release, species)

    # check if the required cache files exist
    if ensembl.required_local_files_exist():
        logger.debug(f"Using Ensembl files in {ensembl.download_cache.cache_directory_path}")
    elif download_if_missing:
        logger.info("Downloading required Ensembl files (this may take some time)...")
        ensembl.download()
        ensembl.index()
    else:
        raise FileNotFoundError("Missing required Ensembl files")

    return ensembl


# --------------------------------------
#  Genome FASTA functions
# --------------------------------------
class Genome:
    """Class for interacting with an Ensembl reference genome."""

    def __init__(
        self,
        release: int,
        species: str,
        reference: str,
        cache_dir: str,
        download_if_missing: bool = False,
    ):
        self.release = release
        self.species = species
        self.reference = reference
        self.cache_dir = cache_dir
        self.fasta_file = self._download_if_missing(download_if_missing)
        self.index = self._get_index()
        self.fasta = self._load_fasta()

    def sequence(
        self, contig: str, start: int, end: Optional[int] = None, strand: str = "+"
    ) -> str:
        """Return the nucleotide sequence at the given contig coordinates."""
        contiglen = len(self.fasta[contig])
        end = end if end else start
        start = start - 1  # 0-based coordiantes
        if not (0 <= start < contiglen):
            raise ValueError(f"Start must be from 1 to {contiglen} ({start})")
        if not (1 < end <= contiglen):
            raise ValueError(f"End must be from 2 to {contiglen + 1} ({end})")
        if end <= start:
            raise ValueError(f"End must be > start ({end} <= {start})")

        seq = self.fasta[contig][start:end]
        if strand == "-":
            seq = reverse_complement(seq)

        return seq

    def _download_if_missing(self, download_if_missing: bool) -> str:
        """Download the contig sequences from Ensembl, if they don't already exist."""
        logger.debug(f"Loading Ensembl genome {self.release} (species: {self.species})")

        cache = DownloadCache(
            reference_name=self.reference,
            annotation_name="ensembl",
            annotation_version=self.release,
            cache_directory_path=self.cache_dir,
            decompress_on_download=False,
        )

        fasta_file_url = self._get_fasta_url()
        fasta_file = cache.download_or_copy_if_necessary(
            fasta_file_url,
            download_if_missing=download_if_missing,
        )

        # pyfaidx only supports compressed FASTA in BGZF format. Ensembl FASTA comes in GZ format,
        # so we need to compress the data with bgzip.
        if not is_bgzipped(fasta_file):
            logger.info(f"Compressing {fasta_file} with bgzip (this may take some time)...")
            fasta_file = bgzip(fasta_file)

        return fasta_file

    def _get_fasta_url(self) -> str:
        """Get the URL to the Ensembl genome FASTA file(s)."""
        if self.release <= 75:
            fasta_file_name = OLD_FASTA_FILENAME_TEMPLATE_DNA.format(
                release=self.release,
                reference=self.reference,
                species=self.species.capitalize(),
            )
        else:
            fasta_file_name = NEW_FASTA_FILENAME_TEMPLATE_DNA.format(
                reference=self.reference,
                species=self.species.capitalize(),
            )

        fasta_file_url = (
            ENSEMBL_FTP_SERVER
            + FASTA_SUBDIR_TEMPLATE
            % {
                "release": self.release,
                "species": self.species,
                "type": "dna",
            }
            + fasta_file_name
        )

        return fasta_file_url

    def _get_index(self) -> str:
        """Return the path to the FASTA index file."""
        if self.fasta_file.endswith(".gz"):
            return self.fasta_file + ".tbi"
        else:
            return self.fasta_file + ".fai"

    def _load_fasta(self) -> Fasta:
        """Return a `pyfaidx.Fasta` object for the Ensembl genome."""
        if not os.path.exists(self.index):
            # index is generated when the `Fasta` object is initialized
            logger.warning(f"Indexing {self.fasta_file} (this may take some time)...")

        return Fasta(self.fasta_file, as_raw=True, sequence_always_upper=True)
