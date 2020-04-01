import logging
import os

from pyensembl import EnsemblRelease

from .util import singleton


DEFAULT_RELEASE = 99
DEFAULT_SPECIES = "homo_sapiens"


def set_cache_dir(path):
    """Set the root directory where cache files are stored."""
    Ensembl()._set_cache_dir(path)


def set_ensembl_release(release, species=None, cache_dir=None, download_if_missing=False):
    """Set the Ensembl release to use."""
    return Ensembl().load(release, species, cache_dir, download_if_missing)


def release():
    """Return the Ensembl release being used."""
    return Ensembl().data.release


def species():
    """Return the species being used."""
    return Ensembl().data.species.latin_name


@singleton
class Ensembl:
    """Container for a `pyensembl.EnsemblRelease` object.
    
    Class handles the required cache files. Cache files are not loaded until needed (as
    opposed to loading them when the package is imported). If the cache files do not 
    exist locally, they will be downloaded and indexed automatically.

    Attributes:
        cache_dir (str): root directory where the `pyensembl` cache files are located
        data (`EnsemblRelease`): object for accessing annotation data and queries
    """

    def __init__(self):
        self._data = None
        self.cache_dir = None

    @property
    def data(self):
        if self._data is None:
            logging.debug("Loading default Ensembl release")
            self.load()
        return self._data

    def load(
        self,
        release=DEFAULT_RELEASE,
        species=DEFAULT_SPECIES,
        cache_dir=None,
        download_if_missing=False,
    ):
        if not release or not species:
            raise ValueError("A release number or species must be given")

        if cache_dir:
            self._set_cache_dir(cache_dir)
        self._set_ensembl_release(release, species)
        self._download_cache(download_if_missing)

    def _set_cache_dir(self, path):
        """Set the root directory where cache files are stored."""
        self.cache_dir = path
        os.environ["PYENSEMBL_CACHE_DIR"] = self.cache_dir
        logging.debug(f"Set PYENSEMBL_CACHE_DIR to {self.cache_dir}")

    def _set_ensembl_release(self, release, species):
        """Set the Ensembl release."""
        self._data = EnsemblRelease(release=release, species=species)
        logging.debug(f"Using Ensembl {release} (species={species})")

    def _download_cache(self, download_if_missing=True):
        """Download any missing cache files."""
        if self._data.required_local_files_exist():
            logging.info(f"Using cache files in {self._data.download_cache.cache_directory_path}")
        elif download_if_missing:
            logging.warning("Downloading required cache files (this may take a while)...")
            self._data.download()
            self._data.index()
        else:
            raise FileNotFoundError("Missing required cache files")
