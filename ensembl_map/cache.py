import logging
import os

from pyensembl import EnsemblRelease


DEFAULT_RELEASE = 99
DEFAULT_SPECIES = "homo_sapiens"


class Cache:
    """
    Class handles the required cache files. Cache files are not loaded until needed (as
    opposed to loading them when the package is imported). If the cache files do not 
    exist locally, they will be downloaded and indexed automatically.
    """

    _cache = None

    @classmethod
    def get_cache(cls):
        if not cls._cache:
            cls._cache = cls.set_cache()
        return cls._cache

    @classmethod
    def set_cache(
        cls, release=DEFAULT_RELEASE, species=DEFAULT_SPECIES, download_if_missing=False
    ):
        cls._cache = EnsemblRelease(release, species)
        if cls._cache.required_local_files_exist():
            logging.info(
                f"Using cache files in {cls._cache.download_cache.cache_directory_path}"
            )
        elif download_if_missing:
            logging.warning(
                "Downloading required cache files (this may take a while)..."
            )
            cls._cache.download()
            cls._cache.index()
        else:
            raise FileNotFoundError("Missing required cache files")

        return cls._cache

    @classmethod
    def set_cache_dir(cls, path):
        """Set the root directory where cache files are stored."""
        os.environ["PYENSEMBL_CACHE_DIR"] = path

