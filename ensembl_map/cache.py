import logging

from pyensembl import EnsemblRelease


DEFAULT_RELEASE = 99
DEFAULT_SPECIES = "homo_sapiens"


class Cache:
    """
    Class handles downloading and indexing the required cache files, if they do not 
    exist locally. The intention is that cache files are not loaded until needed (as 
    opposed to loading them on import), to give the user the chance to specify which 
    release to use before loading.
    """

    _cache = None

    @classmethod
    def get_cache(cls):
        if not cls._cache:
            cls._cache = cls.set_cache()
        return cls._cache

    @classmethod
    def set_cache(cls, release=DEFAULT_RELEASE, species=DEFAULT_SPECIES):
        cls._cache = EnsemblRelease(release, species)
        if cls._cache.required_local_files_exist():
            logging.info(
                f"Using cache files in {cls._cache.download_cache.cache_directory_path}"
            )
        else:
            logging.warning(
                "Downloading required cache files (this may take a while)..."
            )
            cls._cache.download()
            cls._cache.index()

        return cls._cache


def set_ensembl_release(release, species=DEFAULT_SPECIES):
    """Set the Ensembl release to use."""
    args = []
    if release is not None:
        args.append(release)
    if species is not None:
        args.append(species)

    return Cache.set_cache(release, species)

