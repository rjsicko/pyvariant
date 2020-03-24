from .cache import Ensembl


def set_cache_dir(path):
    """Set the root directory where cache files are stored."""
    Ensembl().set_cache_dir(path)


def set_ensembl_release(release, species, cache_dir=None, download_if_missing=False):
    """Set the Ensembl release to use."""
    return Ensembl().load(release, species, cache_dir, download_if_missing)
