from .cache import Cache


def set_cache_dir(path):
    """Set the root directory where cache files are stored."""
    Cache.set_cache_dir(path)


def set_ensembl_release(release=None, species=None, download_if_missing=False):
    """Set the Ensembl release to use."""
    args = []
    if release is not None:
        args.append(release)
    if species is not None:
        args.append(species)
    if not release and not species:
        raise ValueError("A release number of species must be given")

    return Cache.set_cache(release, species, download_if_missing)

