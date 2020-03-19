from .cache import Cache


def set_ensembl_release(release=None, species=None):
    """Set the Ensembl release to use."""
    args = []
    if release is not None:
        args.append(release)
    if species is not None:
        args.append(species)
    if not release or species:
        raise ValueError("A release number of species must be given")

    return Cache.set_cache(release, species)
