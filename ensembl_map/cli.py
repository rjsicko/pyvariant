import argparse

from .cache import EnsemblCache, get_cache_dir


def parse_args() -> argparse.Namespace:
    """Parse commmand line arguments."""
    parser = argparse.ArgumentParser()
    subparser = parser.add_subparsers(dest="subcommand")

    # 'install' subcommand ------------------------------------------------------------------------
    install = subparser.add_parser("install", add_help=False)

    install_required = install.add_argument_group("required inputs")
    install_required.add_argument(
        "-r", "--release", type=int, required=True, help="Ensembl release number"
    )
    install_required.add_argument("-s", "--species", required=True, help="species name")

    install_optional = install.add_argument_group("optional inputs")
    install_optional.add_argument(
        "-c", "--cache", default=get_cache_dir(), help="cache directory (default = %(default)s)"
    )
    install_optional.add_argument(
        "-k",
        "--keep",
        dest="clean",
        action="store_false",
        help="do not delete original annotation files after caching",
    )
    install_optional.add_argument(
        "--recache", action="store_true", help="rebuild cache files (default = false)"
    )
    install_optional.add_argument(
        "--redownload", action="store_true", help="redownload annotation files (default = false)"
    )
    install_optional.add_argument(
        "-h", "--help", action="help", help="show this help message and exit"
    )

    return parser.parse_args()


def main():
    """Entrypoint for `ensembl_map`."""
    args = parse_args()

    if args.subcommand == "install":
        cache = EnsemblCache(species=args.species, release=args.release, cache_dir=args.cache)
        cache.install(clean=args.clean, recache=args.recache, redownload=args.redownload)
