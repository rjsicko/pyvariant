import csv
from typing import Dict, Optional

from logzero import logger


class ContigNames:
    """Helper class provides functions for dealing with chromsome/contig IDs."""

    _data: Dict[str, str] = {}
    _loaded = False

    @classmethod
    def load(cls, path: str = ""):
        """Load one-to-one mappings of known aliases to Ensembl contig IDs.

        Expects a tab-seperated file like:
            chr1    1
            chr2    2
            chr3    3
            ...

        Args:
            path (str): path to the file

        Example:
            >>> ContigNames.load()
            >>> ContigNames.normalize("chr12")
            '12'
        """
        cls._data = {}
        cls._loaded = False

        if path:
            logger.debug(f"Loading contig aliases from '{path}'")
            with open(path, "r") as fh:
                reader = csv.DictReader(fh, delimiter="\t", fieldnames=["alias", "name"])
                for row in reader:
                    cls._data[row["alias"]] = row["name"]
        else:
            logger.warning("No contig alias file")

        cls._loaded = True

    @classmethod
    def normalize(cls, key: str) -> Optional[str]:
        """Try to convert the given contig ID to the one used in the given release."""
        if not cls._loaded:
            cls.load()

        if key in cls._data:
            return cls._data[key]
        elif key in cls._data.values():
            return key
        else:
            return None
