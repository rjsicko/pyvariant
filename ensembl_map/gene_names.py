import csv
from typing import Dict, Optional

from logzero import logger

from .ensembl import release
from .utils import is_ensembl_id, iter_normalized_ids


class GeneNames:
    """Helper class provides functions for dealing with gene symbols."""

    _name_to_id: Dict[str, Dict[str, str]] = {}
    _id_to_name: Dict[str, Dict[str, str]] = {}
    _loaded = False

    @classmethod
    def load(cls, path: str = ""):
        """Load the annotations needed by the class.

        The expected annotation file has one column labelled 'gene_id' that lists the Ensembl gene
        ID, then additional columns where the label is the Ensembl release version and the value is
        the HUGO gene symbol. Example:
            gene_id         69      100
            ENSG00000000419 DPM1    DPM1
            ENSG00000002726 ABP1    AOC1
            ENSG00000003509 C2orf56 NDUFAF7
        """
        cls._name_to_id = {}
        cls._id_to_name = {}
        cls._loaded = False

        if path:
            with open(path, "r") as fh:
                reader = csv.DictReader(fh, dialect="excel-tab")
                for row in reader:
                    gene_id = row["gene_id"]
                    name_by_release = {k: row[k] for k in row if k != "gene_id" and row[k]}
                    id_by_release = {k: gene_id for k in row if k != "gene_id" and row[k]}
                    cls._id_to_name[gene_id] = name_by_release
                    for _, gene_name in name_by_release.items():
                        cls._name_to_id.setdefault(gene_name, {}).update(id_by_release)
        else:
            logger.warning("No gene name change file")

        cls._loaded = True

    @classmethod
    def normalize(cls, key: str) -> Optional[str]:
        """Convert a gene name or ID to the equivalent representation in the given Ensembl release."""
        if not cls._loaded:
            cls.load()

        if is_ensembl_id(key):
            return cls.to_id(key)
        else:
            return cls.to_name(key)

    @classmethod
    def to_id(cls, key: str) -> Optional[str]:
        """Convert a gene symbol to the Ensembl gene ID used in the given Ensembl release."""
        if not cls._loaded:
            cls.load()

        for k in iter_normalized_ids(key):
            if k in cls._id_to_name:
                return key
            elif k in cls._name_to_id:
                id_by_release = cls._name_to_id[k]
                return id_by_release.get(str(release()), None)
        else:
            return None

    @classmethod
    def to_name(cls, key: str) -> Optional[str]:
        """Convert a gene symbol to the HUGO gene symbol used in the given Ensembl release."""
        if not cls._loaded:
            cls.load()

        if (gene_id := cls.to_id(key)) in cls._id_to_name:
            name_by_release = cls._id_to_name[gene_id]
            return name_by_release.get(str(release()), None)
        else:
            return None
