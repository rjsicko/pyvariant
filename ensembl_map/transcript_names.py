import csv
from typing import Dict, List, Set

from logzero import logger

from .utils import iter_normalized_ids, skip_header


class TranscriptNames:
    """Helper class provides functions for dealing with transcript symbols."""

    _ens_to_ens: Dict[str, Set[str]] = {}
    _ens_to_refseq: Dict[str, Set[str]] = {}
    _refseq_to_ens: Dict[str, Set[str]] = {}
    _refseq_to_refseq: Dict[str, Set[str]] = {}
    _loaded = False

    @classmethod
    def load(cls, path: str = ""):
        """Load the annotations needed by the class."""
        cls._ens_to_ens = {}
        cls._ens_to_refseq = {}
        cls._refseq_to_ens = {}
        cls._refseq_to_refseq = {}
        cls._loaded = False

        if path:
            logger.debug(f"Loading transcript IDs from '{path}'")
            with open(path, "r") as fh:
                reader = csv.DictReader(skip_header(fh), delimiter="\t")
                for row in reader:
                    refseq_id_list = _get_refseq_ids(row)
                    ensembl_id_list = _get_ensembl_ids(row)
                    for refseq_id in refseq_id_list:
                        refseq_id_lower = refseq_id.lower()
                        cls._refseq_to_refseq.setdefault(refseq_id_lower, set())
                        cls._refseq_to_refseq[refseq_id_lower].add(refseq_id)
                        for ensembl_id in ensembl_id_list:
                            cls._refseq_to_ens.setdefault(refseq_id_lower, set())
                            cls._refseq_to_ens[refseq_id_lower].add(ensembl_id)
                    for ensembl_id in ensembl_id_list:
                        ensembl_id_lower = ensembl_id.lower()
                        cls._ens_to_ens.setdefault(ensembl_id_lower, set())
                        cls._ens_to_ens[ensembl_id_lower].add(ensembl_id)
                        for refseq_id in refseq_id_list:
                            cls._ens_to_refseq.setdefault(ensembl_id_lower, set())
                            cls._ens_to_refseq[ensembl_id_lower].add(refseq_id)
        else:
            logger.warning("No RefSeq-to-Ensembl mapping file")

        cls._loaded = True

    @classmethod
    def to_ensembl(cls, key: str) -> List[str]:
        """Map a transcript symbol to the corresponding Ensembl transcript ID."""
        if not cls._loaded:
            cls.load()

        for k in iter_normalized_ids(key):
            if k in cls._refseq_to_ens:
                return sorted(cls._refseq_to_ens[k])
            elif k in cls._ens_to_ens:
                return sorted(cls._ens_to_ens[k])
        else:
            return []

    @classmethod
    def to_refseq(cls, key: str) -> List[str]:
        """Map a transcript symbol to the corresponding RefSeq transcript ID."""
        if not cls._loaded:
            cls.load()

        for k in iter_normalized_ids(key):
            if k in cls._ens_to_refseq:
                return sorted(cls._ens_to_refseq[k])
            elif k in cls._refseq_to_refseq:
                return sorted(cls._refseq_to_refseq[k])
        else:
            return []


def _get_ensembl_ids(row: Dict) -> List[str]:
    """Select the value corresponding to the Ensembl transcript ID."""
    columns = [
        "Transcript stable ID",  # GRCh38
        "Ensembl Transcript ID",  # GRCh37
    ]
    return _get_non_empty(row, columns)


def _get_refseq_ids(row: Dict) -> List[str]:
    """Select the value corresponding to the RefSeq ID."""
    columns = [
        "RefSeq mRNA ID",  # GRCh38
        "RefSeq ncRNA ID",  # GRCh38
        "RefSeq mRNA",  # GRCh37
        "RefSeq ncRNA",  # GRCh37
    ]
    return _get_non_empty(row, columns)


def _get_non_empty(row: Dict, columns: List[str]) -> List[str]:
    """Return non-empty values from a list of columns."""
    assert any(
        [key in row for key in columns]
    ), f"Expected at least one of the following columns: {columns}"

    values = []
    for key in columns:
        if key in row and row[key]:
            values.append(row[key])

    return values
