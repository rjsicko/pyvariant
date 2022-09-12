from functools import lru_cache
from typing import Optional, Union

from .constants import CDNA, CONTIG, GENE, PROTEIN, TRANSCRIPT
from .features import Cdna, Contig, Gene, Protein, Transcript
from .map import get_map_func


@lru_cache()
def cdna_sequence(
    cdna: str, start: int, end: Optional[int] = None, raise_error: bool = True
) -> str:
    """Return the sequence of a CDNA at the given position(s)."""
    return _get_feature(cdna, start, end, CDNA, raise_error).sequence


@lru_cache()
def contig_sequence(
    contig: str, start: int, end: Optional[int] = None, raise_error: bool = True
) -> str:
    """Return the sequence of a contig at the given position(s)."""
    return _get_feature(contig, start, end, CONTIG, raise_error).sequence


@lru_cache()
def gene_sequence(
    gene: str, start: int, end: Optional[int] = None, raise_error: bool = True
) -> str:
    """Return the sequence of a gene at the given position(s)."""
    return _get_feature(gene, start, end, GENE, raise_error).sequence


@lru_cache()
def protein_sequence(
    protein: str, start: int, end: Optional[int] = None, raise_error: bool = True
) -> str:
    """Return the sequence of a protein at the given position(s)."""
    return _get_feature(protein, start, end, PROTEIN, raise_error).sequence


@lru_cache()
def transcript_sequence(
    transcript: str, start: int, end: Optional[int] = None, raise_error: bool = True
) -> str:
    """Return the sequence of a transcript at the given position(s)."""
    return _get_feature(transcript, start, end, TRANSCRIPT, raise_error).sequence


def _get_feature(
    feature: str, start: int, end: Optional[int], feature_type: str, raise_error: bool
) -> Union[Cdna, Contig, Gene, Protein, Transcript]:
    func = get_map_func(feature_type, feature_type)
    result = func(feature, start, end, raise_error=raise_error, feature_type=feature_type)
    if len(result) != 1:
        raise ValueError(f"Unable to uniquely map {feature_type} {feature, start, end}: {result}")
    else:
        return result[0]
