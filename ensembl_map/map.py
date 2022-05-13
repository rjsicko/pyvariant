from functools import lru_cache
from math import floor
from typing import Any, Callable, List, Optional, Set, Tuple

import pyensembl
from logzero import logger

from .constants import CDS, CONTIG, EXON, GENE, PROTEIN, STRAND_SYMBOLS, TRANSCRIPT
from .core import instance as CM
from .exceptions import CdsOutOfRange, ExonOutOfRange, GeneOutOfRange, TranscriptOutOfRange
from .features import Contig, get_load_function
from .utils import assert_valid_position


# --------------------------------------
#  General mapping functions
# --------------------------------------
@lru_cache()
def cds_to_cds(
    feature: str,
    start: int,
    end: Optional[int] = None,
    raise_error: bool = True,
    best_only: bool = False,
    feature_type: Optional[str] = None,
) -> List:
    """Map CDS coordinates."""
    result = set()

    for feature_id in CM().transcript_ids(feature, feature_type, best_only):
        result |= _map(feature_id, start, end, CDS, CDS, best_only=best_only)

    return _validate_result(result, "CDS", "CDS", (feature, start, end), raise_error)


@lru_cache()
def cds_to_contig(
    feature: str,
    start: int,
    end: Optional[int] = None,
    raise_error: bool = True,
    best_only: bool = False,
    feature_type: Optional[str] = None,
) -> List:
    """Map CDS coordinates to contig coordinates."""
    result = set()

    for feature_id in CM().transcript_ids(feature, feature_type, best_only):
        for pos in cds_to_gene(
            feature_id,
            start,
            end,
            raise_error=False,
            best_only=best_only,
            feature_type=TRANSCRIPT,
        ):
            result |= _map(pos.gene_id, pos.start, pos.end, GENE, CONTIG, best_only=best_only)

    return _validate_result(result, "CDS", "contig", (feature, start, end), raise_error)


@lru_cache()
def cds_to_exon(
    feature: str,
    start: int,
    end: Optional[int] = None,
    raise_error: bool = True,
    best_only: bool = False,
    feature_type: Optional[str] = None,
) -> List:
    """Map CDS coordinates to exon coordinates."""
    result = set()

    for feature_id in CM().transcript_ids(feature, feature_type, best_only):
        for pos in cds_to_transcript(
            feature_id,
            start,
            end,
            raise_error=False,
            best_only=best_only,
            feature_type=TRANSCRIPT,
        ):
            for pos2 in transcript_to_exon(
                pos.transcript_id,
                pos.start,
                pos.end,
                raise_error=False,
                best_only=best_only,
                feature_type=TRANSCRIPT,
            ):
                if pos2.transcript_id == feature_id:
                    result.add(pos2)
    return _validate_result(result, "CDS", "exon", (feature, start, end), raise_error)


@lru_cache()
def cds_to_gene(
    feature: str,
    start: int,
    end: Optional[int] = None,
    raise_error: bool = True,
    best_only: bool = False,
    feature_type: Optional[str] = None,
) -> List:
    """Map CDS coordinates to gene coordinates."""
    result = set()

    for feature_id in CM().transcript_ids(feature, feature_type, best_only):
        for pos in cds_to_transcript(
            feature_id,
            start,
            end,
            raise_error=False,
            best_only=best_only,
            feature_type=TRANSCRIPT,
        ):
            for pos2 in transcript_to_gene(
                pos.transcript_id,
                pos.start,
                pos.end,
                raise_error=False,
                best_only=best_only,
                feature_type=TRANSCRIPT,
            ):
                result.add(pos2)

    return _validate_result(result, "CDS", "gene", (feature, start, end), raise_error)


@lru_cache()
def cds_to_protein(
    feature: str,
    start: int,
    end: Optional[int] = None,
    raise_error: bool = True,
    best_only: bool = False,
    feature_type: Optional[str] = None,
) -> List:
    """Map CDS coordinates to protein coordinates."""
    result = set()

    for feature_id in CM().transcript_ids(feature, feature_type, best_only):
        result |= _map(feature_id, start, end, CDS, PROTEIN, best_only=best_only)

    return _validate_result(result, "CDS", "protein", (feature, start, end), raise_error)


@lru_cache()
def cds_to_transcript(
    feature: str,
    start: int,
    end: Optional[int] = None,
    raise_error: bool = True,
    best_only: bool = False,
    feature_type: Optional[str] = None,
) -> List:
    """Map CDS coordinates to transcript coordinates."""
    result = set()

    for feature_id in CM().transcript_ids(feature, feature_type, best_only):
        result |= _map(feature_id, start, end, CDS, TRANSCRIPT, best_only=best_only)

    return _validate_result(result, "CDS", "transcript", (feature, start, end), raise_error)


@lru_cache()
def contig_to_cds(
    feature: str,
    start: int,
    end: Optional[int] = None,
    raise_error: bool = True,
    best_only: bool = False,
    feature_type: Optional[str] = None,
    strand: Optional[str] = None,
) -> List:
    """Map contig coordinates to CDS coordinates."""
    result = set()

    if strand not in STRAND_SYMBOLS:
        raise ValueError(f"strand must be one of {STRAND_SYMBOLS}")

    for feature_id in CM().contig_ids(feature, feature_type, best_only):
        for pos in contig_to_gene(
            feature_id,
            start,
            end,
            raise_error=False,
            best_only=best_only,
            feature_type=CONTIG,
            strand=strand,
        ):
            for pos2 in gene_to_cds(
                pos.gene_id,
                pos.start,
                pos.end,
                raise_error=False,
                best_only=best_only,
                feature_type=GENE,
            ):
                result.add(pos2)

    return _validate_result(result, "contig", "CDS", (feature, start, end), raise_error)


@lru_cache()
def contig_to_contig(
    feature: str,
    start: int,
    end: Optional[int] = None,
    raise_error: bool = True,
    best_only: bool = False,
    feature_type: Optional[str] = None,
    strand: Optional[str] = None,
) -> List:
    """Map contig coordinates."""
    result = set()

    strand = strand if strand is not None else "+"
    if strand not in STRAND_SYMBOLS:
        raise ValueError(f"strand must be one of {STRAND_SYMBOLS}")

    # A contig/chromosomal position may not always map to a gene, therefore don't try to map the given
    # coordinates. Instead just format the values as a Contig object and return. TODO: validate the
    # given coordinates against a reference fasta?
    for feature_id in CM().contig_ids(feature, feature_type, best_only):
        result.add(Contig(contig=feature_id, start=start, end=end or start, strand=strand))

    return _validate_result(result, "contig", "contig", (feature, start, end), raise_error)


@lru_cache()
def contig_to_exon(
    feature: str,
    start: int,
    end: Optional[int] = None,
    raise_error: bool = True,
    best_only: bool = False,
    feature_type: Optional[str] = None,
    strand: Optional[str] = None,
) -> List:
    """Map contig coordinates to exon coordinates."""
    result = set()

    if strand not in STRAND_SYMBOLS:
        raise ValueError(f"strand must be one of {STRAND_SYMBOLS}")

    for feature_id in CM().contig_ids(feature, feature_type, best_only):
        for pos in contig_to_transcript(
            feature_id,
            start,
            end,
            raise_error=False,
            best_only=best_only,
            feature_type=CONTIG,
            strand=strand,
        ):
            for pos2 in transcript_to_exon(
                pos.transcript_id,
                pos.start,
                pos.end,
                raise_error=False,
                best_only=best_only,
                feature_type=TRANSCRIPT,
            ):
                result.add(pos2)

    return _validate_result(result, "contig", "exon", (feature, start, end), raise_error)


@lru_cache()
def contig_to_gene(
    feature: str,
    start: int,
    end: Optional[int] = None,
    raise_error: bool = True,
    best_only: bool = False,
    feature_type: Optional[str] = None,
    strand: Optional[str] = None,
) -> List:
    """Map contig coordinates to gene coordinates."""
    result = set()

    end = end if end is not None else start
    if strand not in STRAND_SYMBOLS:
        raise ValueError(f"strand must be one of {STRAND_SYMBOLS}")

    for feature_id in CM().contig_ids(feature, feature_type, best_only):
        for gene in CM().ensembl.genes_at_locus(feature_id, start, end, strand):
            gstart = start if start >= gene.start else gene.start
            gend = end if end <= gene.end else gene.end
            result |= _map(gene.gene_id, gstart, gend, GENE, GENE, best_only=best_only)

    return _validate_result(result, "contig", "gene", (feature, start, end), raise_error)


@lru_cache()
def contig_to_protein(
    feature: str,
    start: int,
    end: Optional[int] = None,
    raise_error: bool = True,
    best_only: bool = False,
    feature_type: Optional[str] = None,
    strand: Optional[str] = None,
) -> List:
    """Map contig coordinates to protein coordinates."""
    result = set()

    if strand not in STRAND_SYMBOLS:
        raise ValueError(f"strand must be one of {STRAND_SYMBOLS}")

    for feature_id in CM().contig_ids(feature, feature_type, best_only):
        for pos in contig_to_gene(
            feature_id,
            start,
            end,
            raise_error=False,
            best_only=best_only,
            feature_type=CONTIG,
            strand=strand,
        ):
            for pos2 in gene_to_protein(
                pos.gene_id,
                pos.start,
                pos.end,
                raise_error=False,
                best_only=best_only,
                feature_type=GENE,
            ):
                result.add(pos2)

    return _validate_result(result, "contig", "protein", (feature, start, end), raise_error)


@lru_cache()
def contig_to_transcript(
    feature: str,
    start: int,
    end: Optional[int] = None,
    raise_error: bool = True,
    best_only: bool = False,
    feature_type: Optional[str] = None,
    strand: Optional[str] = None,
) -> List:
    """Map contig coordinates to transcript coordinates."""
    result = set()

    if strand not in STRAND_SYMBOLS:
        raise ValueError(f"strand must be one of {STRAND_SYMBOLS}")

    for feature_id in CM().contig_ids(feature, feature_type, best_only):
        for pos in contig_to_gene(
            feature_id,
            start,
            end,
            raise_error=False,
            best_only=best_only,
            feature_type=CONTIG,
            strand=strand,
        ):
            for pos2 in gene_to_transcript(
                pos.gene_id,
                pos.start,
                pos.end,
                raise_error=False,
                best_only=best_only,
                feature_type=GENE,
            ):
                result.add(pos2)

    return _validate_result(result, "contig", "transcript", (feature, start, end), raise_error)


@lru_cache()
def exon_to_cds(
    feature: str,
    position: Optional[int] = None,
    _: Optional[int] = None,
    raise_error: bool = True,
    best_only: bool = False,
    feature_type: Optional[str] = None,
) -> List:
    """Map an exon to CDS coordinates."""
    result = set()

    for feature_id in CM().exon_ids(feature, feature_type, best_only):
        for pos in exon_to_gene(
            feature_id, position, raise_error=False, best_only=best_only, feature_type=EXON
        ):
            result |= _map(pos.gene_id, pos.start, pos.end, GENE, CDS, best_only=best_only)

    return _validate_result(result, "exon", "CDS", (feature, position), raise_error)


@lru_cache()
def exon_to_contig(
    feature: str,
    position: Optional[int] = None,
    _: Optional[int] = None,
    raise_error: bool = True,
    best_only: bool = False,
    feature_type: Optional[str] = None,
) -> List:
    """Map exon coordinates to contig coordinates."""
    result = set()

    for feature_id in CM().exon_ids(feature, feature_type, best_only):
        for pos in exon_to_gene(
            feature_id, position, raise_error=False, best_only=best_only, feature_type=EXON
        ):
            result |= _map(pos.gene_id, pos.start, pos.end, GENE, CONTIG, best_only=best_only)

    return _validate_result(result, "exon", "contig", (feature, position), raise_error)


@lru_cache()
def exon_to_exon(
    feature: str,
    position: Optional[int] = None,
    _: Optional[int] = None,
    raise_error: bool = True,
    best_only: bool = False,
    feature_type: Optional[str] = None,
) -> List:
    """Map an exon."""
    result = set()

    for feature_id in CM().exon_ids(feature, feature_type, best_only):
        try:
            exon = CM().ensembl.exon_by_id(feature_id)
        except TypeError:
            continue
        else:
            for pos in _map(exon.exon_id, exon.start, exon.end, EXON, EXON, best_only=best_only):
                if not position or pos.index == position:
                    result.add(pos)

    return _validate_result(result, "exon", "exon", (feature, position), raise_error)


@lru_cache()
def exon_to_gene(
    feature: str,
    position: Optional[int] = None,
    _: Optional[int] = None,
    raise_error: bool = True,
    best_only: bool = False,
    feature_type: Optional[str] = None,
) -> List:
    """Map an exon to gene coordinates."""
    result = set()

    for feature_id in CM().exon_ids(feature, feature_type, best_only):
        for pos in exon_to_transcript(
            feature_id, position, raise_error=False, best_only=best_only, feature_type=EXON
        ):
            for pos2 in transcript_to_gene(
                pos.transcript_id,
                pos.start,
                pos.end,
                raise_error=False,
                best_only=best_only,
                feature_type=TRANSCRIPT,
            ):
                result.add(pos2)

    return _validate_result(result, "exon", "gene", (feature, position), raise_error)


@lru_cache()
def exon_to_protein(
    feature: str,
    position: Optional[int] = None,
    _: Optional[int] = None,
    raise_error: bool = True,
    best_only: bool = False,
    feature_type: Optional[str] = None,
) -> List:
    """Map an exon to protein coordinates."""
    result = set()

    for feature_id in CM().exon_ids(feature, feature_type, best_only):
        for pos in exon_to_cds(
            feature_id, position, raise_error=False, best_only=best_only, feature_type=EXON
        ):
            for pos2 in cds_to_protein(
                pos.transcript_id,
                pos.start,
                pos.end,
                raise_error=False,
                best_only=best_only,
                feature_type=TRANSCRIPT,
            ):
                result.add(pos2)

    return _validate_result(result, "exon", "protein", (feature, position), raise_error)


@lru_cache()
def exon_to_transcript(
    feature: str,
    position: Optional[int] = None,
    _: Optional[int] = None,
    raise_error: bool = True,
    best_only: bool = False,
    feature_type: Optional[str] = None,
) -> List:
    """Map an exon to transcript coordinates."""
    result = set()

    for feature_id in CM().exon_ids(feature, feature_type, best_only):
        try:
            exon = CM().ensembl.exon_by_id(feature_id)
        except TypeError:
            continue
        else:
            # Filter out transcripts that don't contain the query exon.
            valid_transcripts = _transcript_ids_with_exon(feature_id, position)
            for pos in gene_to_transcript(
                exon.gene_id,
                exon.start,
                exon.end,
                raise_error=False,
                best_only=best_only,
                feature_type=GENE,
            ):
                if pos.transcript_id in valid_transcripts:
                    result.add(pos)
                else:
                    logger.debug(f"{pos.transcript_id} does not contain exon {exon.exon_id}")

    return _validate_result(result, "exon", "transcript", (feature, position), raise_error)


@lru_cache()
def gene_to_cds(
    feature: str,
    start: int,
    end: Optional[int] = None,
    raise_error: bool = True,
    best_only: bool = False,
    feature_type: Optional[str] = None,
) -> List:
    """Map gene coordinates to CDS coordinates."""
    result = set()

    if feature_type == CONTIG or CM().is_contig(feature):
        return contig_to_cds((feature, start, end), raise_error, best_only, feature_type)

    for feature_id in CM().gene_ids(feature, feature_type):
        for pos in gene_to_transcript(
            feature_id, start, end, raise_error=False, best_only=best_only, feature_type=GENE
        ):
            for pos2 in transcript_to_cds(
                pos.transcript_id,
                pos.start,
                pos.end,
                raise_error=False,
                best_only=best_only,
                feature_type=TRANSCRIPT,
            ):
                result.add(pos2)

    return _validate_result(result, "gene", "CDS", (feature, start, end), raise_error)


@lru_cache()
def gene_to_contig(
    feature: str,
    start: int,
    end: Optional[int] = None,
    raise_error: bool = True,
    best_only: bool = False,
    feature_type: Optional[str] = None,
) -> List:
    """Map gene coordinates to contig coordinates."""
    result = set()

    if feature_type == CONTIG or CM().is_contig(feature):
        return contig_to_contig((feature, start, end), raise_error, best_only, feature_type)

    for feature_id in CM().gene_ids(feature, feature_type):
        result |= _map(feature_id, start, end, GENE, CONTIG, best_only=best_only)

    return _validate_result(result, "gene", "contig", (feature, start, end), raise_error)


@lru_cache()
def gene_to_exon(
    feature: str,
    start: int,
    end: Optional[int] = None,
    raise_error: bool = True,
    best_only: bool = False,
    feature_type: Optional[str] = None,
) -> List:
    """Map gene coordinates to exon coordinates."""
    result = set()

    if feature_type == CONTIG or CM().is_contig(feature):
        return contig_to_exon((feature, start, end), raise_error, best_only, feature_type)

    for feature_id in CM().gene_ids(feature, feature_type):
        result |= _map(feature_id, start, end, GENE, EXON, best_only=best_only)

    return _validate_result(result, "gene", "exon", (feature, start, end), raise_error)


@lru_cache()
def gene_to_gene(
    feature: str,
    start: int,
    end: Optional[int] = None,
    raise_error: bool = True,
    best_only: bool = False,
    feature_type: Optional[str] = None,
) -> List:
    """Map gene coordinates."""
    result = set()

    if feature_type == CONTIG or CM().is_contig(feature):
        return contig_to_gene((feature, start, end), raise_error, best_only, feature_type)

    for feature_id in CM().gene_ids(feature, feature_type):
        result |= _map(feature_id, start, end, GENE, GENE, best_only=best_only)

    return _validate_result(result, "gene", "gene", (feature, start, end), raise_error)


@lru_cache()
def gene_to_protein(
    feature: str,
    start: int,
    end: Optional[int] = None,
    raise_error: bool = True,
    best_only: bool = False,
    feature_type: Optional[str] = None,
) -> List:
    """Map gene coordinates to protein coordinates."""
    result = set()

    if feature_type == CONTIG or CM().is_contig(feature):
        return contig_to_protein((feature, start, end), raise_error, best_only, feature_type)

    for feature_id in CM().gene_ids(feature, feature_type):
        for pos in gene_to_cds(
            feature_id, start, end, raise_error=False, best_only=best_only, feature_type=GENE
        ):
            for pos2 in cds_to_protein(
                pos.transcript_id,
                pos.start,
                pos.end,
                raise_error=False,
                best_only=best_only,
                feature_type=TRANSCRIPT,
            ):
                result.add(pos2)

    return _validate_result(result, "gene", "protein", (feature, start, end), raise_error)


@lru_cache()
def gene_to_transcript(
    feature: str,
    start: int,
    end: Optional[int] = None,
    raise_error: bool = True,
    best_only: bool = False,
    feature_type: Optional[str] = None,
) -> List:
    """Map gene coordinates to transcript coordinates."""
    result = set()

    if feature_type == CONTIG or CM().is_contig(feature):
        return contig_to_transcript((feature, start, end), raise_error, best_only, feature_type)

    for feature_id in CM().gene_ids(feature, feature_type):
        result |= _map(feature_id, start, end, GENE, TRANSCRIPT, best_only=best_only)

    return _validate_result(result, "gene", "transcript", (feature, start, end), raise_error)


@lru_cache()
def protein_to_cds(
    feature: str,
    start: int,
    end: Optional[int] = None,
    raise_error: bool = True,
    best_only: bool = False,
    feature_type: Optional[str] = None,
) -> List:
    """Map protein coordinates to CDS coordinates."""
    result = set()

    for feature_id in CM().protein_ids(feature, feature_type, best_only):
        result |= _map(feature_id, start, end, PROTEIN, CDS, best_only=best_only)

    return _validate_result(result, "protein", "CDS", (feature, start, end), raise_error)


@lru_cache()
def protein_to_contig(
    feature: str,
    start: int,
    end: Optional[int] = None,
    raise_error: bool = True,
    best_only: bool = False,
    feature_type: Optional[str] = None,
) -> List:
    """Map protein coordinates to contig coordinates."""
    result = set()

    for feature_id in CM().protein_ids(feature, feature_type, best_only):
        for pos in protein_to_gene(
            feature_id, start, end, raise_error=False, best_only=best_only, feature_type=PROTEIN
        ):
            result |= _map(pos.gene_id, pos.start, pos.end, GENE, CONTIG, best_only=best_only)

    return _validate_result(result, "protein", "contig", (feature, start, end), raise_error)


@lru_cache()
def protein_to_exon(
    feature: str,
    start: int,
    end: Optional[int] = None,
    raise_error: bool = True,
    best_only: bool = False,
    feature_type: Optional[str] = None,
) -> List:
    """Map protein coordinates to exon coordinates."""
    result = set()

    for feature_id in CM().protein_ids(feature, feature_type, best_only):
        for pos in protein_to_transcript(
            feature_id, start, end, raise_error=False, best_only=best_only, feature_type=PROTEIN
        ):
            for pos2 in transcript_to_exon(
                pos.transcript_id,
                pos.start,
                pos.end,
                raise_error=False,
                best_only=best_only,
                feature_type=TRANSCRIPT,
            ):
                result.add(pos2)

    return _validate_result(result, "protein", "exon", (feature, start, end), raise_error)


@lru_cache()
def protein_to_gene(
    feature: str,
    start: int,
    end: Optional[int] = None,
    raise_error: bool = True,
    best_only: bool = False,
    feature_type: Optional[str] = None,
) -> List:
    """Map protein coordinates to gene coordinates."""
    result = set()

    for feature_id in CM().protein_ids(feature, feature_type, best_only):
        for pos in protein_to_transcript(
            feature_id, start, end, raise_error=False, best_only=best_only, feature_type=PROTEIN
        ):
            for pos2 in transcript_to_gene(
                pos.transcript_id,
                pos.start,
                pos.end,
                raise_error=False,
                best_only=best_only,
                feature_type=TRANSCRIPT,
            ):
                result.add(pos2)

    return _validate_result(result, "protein", "gene", (feature, start, end), raise_error)


@lru_cache()
def protein_to_protein(
    feature: str,
    start: int,
    end: Optional[int] = None,
    raise_error: bool = True,
    best_only: bool = False,
    feature_type: Optional[str] = None,
) -> List:
    """Map protein coordinates."""
    result = set()

    for feature_id in CM().protein_ids(feature, feature_type, best_only):
        result |= _map(feature_id, start, end, PROTEIN, PROTEIN, best_only=best_only)

    return _validate_result(result, "protein", "protein", (feature, start, end), raise_error)


@lru_cache()
def protein_to_transcript(
    feature: str,
    start: int,
    end: Optional[int] = None,
    raise_error: bool = True,
    best_only: bool = False,
    feature_type: Optional[str] = None,
) -> List:
    """Map protein coordinates to transcript coordinates."""
    result = set()

    for feature_id in CM().protein_ids(feature, feature_type, best_only):
        for pos in protein_to_cds(
            feature_id, start, end, raise_error=False, best_only=best_only, feature_type=PROTEIN
        ):
            for pos2 in cds_to_transcript(
                pos.transcript_id,
                pos.start,
                pos.end,
                raise_error=False,
                best_only=best_only,
                feature_type=TRANSCRIPT,
            ):
                result.add(pos2)

    return _validate_result(result, "protein", "transcript", (feature, start, end), raise_error)


@lru_cache()
def transcript_to_cds(
    feature: str,
    start: int,
    end: Optional[int] = None,
    raise_error: bool = True,
    best_only: bool = False,
    feature_type: Optional[str] = None,
) -> List:
    """Map transcript coordinates to CDS coordinates."""
    result = set()

    for feature_id in CM().transcript_ids(feature, feature_type, best_only):
        result |= _map(feature_id, start, end, TRANSCRIPT, CDS, best_only=best_only)

    return _validate_result(result, "transcript", "CDS", (feature, start, end), raise_error)


@lru_cache()
def transcript_to_contig(
    feature: str,
    start: int,
    end: Optional[int] = None,
    raise_error: bool = True,
    best_only: bool = False,
    feature_type: Optional[str] = None,
) -> List:
    """Map transcript coordinates to contig coordinates."""
    result = set()

    for feature_id in CM().transcript_ids(feature, feature_type, best_only):
        for pos in transcript_to_gene(
            feature_id,
            start,
            end,
            raise_error=False,
            best_only=best_only,
            feature_type=TRANSCRIPT,
        ):
            result |= _map(pos.gene_id, pos.start, pos.end, GENE, CONTIG, best_only=best_only)

    return _validate_result(result, "transcript", "contig", (feature, start, end), raise_error)


@lru_cache()
def transcript_to_exon(
    feature: str,
    start: int,
    end: Optional[int] = None,
    raise_error: bool = True,
    best_only: bool = False,
    feature_type: Optional[str] = None,
) -> List:
    """Map transcript coordinates to exon coordinates."""
    result = set()

    for feature_id in CM().transcript_ids(feature, feature_type, best_only):
        for pos in transcript_to_gene(
            feature_id,
            start,
            end,
            raise_error=False,
            best_only=best_only,
            feature_type=TRANSCRIPT,
        ):
            for pos2 in gene_to_exon(
                pos.gene_id,
                pos.start,
                pos.end,
                raise_error=False,
                best_only=best_only,
                feature_type=GENE,
            ):
                if pos2.transcript_id == feature_id:
                    result.add(pos2)

    return _validate_result(result, "transcript", "exon", (feature, start, end), raise_error)


@lru_cache()
def transcript_to_gene(
    feature: str,
    start: int,
    end: Optional[int] = None,
    raise_error: bool = True,
    best_only: bool = False,
    feature_type: Optional[str] = None,
) -> List:
    """Map transcript coordinates to gene coordinates."""
    result = set()

    for feature_id in CM().transcript_ids(feature, feature_type, best_only):
        result |= _map(feature_id, start, end, TRANSCRIPT, GENE, best_only=best_only)

    return _validate_result(result, "transcript", "gene", (feature, start, end), raise_error)


@lru_cache()
def transcript_to_protein(
    feature: str,
    start: int,
    end: Optional[int] = None,
    raise_error: bool = True,
    best_only: bool = False,
    feature_type: Optional[str] = None,
) -> List:
    """Map transcript coordinates to protein coordinates."""
    result = set()

    for feature_id in CM().transcript_ids(feature, feature_type, best_only):
        for pos in transcript_to_cds(
            feature_id,
            start,
            end,
            raise_error=False,
            best_only=best_only,
            feature_type=TRANSCRIPT,
        ):
            for pos2 in cds_to_protein(
                pos.transcript_id,
                pos.start,
                pos.end,
                raise_error=False,
                best_only=best_only,
                feature_type=TRANSCRIPT,
            ):
                result.add(pos2)

    return _validate_result(result, "transcript", "protein", (feature, start, end), raise_error)


@lru_cache()
def transcript_to_transcript(
    feature: str,
    start: int,
    end: Optional[int] = None,
    raise_error: bool = True,
    best_only: bool = False,
    feature_type: Optional[str] = None,
) -> List:
    """Map transcript coordinates."""
    result = set()

    for feature_id in CM().transcript_ids(feature, feature_type, best_only):
        result |= _map(feature_id, start, end, TRANSCRIPT, TRANSCRIPT, best_only=best_only)

    return _validate_result(result, "transcript", "transcript", (feature, start, end), raise_error)


def get_map_func(from_type: str, to_type: str) -> Callable:
    """Return a function for mapping between the specified feature types."""
    if from_type == CDS and to_type == CDS:
        return cds_to_cds
    elif from_type == CDS and to_type == CONTIG:
        return cds_to_contig
    elif from_type == CDS and to_type == EXON:
        return cds_to_exon
    elif from_type == CDS and to_type == GENE:
        return cds_to_gene
    elif from_type == CDS and to_type == PROTEIN:
        return cds_to_protein
    elif from_type == CDS and to_type == TRANSCRIPT:
        return cds_to_transcript
    elif from_type == CONTIG and to_type == CDS:
        return contig_to_cds
    elif from_type == CONTIG and to_type == CONTIG:
        return contig_to_contig
    elif from_type == CONTIG and to_type == EXON:
        return contig_to_exon
    elif from_type == CONTIG and to_type == GENE:
        return contig_to_gene
    elif from_type == CONTIG and to_type == PROTEIN:
        return contig_to_protein
    elif from_type == CONTIG and to_type == TRANSCRIPT:
        return contig_to_transcript
    elif from_type == EXON and to_type == CDS:
        return exon_to_cds
    elif from_type == EXON and to_type == CONTIG:
        return exon_to_contig
    elif from_type == EXON and to_type == EXON:
        return exon_to_exon
    elif from_type == EXON and to_type == GENE:
        return exon_to_gene
    elif from_type == EXON and to_type == PROTEIN:
        return exon_to_protein
    elif from_type == EXON and to_type == TRANSCRIPT:
        return exon_to_transcript
    elif from_type == GENE and to_type == CDS:
        return gene_to_cds
    elif from_type == GENE and to_type == CONTIG:
        return gene_to_contig
    elif from_type == GENE and to_type == EXON:
        return gene_to_exon
    elif from_type == GENE and to_type == GENE:
        return gene_to_gene
    elif from_type == GENE and to_type == PROTEIN:
        return gene_to_protein
    elif from_type == GENE and to_type == TRANSCRIPT:
        return gene_to_transcript
    elif from_type == PROTEIN and to_type == CDS:
        return protein_to_cds
    elif from_type == PROTEIN and to_type == CONTIG:
        return protein_to_contig
    elif from_type == PROTEIN and to_type == EXON:
        return protein_to_exon
    elif from_type == PROTEIN and to_type == GENE:
        return protein_to_gene
    elif from_type == PROTEIN and to_type == PROTEIN:
        return protein_to_protein
    elif from_type == PROTEIN and to_type == TRANSCRIPT:
        return protein_to_transcript
    elif from_type == TRANSCRIPT and to_type == CDS:
        return transcript_to_cds
    elif from_type == TRANSCRIPT and to_type == CONTIG:
        return transcript_to_contig
    elif from_type == TRANSCRIPT and to_type == EXON:
        return transcript_to_exon
    elif from_type == TRANSCRIPT and to_type == GENE:
        return transcript_to_gene
    elif from_type == TRANSCRIPT and to_type == PROTEIN:
        return transcript_to_protein
    elif from_type == TRANSCRIPT and to_type == TRANSCRIPT:
        return transcript_to_transcript
    elif from_type == to_type:
        return _no_convert
    else:
        raise ValueError(f"Cannot map '{from_type}' to '{to_type}'")


def _map(
    feature: str,
    start: int,
    end: Optional[int],
    from_type: str,
    to_type: str,
    best_only: bool = False,
) -> Set:
    """
    Template function for mapping a feature to the associated transcript(s), then
    converting the given coordinates.

    Args:
        feature (str): feature name or Ensembl ID
        start (int): first position relative to `feature`
        end (int or None): second position relative to `feature`
        from_type (str): coordinates relative to this type of feature (e.g. 'gene')
        to_type (str): map coordinates to this type of feature (e.g 'transcript')
        best_only (bool): only return results that map to a best transcript

    Returns:
        list: of converted coordinates, mapped to a `features` instance
    """
    result = set()
    skipped = set()

    logger.debug(
        f"Mapping {from_type} ({feature}, {start}, {end}) to {to_type}{' (best_only)' if best_only else ''}"
    )
    assert_valid_position(start, end)
    map_func = _get_pos_func(from_type, to_type)
    load_func = get_load_function(to_type)

    for transcript in CM().get_transcripts(feature, from_type):
        if best_only and not CM().is_canonical_transcript(transcript.transcript_id):
            skipped.add(transcript.transcript_id)
            continue

        try:
            position = map_func(transcript, start, end)
            feature = load_func(transcript, *position)
            result.add(feature)
        except ValueError as exc:
            logger.debug(exc)
            continue

    if skipped:
        logger.debug(f"Ignoring non-best transcripts: {', '.join(skipped)}")

    return result


# --------------------------------------
#  Position-mapping functions
# --------------------------------------
def _cpos_to_cpos(
    transcript: pyensembl.Transcript, start: int, end: Optional[int] = None
) -> Tuple[int, int]:
    """Assert the given position falls actually falls on the given transcript.

    Args:
        start (int): position relative to the CDS
        end (int): optional, second position relative to the CDS

    Returns:
        tuple of int: position relative to the CDS
    """

    def check(y):
        if not (1 <= y <= len(transcript.coding_sequence)):
            raise CdsOutOfRange(transcript, y)

    check(start)
    if end:
        check(end)
        return start, end
    else:
        return start, start


def _cpos_to_ppos(_: Any, start: int, end: Optional[int] = None) -> Tuple[int, int]:
    """Compute the equivalent protein position for a CDS position.

    Args:
        start (int): position relative to the CDS
        end (int): optional, second position relative to the CDS

    Returns:
        tuple of int: amino acid position
    """

    def convert(x):
        return floor((x - 1) / 3 + 1)

    pstart = convert(start)
    pend = convert(end) if end else pstart

    return pstart, pend


def _cpos_to_tpos(
    transcript: pyensembl.Transcript, start: int, end: Optional[int] = None
) -> Tuple[int, int]:
    """Compute the equivalent CDS position for a transcript position.

    Args:
        transcript: `pyensembl.Transcript` instance
        start (int): position relative to the CDS
        end (int): optional, second position relative to the CDS

    Returns:
        tuple of int: position relative to the transcript
    """

    def convert(x):
        y = transcript.first_start_codon_spliced_offset + x
        if not (1 <= y <= len(transcript.sequence)):
            raise TranscriptOutOfRange(transcript, x)
        return y

    tstart = convert(start)
    tend = convert(end) if end else tstart

    return tstart, tend


def _gpos_to_cpos(transcript: pyensembl.Transcript, start: int, end: int) -> Tuple[int, int]:
    """Compute the equivalent CDS position for an exon.

    Args:
        transcript: `pyensembl.Transcript` instance
        start (int): genomic coordinate
        end (int): optional, second genomic coordinate

    Returns:
        tuple of int: position relative to the CDS
    """
    for cds in sorted(transcript.coding_sequence_position_ranges):
        if end:
            if (cds[0] <= start <= cds[1]) or (cds[0] <= end <= cds[1]):
                gstart = max(start, cds[0])
                gend = min(end, cds[1])
                if start != gstart:
                    logger.debug(f"CDS start adjusted from {start} to {gstart}")
                if end != gend:
                    logger.debug(f"CDS end adjusted from {end} to {gend}")
                break
        else:
            if cds[0] <= start <= cds[1]:
                gstart = start
                gend = start
                break
    else:
        raise CdsOutOfRange(transcript, start)

    tstart, tend = _gpos_to_tpos(transcript, gstart, gend)
    cstart, cend = _tpos_to_cpos(transcript, tstart, tend)

    return cstart, cend


def _gpos_to_gpos(
    transcript: pyensembl.Transcript, start: int, end: Optional[int] = None
) -> Tuple[int, int]:
    """Assert the given position falls actually falls on the given gene.

    Args:
        transcript: `pyensembl.Transcript` instance
        start (int): genomic coordinate
        end (int): optional, second genomic coordinate

    Returns:
        tuple of int: genomic coordinate
    """
    gene = CM().get_genes(transcript.gene_id, GENE)[0]

    def check(y):
        if not (gene.start <= y <= gene.end):
            raise GeneOutOfRange(gene, y)

    check(start)
    if end:
        check(end)
        return start, end
    else:
        return start, start


def _gpos_to_tpos(
    transcript: pyensembl.Transcript, start: int, end: Optional[int] = None
) -> Tuple[int, int]:
    """Compute the equivalent transcript position for a gene position.

    Args:
        transcript: `pyensembl.Transcript` instance
        start (int): genomic coordinate
        end (int): optional, second genomic coordinate

    Returns:
        tuple of int: position relative to the transcript
    """

    def convert(x):
        return transcript.spliced_offset(x) + 1

    tstart = convert(start)
    tend = convert(end) if end else tstart

    return tstart, tend


def _ppos_to_cpos(
    transcript: pyensembl.Transcript, start: int, end: Optional[int] = None
) -> Tuple[int, int]:
    """Compute the equivalent protein position for a CDS position.

    Args:
        transcript: `pyensembl.Transcript` instance
        start (int): amino acid position
        end (int): optional, second amino acid position

    Returns:
        tuple of int: CDS position of the first base of the codon
    """

    def convert(x):
        y = (x - 1) * 3 + 1
        if not (1 <= y <= len(transcript.coding_sequence)):
            raise CdsOutOfRange(transcript, y)
        else:
            return y

    cstart = convert(start)
    # add 2 to the end position to get the end of the codon
    cend = (convert(end) if end else cstart) + 2

    return cstart, cend


def _ppos_to_ppos(
    transcript: pyensembl.Transcript, start: int, end: Optional[int] = None
) -> Tuple[int, int]:
    """Assert the given position falls actually falls on the given protein.

    Args:
        transcript: `pyensembl.Transcript` instance
        start (int): amino acid position
        end (int): optional, second amino acid position

    Returns:
        tuple of int: amino acid position
    """

    def check(y):
        if not (1 <= y <= len(transcript.coding_sequence) // 3):
            raise TranscriptOutOfRange(transcript, y)

    check(start)
    if end:
        check(end)
        return start, end
    else:
        return start, start


def _gpos_to_epos(
    transcript: pyensembl.Transcript, start: int, end: Optional[int] = None
) -> Tuple[int, int]:
    """Return the genomic coordinates of the exon that contains a genomic coordinate.

    Args:
        transcript: `pyensembl.Transcript` instance
        start (int): genomic coordinate
        end (int): optional, second genomic coordinate

    Returns:
        tuple of int: genomic coordinates of the exon
    """

    def convert(x):
        for exon in sorted(transcript.exons, key=lambda x: x.start):
            if exon.start <= x <= exon.end:
                return exon.start, exon.end
        else:
            raise ExonOutOfRange(transcript, x)

    exon1 = convert(start)
    exon2 = convert(end) if end else exon1

    if exon1 != exon2:
        raise ValueError(
            f"{start} and {end} map to different exons of {transcript.transcript_id} ({exon1}, {exon2})"
        )
    else:
        return exon1


def _tpos_to_cpos(
    transcript: pyensembl.Transcript, start: int, end: Optional[int] = None
) -> Tuple[int, int]:
    """Compute the equivalent transcript position for a CDS position.

    Args:
        transcript: `pyensembl.Transcript` instance
        start (int): position relative to the transcript
        end (int): optional, second position relative to the transcript

    Returns:
        tuple of int: position relative to the CDS
    """

    def convert(x):
        y = x - transcript.first_start_codon_spliced_offset
        if not (1 <= y <= len(transcript.coding_sequence)):
            raise CdsOutOfRange(transcript, x)
        return y

    cstart = convert(start)
    cend = convert(end) if end else cstart

    return cstart, cend


def _tpos_to_gpos(
    transcript: pyensembl.Transcript, start: int, end: Optional[int] = None
) -> Tuple[int, int]:
    """Compute the equivalent gene position for a transcript position.

    Args:
        transcript: `pyensembl.Transcript` instance
        start (int): position relative to the transcript
        end (int): optional, second position relative to the transcript

    Returns:
        tuple of int: genomic coordinate
    """
    # make sure all the ranges are sorted from smallest to biggest
    ranges = sorted([sorted(i) for i in transcript.exon_intervals])

    # for transcripts on the negative strand, start counting from the last position
    if transcript.on_negative_strand:
        ranges = ranges[::-1]

    def convert(x):
        remain = x - 1
        for i in ranges:
            length = i[1] - i[0] + 1
            if remain >= length:
                remain -= length
            else:
                if transcript.on_positive_strand:
                    return i[0] + remain
                elif transcript.on_negative_strand:
                    return i[1] - remain
        else:
            raise TranscriptOutOfRange(transcript, x)

    gstart = convert(start)
    gend = convert(end) if end else gstart

    return gstart, gend


def _tpos_to_tpos(
    transcript: pyensembl.Transcript, start: int, end: Optional[int] = None
) -> Tuple[int, int]:
    """Assert the given position falls actually falls on the given transcript.

    Args:
        transcript: `pyensembl.Transcript` instance
        start (int): position relative to the transcript
        end (int): optional, second position relative to the transcript

    Returns:
        tuple of int: position relative to the transcript
    """

    def check(y):
        if not (1 <= y <= len(transcript.sequence)):
            raise TranscriptOutOfRange(transcript, y)

    check(start)
    if end:
        check(end)
        return start, end
    else:
        return start, start


def _no_convert(_: Any, start: int, end: Optional[int] = None) -> Tuple[int, int]:
    """Dummy function for when no conversion is needed."""
    if end:
        return start, end
    else:
        return start, start


def _get_pos_func(from_type: str, to_type: str) -> Callable:
    if from_type == CDS and to_type == CDS:
        return _cpos_to_cpos
    elif from_type == CDS and to_type == PROTEIN:
        return _cpos_to_ppos
    elif from_type == CDS and to_type == TRANSCRIPT:
        return _cpos_to_tpos
    elif from_type == EXON and to_type == TRANSCRIPT:
        return _gpos_to_tpos
    elif from_type == GENE and to_type == CDS:
        return _gpos_to_cpos
    elif from_type == GENE and to_type == CONTIG:
        return _no_convert
    elif from_type == GENE and to_type == EXON:
        return _gpos_to_epos
    elif from_type == GENE and to_type == GENE:
        return _gpos_to_gpos
    elif from_type == GENE and to_type == TRANSCRIPT:
        return _gpos_to_tpos
    elif from_type == PROTEIN and to_type == CDS:
        return _ppos_to_cpos
    elif from_type == PROTEIN and to_type == PROTEIN:
        return _ppos_to_ppos
    elif from_type == TRANSCRIPT and to_type == CDS:
        return _tpos_to_cpos
    elif from_type == TRANSCRIPT and to_type == GENE:
        return _tpos_to_gpos
    elif from_type == TRANSCRIPT and to_type == TRANSCRIPT:
        return _tpos_to_tpos
    elif from_type == to_type:
        return _no_convert
    else:
        raise ValueError(f"Cannot map '{from_type}' position to '{to_type}'")


# --------------------------------------
#  Extra functions
# --------------------------------------
def _transcript_ids_with_exon(exon_id: str, position: Optional[int] = None) -> List[str]:
    all_transcript_ids = []
    transcripts = CM().get_transcripts(exon_id, EXON)

    for tobj in transcripts:
        # assert the exon is the nth exon of the transcript
        if position:
            try:
                if exon_id == tobj.exons[position - 1].exon_id:
                    all_transcript_ids.append(tobj.transcript_id)
            except IndexError:
                continue
        # otherwise, add all the transcripts
        else:
            all_transcript_ids.append(tobj.transcript_id)

    return all_transcript_ids


def _validate_result(
    result: Set,
    from_key: str,
    to_key: str,
    feature: Tuple,
    raise_error: bool = True,
) -> List:
    if not result and raise_error:
        raise ValueError(f"Could not map {from_key} {feature} to {to_key}")

    return sorted(result)
