from functools import lru_cache
from math import floor
from typing import Any, Callable, List, Optional, Set, Tuple

from logzero import logger

from .constants import CDNA, CONTIG, EXON, GENE, PROTEIN, STRAND_SYMBOLS, TRANSCRIPT
from .core import (
    CdnaRecord,
    ContigRecord,
    ExonRecord,
    GeneRecord,
    ProteinRecord,
    Record,
    RecordType,
    TranscriptRecord,
)
from .core import instance as CM
from .exceptions import CdnaOutOfRange, ExonOutOfRange, GeneOutOfRange, TranscriptOutOfRange
from .utils import assert_valid_position


# --------------------------------------
#  General mapping functions
# --------------------------------------
@lru_cache()
def cdna_to_cdna(
    feature: str,
    start: int,
    end: Optional[int] = None,
    raise_error: bool = True,
    best_only: bool = False,
    feature_type: Optional[str] = None,
) -> List:
    """Map cDNA coordinates."""
    result = set()
    end = end if end is not None else start

    for cdna in CM().cdna_info(feature, feature_type):
        if 1 <= start <= cdna.length and 1 <= end <= cdna.length:
            new_cdna = _convert_record(CdnaRecord, cdna, start, end)
            result.add(new_cdna)

    return _validate_result(result, "cDNA", "cDNA", (feature, start, end), raise_error)


@lru_cache()
def cdna_to_contig(
    feature: str,
    start: int,
    end: Optional[int] = None,
    raise_error: bool = True,
    best_only: bool = False,
    feature_type: Optional[str] = None,
) -> List:
    """Map cDNA coordinates to contig coordinates."""
    result = set()
    end = end if end is not None else start

    for cdna in cdna_to_cdna(
        feature,
        start,
        end=end,
        raise_error=False,
        best_only=best_only,
        feature_type=feature_type,
    ):
        contig_start = cdna.unspliced_position(start)
        contig_end = cdna.unspliced_position(end)
        contig_start, contig_end = sorted((contig_start, contig_end))
        for contig in CM().contig_info(cdna.transcript_id, TRANSCRIPT):
            new_contig = _convert_record(ContigRecord, contig, contig_start, contig_end)
            result.add(new_contig)

    return _validate_result(result, "cDNA", "contig", (feature, start, end), raise_error)


@lru_cache()
def cdna_to_exon(
    feature: str,
    start: int,
    end: Optional[int] = None,
    raise_error: bool = True,
    best_only: bool = False,
    feature_type: Optional[str] = None,
) -> List:
    """Map cDNA coordinates to exon coordinates."""
    result = set()
    end = end if end is not None else start

    for transcript in cdna_to_transcript(
        feature,
        start,
        end=end,
        raise_error=False,
        best_only=best_only,
        feature_type=feature_type,
    ):
        for exon in transcript_to_exon(
            transcript.transcript_id,
            transcript.start,
            end=transcript.end,
            raise_error=False,
            best_only=best_only,
            feature_type=feature_type,
        ):
            result.add(exon)

    return _validate_result(result, "cDNA", "exon", (feature, start, end), raise_error)


@lru_cache()
def cdna_to_gene(
    feature: str,
    start: int,
    end: Optional[int] = None,
    raise_error: bool = True,
    best_only: bool = False,
    feature_type: Optional[str] = None,
) -> List:
    """Map cDNA coordinates to gene coordinates."""
    result = set()
    end = end if end is not None else start

    for transcript in cdna_to_transcript(
        feature,
        start,
        end=end,
        raise_error=False,
        best_only=best_only,
        feature_type=feature_type,
    ):
        for gene in transcript_to_gene(
            transcript.transcript_id,
            transcript.start,
            end=transcript.end,
            raise_error=False,
            best_only=best_only,
            feature_type=feature_type,
        ):
            result.add(gene)

    return _validate_result(result, "cDNA", "gene", (feature, start, end), raise_error)


@lru_cache()
def cdna_to_protein(
    feature: str,
    start: int,
    end: Optional[int] = None,
    raise_error: bool = True,
    best_only: bool = False,
    feature_type: Optional[str] = None,
) -> List:
    """Map cDNA coordinates to protein coordinates."""
    result = set()
    end = end if end is not None else start

    for cdna in cdna_to_cdna(
        feature,
        start,
        end=end,
        raise_error=False,
        best_only=best_only,
        feature_type=feature_type,
    ):
        for protein in CM().protein_info(cdna.transcript_id, CDNA):
            protein_start = floor((start - 1) / 3 + 1)
            protein_end = floor((end - 1) / 3 + 1)
            new_protein = _convert_record(ProteinRecord, protein, protein_start, protein_end)
            result.add(new_protein)

    return _validate_result(result, "cDNA", "protein", (feature, start, end), raise_error)


@lru_cache()
def cdna_to_transcript(
    feature: str,
    start: int,
    end: Optional[int] = None,
    raise_error: bool = True,
    best_only: bool = False,
    feature_type: Optional[str] = None,
) -> List:
    """Map cDNA coordinates to transcript coordinates."""
    result = set()
    end = end if end is not None else start

    for cdna in cdna_to_cdna(
        feature,
        start,
        end=end,
        raise_error=False,
        best_only=best_only,
        feature_type=feature_type,
    ):
        for transcript in CM().transcript_info(cdna.transcript_id, TRANSCRIPT):
            transcript_start = start + transcript.first_start_codon_spliced_offset()
            transcript_end = end + transcript.first_start_codon_spliced_offset()
            new_transcript = _convert_record(
                TranscriptRecord, transcript, transcript_start, transcript_end
            )
            result.add(new_transcript)

    return _validate_result(result, "cDNA", "transcript", (feature, start, end), raise_error)


@lru_cache()
def contig_to_cdna(
    feature: str,
    start: int,
    end: Optional[int] = None,
    raise_error: bool = True,
    best_only: bool = False,
    feature_type: Optional[str] = None,
    strand: Optional[str] = None,
) -> List:
    """Map contig coordinates to cDNA coordinates."""
    result = set()

    if strand not in STRAND_SYMBOLS:
        raise ValueError(f"strand must be one of {STRAND_SYMBOLS}")

    for feature_id in CM().contig_ids(feature, feature_type):
        for pos in contig_to_gene(
            feature_id,
            start,
            end,
            raise_error=False,
            best_only=best_only,
            feature_type=CONTIG,
            strand=strand,
        ):
            for pos2 in gene_to_cdna(
                pos.gene_id,
                pos.start,
                pos.end,
                raise_error=False,
                best_only=best_only,
                feature_type=GENE,
            ):
                result.add(pos2)

    return _validate_result(result, "contig", "cDNA", (feature, start, end), raise_error)


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
    end = end if end is not None else start

    for contig in CM().contig_info(feature, feature_type):
        if 1 <= start <= contig.length and 1 <= end <= contig.length:
            new_contig = _convert_record(ContigRecord, contig, start, end)
            result.add(new_contig)

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

    for feature_id in CM().contig_ids(feature, feature_type):
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

    for gene in CM().gene_info(feature, feature_type):
        if gene.start <= start <= gene.end and gene.start <= end <= gene.end:
            new_gene = _convert_record(GeneRecord, gene, start, end)
            result.add(new_gene)

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

    for feature_id in CM().contig_ids(feature, feature_type):
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

    for feature_id in CM().contig_ids(feature, feature_type):
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
def exon_to_cdna(
    feature: str,
    position: Optional[int] = None,
    _: Optional[int] = None,
    raise_error: bool = True,
    best_only: bool = False,
    feature_type: Optional[str] = None,
) -> List:
    """Map an exon to cDNA coordinates."""
    result = set()

    for feature_id in CM().exon_ids(feature, feature_type):
        for pos in exon_to_gene(
            feature_id, position, raise_error=False, best_only=best_only, feature_type=EXON
        ):
            result |= _map(pos.gene_id, pos.start, pos.end, GENE, CDNA, best_only=best_only)

    return _validate_result(result, "exon", "cDNA", (feature, position), raise_error)


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

    for feature_id in CM().exon_ids(feature, feature_type):
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

    for exon in CM().exon_info(feature, feature_type):
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

    for feature_id in CM().exon_ids(feature, feature_type):
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

    for feature_id in CM().exon_ids(feature, feature_type):
        for pos in exon_to_cdna(
            feature_id, position, raise_error=False, best_only=best_only, feature_type=EXON
        ):
            for pos2 in cdna_to_protein(
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

    for exon in CM().exon_info(feature, feature_type):
        # Filter out transcripts that don't contain the query exon.
        valid_transcripts = _transcript_ids_with_exon(exon.exon_id, position)
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
def gene_to_cdna(
    feature: str,
    start: int,
    end: Optional[int] = None,
    raise_error: bool = True,
    best_only: bool = False,
    feature_type: Optional[str] = None,
) -> List:
    """Map gene coordinates to cDNA coordinates."""
    result = set()

    if feature_type == CONTIG or CM().is_contig(feature):
        return contig_to_cdna((feature, start, end), raise_error, best_only, feature_type)

    for feature_id in CM().gene_ids(feature, feature_type):
        for pos in gene_to_transcript(
            feature_id, start, end, raise_error=False, best_only=best_only, feature_type=GENE
        ):
            for pos2 in transcript_to_cdna(
                pos.transcript_id,
                pos.start,
                pos.end,
                raise_error=False,
                best_only=best_only,
                feature_type=TRANSCRIPT,
            ):
                result.add(pos2)

    return _validate_result(result, "gene", "cDNA", (feature, start, end), raise_error)


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
        for pos in gene_to_cdna(
            feature_id, start, end, raise_error=False, best_only=best_only, feature_type=GENE
        ):
            for pos2 in cdna_to_protein(
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
def protein_to_cdna(
    feature: str,
    start: int,
    end: Optional[int] = None,
    raise_error: bool = True,
    best_only: bool = False,
    feature_type: Optional[str] = None,
) -> List:
    """Map protein coordinates to cDNA coordinates."""
    result = set()

    for feature_id in CM().protein_ids(feature, feature_type):
        result |= _map(feature_id, start, end, PROTEIN, CDNA, best_only=best_only)

    return _validate_result(result, "protein", "cDNA", (feature, start, end), raise_error)


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

    for feature_id in CM().protein_ids(feature, feature_type):
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

    for feature_id in CM().protein_ids(feature, feature_type):
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

    for feature_id in CM().protein_ids(feature, feature_type):
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

    for feature_id in CM().protein_ids(feature, feature_type):
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

    for feature_id in CM().protein_ids(feature, feature_type):
        for pos in protein_to_cdna(
            feature_id, start, end, raise_error=False, best_only=best_only, feature_type=PROTEIN
        ):
            for pos2 in cdna_to_transcript(
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
def transcript_to_cdna(
    feature: str,
    start: int,
    end: Optional[int] = None,
    raise_error: bool = True,
    best_only: bool = False,
    feature_type: Optional[str] = None,
) -> List:
    """Map transcript coordinates to cDNA coordinates."""
    result = set()
    end = end if end is not None else start

    for transcript in transcript_to_transcript(
        feature,
        start,
        end=end,
        raise_error=False,
        best_only=best_only,
        feature_type=feature_type,
    ):
        for cdna in CM().cdna_info(transcript.transcript_id, TRANSCRIPT):
            offset = cdna.first_start_codon_spliced_offset()
            cdna_start = start - offset
            cdna_end = end - offset
            new_cdna = _convert_record(CdnaRecord, cdna, cdna_start, cdna_end)
            result.add(new_cdna)

    return _validate_result(result, "transcript", "cDNA", (feature, start, end), raise_error)


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
    end = end if end is not None else start

    for transcript in transcript_to_transcript(
        feature,
        start,
        end=end,
        raise_error=False,
        best_only=best_only,
        feature_type=feature_type,
    ):
        contig_start = transcript.unspliced_position(start)
        contig_end = transcript.unspliced_position(end)
        contig_start, contig_end = sorted((contig_start, contig_end))
        for contig in CM().contig_info(transcript.transcript_id, TRANSCRIPT):
            new_contig = _convert_record(ContigRecord, contig, contig_start, contig_end)
            result.add(new_contig)

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
    end = end if end is not None else start

    for transcript in transcript_to_transcript(
        feature,
        start,
        end=end,
        raise_error=False,
        best_only=best_only,
        feature_type=feature_type,
    ):
        remain_start = start - 1
        remain_end = end - 1
        exons = CM().exon_info(transcript.transcript_id, TRANSCRIPT)
        exons = sorted(exons, key=lambda x: (x.start, x.end))
        exons = exons[::-1] if transcript.on_negative_strand() else exons
        for i in exons:
            # if both start and end fall inside the same exons, add the exon
            if remain_start - i.length <= 0 and remain_end - i.length <= 0:
                result.add(i)
                break
            # if the start and end fall inside different exons, raise a warning
            elif remain_start - i.length <= 0 or remain_end - i.length <= 0:
                logger.debug(f"{transcript.transcript_id} ({start}, {end}) are on different exons")
                break
            # otherwise, keep checking
            else:
                remain_start -= i.length
                remain_end -= i.length

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
    end = end if end is not None else start

    for transcript in transcript_to_transcript(
        feature,
        start,
        end=end,
        raise_error=False,
        best_only=best_only,
        feature_type=feature_type,
    ):
        gene_start = transcript.unspliced_position(start)
        gene_end = transcript.unspliced_position(end)
        gene_start, gene_end = sorted((gene_start, gene_end))
        for gene in CM().gene_info(transcript.transcript_id, TRANSCRIPT):
            new_gene = _convert_record(GeneRecord, gene, gene_start, gene_end)
            result.add(new_gene)

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
    end = end if end is not None else start

    for cdna in transcript_to_cdna(
        feature,
        start,
        end=end,
        raise_error=False,
        best_only=best_only,
        feature_type=feature_type,
    ):
        for protein in CM().protein_info(cdna.transcript_id, CDNA):
            protein_start = floor((cdna.start - 1) / 3 + 1)
            protein_end = floor((cdna.end - 1) / 3 + 1)
            protein = _convert_record(ProteinRecord, protein, protein_start, protein_end)
            result.add(protein)

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
    end = end if end is not None else start

    for transcript in CM().transcript_info(feature, feature_type):
        if 1 <= start <= transcript.length and 1 <= end <= transcript.length:
            new_transcript = _convert_record(TranscriptRecord, transcript, start, end)
            result.add(new_transcript)
        else:
            logger.debug(f"({start}, {end}) is outside transcript (1, {transcript.length})")

    return _validate_result(result, "transcript", "transcript", (feature, start, end), raise_error)


def get_map_func(from_type: str, to_type: str) -> Callable:
    """Return a function for mapping between the specified feature types."""
    if from_type == CDNA and to_type == CDNA:
        return cdna_to_cdna
    elif from_type == CDNA and to_type == CONTIG:
        return cdna_to_contig
    elif from_type == CDNA and to_type == EXON:
        return cdna_to_exon
    elif from_type == CDNA and to_type == GENE:
        return cdna_to_gene
    elif from_type == CDNA and to_type == PROTEIN:
        return cdna_to_protein
    elif from_type == CDNA and to_type == TRANSCRIPT:
        return cdna_to_transcript
    elif from_type == CONTIG and to_type == CDNA:
        return contig_to_cdna
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
    elif from_type == EXON and to_type == CDNA:
        return exon_to_cdna
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
    elif from_type == GENE and to_type == CDNA:
        return gene_to_cdna
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
    elif from_type == PROTEIN and to_type == CDNA:
        return protein_to_cdna
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
    elif from_type == TRANSCRIPT and to_type == CDNA:
        return transcript_to_cdna
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
    info_func, load_class = _get_load_func(to_type)

    for transcript in CM().transcript_info(feature, from_type):
        if best_only and not CM().is_canonical_transcript(transcript.transcript_id):
            skipped.add(transcript.transcript_id)
            continue

        try:
            start_new, end_new = map_func(transcript, start, end)
            for record in info_func(transcript.transcript_id, TRANSCRIPT):
                new_record = _convert_record(load_class, record, start_new, end_new)
                if isinstance(new_record, ExonRecord) and new_record.exon_id == feature:
                    result.add(new_record)
        except ValueError as exc:
            logger.debug(exc)
            continue

    if skipped:
        logger.debug(f"Ignoring non-best transcripts: {', '.join(skipped)}")

    return result


# --------------------------------------
#  Loader functions
# --------------------------------------
def _get_load_func(to_type: str) -> Tuple[Callable, RecordType]:
    info_func: Callable
    load_class: RecordType

    if to_type == CDNA:
        info_func = CM().cdna_info
        load_class = CdnaRecord
    elif to_type == CONTIG:
        info_func = CM().contig_info
        load_class = ContigRecord
    elif to_type == EXON:
        info_func = CM().exon_info
        load_class = ExonRecord
    elif to_type == GENE:
        info_func = CM().gene_info
        load_class = GeneRecord
    elif to_type == PROTEIN:
        info_func = CM().gene_info
        load_class = ProteinRecord
    elif to_type == TRANSCRIPT:
        info_func = CM().transcript_info
        load_class = TranscriptRecord
    else:
        raise ValueError(f"Cannot get info for '{to_type}'")

    return info_func, load_class


def _convert_record(load_class: RecordType, record: Record, start: int, end: int) -> Record:
    values = record.__dict__.copy()
    _ = values.pop("start")
    _ = values.pop("end")
    new_record = load_class(start=start, end=end, **values)

    return new_record


# --------------------------------------
#  Position-mapping functions
# --------------------------------------
def _cpos_to_cpos(
    transcript: TranscriptRecord, start: int, end: Optional[int] = None
) -> Tuple[int, int]:
    """Assert the given position falls actually falls on the given transcript.

    Args:
        start (int): position relative to the cDNA
        end (int): optional, second position relative to the cDNA

    Returns:
        tuple of int: position relative to the cDNA
    """

    def check(y):
        if not (1 <= y <= len(CM().cdna_sequence(transcript.transcript_id))):
            raise CdnaOutOfRange(transcript, y)

    check(start)
    if end:
        check(end)
        return start, end
    else:
        return start, start


def _cpos_to_ppos(_: TranscriptRecord, start: int, end: Optional[int] = None) -> Tuple[int, int]:
    """Compute the equivalent protein position for a cDNA position.

    Args:
        start (int): position relative to the cDNA
        end (int): optional, second position relative to the cDNA

    Returns:
        tuple of int: amino acid position
    """

    def convert(x):
        return floor((x - 1) / 3 + 1)

    pstart = convert(start)
    pend = convert(end) if end else pstart

    return pstart, pend


def _cpos_to_tpos(
    transcript: TranscriptRecord, start: int, end: Optional[int] = None
) -> Tuple[int, int]:
    """Compute the equivalent transcript position for a cDNA position.

    Args:
        transcript: `TranscriptRecord` instance
        start (int): position relative to the cDNA
        end (int): optional, second position relative to the cDNA

    Returns:
        tuple of int: position relative to the transcript
    """

    def convert(x):
        y = x + transcript.first_start_codon_spliced_offset()
        if not (1 <= y <= transcript.length):
            raise TranscriptOutOfRange(transcript, x)
        return y

    tstart = convert(start)
    tend = convert(end) if end else tstart

    return tstart, tend


def _gpos_to_cpos(transcript: TranscriptRecord, start: int, end: int) -> Tuple[int, int]:
    """Compute the equivalent cDNA position for an exon.

    Args:
        transcript: `TranscriptRecord` instance
        start (int): genomic coordinate
        end (int): optional, second genomic coordinate

    Returns:
        tuple of int: position relative to the cDNA
    """
    for cdna in transcript.cds_intervals():
        if end:
            if (cdna[0] <= start <= cdna[1]) or (cdna[0] <= end <= cdna[1]):
                gstart = max(start, cdna[0])
                gend = min(end, cdna[1])
                if start != gstart:
                    logger.debug(f"cDNA start adjusted from {start} to {gstart}")
                if end != gend:
                    logger.debug(f"cDNA end adjusted from {end} to {gend}")
                break
        else:
            if cdna[0] <= start <= cdna[1]:
                gstart = start
                gend = start
                break
    else:
        raise CdnaOutOfRange(transcript, start)

    tstart, tend = _gpos_to_tpos(transcript, gstart, gend)
    cstart, cend = _tpos_to_cpos(transcript, tstart, tend)

    return cstart, cend


def _gpos_to_gpos(
    transcript: TranscriptRecord, start: int, end: Optional[int] = None
) -> Tuple[int, int]:
    """Assert the given position falls actually falls on the given gene.

    Args:
        transcript: `TranscriptRecord` instance
        start (int): genomic coordinate
        end (int): optional, second genomic coordinate

    Returns:
        tuple of int: genomic coordinate
    """
    gene = CM().gene_info(transcript.gene_id, GENE)[0]

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
    transcript: TranscriptRecord, start: int, end: Optional[int] = None
) -> Tuple[int, int]:
    """Compute the equivalent transcript position for a gene position.

    Args:
        transcript: `TranscriptRecord` instance
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
    transcript: TranscriptRecord, start: int, end: Optional[int] = None
) -> Tuple[int, int]:
    """Compute the equivalent protein position for a cDNA position.

    Args:
        transcript: `TranscriptRecord` instance
        start (int): amino acid position
        end (int): optional, second amino acid position

    Returns:
        tuple of int: cDNA position of the first base of the codon
    """

    def convert(x):
        y = (x - 1) * 3 + 1
        if not (1 <= y <= len(CM().cdna_sequence(transcript.transcript_id))):
            raise CdnaOutOfRange(transcript, y)
        else:
            return y

    cstart = convert(start)
    # add 2 to the end position to get the end of the codon
    cend = (convert(end) if end else cstart) + 2

    return cstart, cend


def _ppos_to_ppos(
    transcript: TranscriptRecord, start: int, end: Optional[int] = None
) -> Tuple[int, int]:
    """Assert the given position falls actually falls on the given protein.

    Args:
        transcript: `TranscriptRecord` instance
        start (int): amino acid position
        end (int): optional, second amino acid position

    Returns:
        tuple of int: amino acid position
    """

    def check(y):
        if not (1 <= y <= len(CM().cdna_sequence(transcript.transcript_id)) // 3):
            raise TranscriptOutOfRange(transcript, y)

    check(start)
    if end:
        check(end)
        return start, end
    else:
        return start, start


def _gpos_to_epos(
    transcript: TranscriptRecord, start: int, end: Optional[int] = None
) -> Tuple[int, int]:
    """Return the genomic coordinates of the exon that contains a genomic coordinate.

    Args:
        transcript: `TranscriptRecord` instance
        start (int): genomic coordinate
        end (int): optional, second genomic coordinate

    Returns:
        tuple of int: genomic coordinates of the exon
    """

    def convert(x):
        for exon in transcript.exons():
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
    transcript: TranscriptRecord, start: int, end: Optional[int] = None
) -> Tuple[int, int]:
    """Compute the equivalent cDNA position for a transcript position.

    Args:
        transcript: `TranscriptRecord` instance
        start (int): position relative to the transcript
        end (int): optional, second position relative to the transcript

    Returns:
        tuple of int: position relative to the cDNA
    """

    def convert(x):
        y = x - transcript.first_start_codon_spliced_offset()
        if not (1 <= y <= len(CM().cdna_sequence(transcript.transcript_id))):
            raise CdnaOutOfRange(transcript, x)
        return y

    cstart = convert(start)
    cend = convert(end) if end else cstart

    return cstart, cend


def _tpos_to_gpos(
    transcript: TranscriptRecord, start: int, end: Optional[int] = None
) -> Tuple[int, int]:
    """Compute the equivalent gene position for a transcript position.

    Args:
        transcript: `TranscriptRecord` instance
        start (int): position relative to the transcript
        end (int): optional, second position relative to the transcript

    Returns:
        tuple of int: genomic coordinate
    """
    gstart = transcript.unspliced_position(start)
    gend = transcript.unspliced_position(end) if end else gstart
    if gend < gstart:
        gstart, gend = gend, gstart

    return gstart, gend


def _tpos_to_tpos(
    transcript: TranscriptRecord, start: int, end: Optional[int] = None
) -> Tuple[int, int]:
    """Assert the given position falls actually falls on the given transcript.

    Args:
        transcript: `TranscriptRecord` instance
        start (int): position relative to the transcript
        end (int): optional, second position relative to the transcript

    Returns:
        tuple of int: position relative to the transcript
    """

    def check(y):
        if not (1 <= y <= transcript.length):
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
    if from_type == CDNA and to_type == CDNA:
        return _cpos_to_cpos
    elif from_type == CDNA and to_type == PROTEIN:
        return _cpos_to_ppos
    elif from_type == CDNA and to_type == TRANSCRIPT:
        return _cpos_to_tpos
    elif from_type == EXON and to_type == TRANSCRIPT:
        return _gpos_to_tpos
    elif from_type == GENE and to_type == CDNA:
        return _gpos_to_cpos
    elif from_type == GENE and to_type == CONTIG:
        return _no_convert
    elif from_type == GENE and to_type == EXON:
        return _gpos_to_epos
    elif from_type == GENE and to_type == GENE:
        return _gpos_to_gpos
    elif from_type == GENE and to_type == TRANSCRIPT:
        return _gpos_to_tpos
    elif from_type == PROTEIN and to_type == CDNA:
        return _ppos_to_cpos
    elif from_type == PROTEIN and to_type == PROTEIN:
        return _ppos_to_ppos
    elif from_type == TRANSCRIPT and to_type == CDNA:
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
    transcripts = CM().transcript_info(exon_id, EXON)

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
