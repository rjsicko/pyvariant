from math import floor
from typing import Any, Callable, Optional, Tuple

import pyensembl
from logzero import logger

from .exceptions import CdsOutOfRange, ExonOutOfRange, GeneOutOfRange, TranscriptOutOfRange
from .ids import get_genes


def get_convert_func(from_type: str, to_type: str) -> Callable:
    if from_type == "cds" and to_type == "cds":
        return _cpos_to_cpos
    elif from_type == "cds" and to_type == "protein":
        return _cpos_to_ppos
    elif from_type == "cds" and to_type == "transcript":
        return _cpos_to_tpos
    elif from_type == "exon" and to_type == "transcript":
        return _gpos_to_tpos
    elif from_type == "gene" and to_type == "cds":
        return _gpos_to_cpos
    elif from_type == "gene" and to_type == "contig":
        return _no_convert
    elif from_type == "gene" and to_type == "exon":
        return _gpos_to_epos
    elif from_type == "gene" and to_type == "gene":
        return _gpos_to_gpos
    elif from_type == "gene" and to_type == "transcript":
        return _gpos_to_tpos
    elif from_type == "protein" and to_type == "cds":
        return _ppos_to_cpos
    elif from_type == "protein" and to_type == "protein":
        return _ppos_to_ppos
    elif from_type == "transcript" and to_type == "cds":
        return _tpos_to_cpos
    elif from_type == "transcript" and to_type == "gene":
        return _tpos_to_gpos
    elif from_type == "transcript" and to_type == "transcript":
        return _tpos_to_tpos
    elif from_type == to_type:
        return _no_convert
    else:
        raise ValueError(f"Cannot convert '{from_type}' to '{to_type}'")


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
    if end:
        pend = convert(end)
    else:
        pend = pstart

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
    if end:
        tend = convert(end)
    else:
        tend = tstart

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
    gene = get_genes(transcript.gene_id, "gene")[0]

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
    if end:
        tend = convert(end)
    else:
        tend = tstart

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
    if end:
        cend = convert(end) + 2
    else:
        cend = cstart + 2

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
    if end:
        exon2 = convert(end)
        if exon1 != exon2:
            raise ValueError(
                f"{start} and {end} map to different exons of {transcript.transcript_id} ({exon1}, {exon2})"
            )

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
    if end:
        cend = convert(end)
    else:
        cend = cstart

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
    if end:
        gend = convert(end)
    else:
        gend = gstart

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
