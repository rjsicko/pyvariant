import logging
from math import floor

from .exceptions import CdsOutOfRange, ExonOutOfRange, TranscriptOutOfRange


def get_map_function(from_type, to_type):
    if from_type == to_type:
        raise ValueError("`from_type` and `to_type` must be different!")
    elif from_type == "cds" and to_type == "protein":
        return _cpos_to_ppos
    elif from_type == "cds" and to_type == "transcript":
        return _cpos_to_tpos
    elif from_type == "exon" and to_type == "cds":
        return _epos_to_cpos
    elif from_type == "exon" and to_type == "transcript":
        return _gpos_to_tpos
    elif from_type == "gene" and to_type == "exon":
        return _gpos_to_epos
    elif from_type == "gene" and to_type == "transcript":
        return _gpos_to_tpos
    elif from_type == "protein" and to_type == "cds":
        return _ppos_to_cpos
    elif from_type == "transcript" and to_type == "cds":
        return _tpos_to_cpos
    elif from_type == "transcript" and to_type == "gene":
        return _tpos_to_gpos
    else:
        raise ValueError(f"Cannot convert {from_type} directly to {to_type}")


def _cpos_to_ppos(_, start, end=None):
    """Compute the equivalent protein position for a CDS position.

    Args:
        start (int): position relative to the CDS
        end (int): optional, second position relative to the CDS

    Returns:
        tuple of int: amino acid position
    """
    pstart = floor((start - 1) / 3 + 1)
    if end:
        pend = _cpos_to_ppos(_, end)[1]
    else:
        pend = pstart

    return pstart, pend


def _cpos_to_tpos(transcript, start, end=None):
    """Compute the equivalent CDS position for a transcript position.
    
    Args:
        transcript: `pyensembl.Transcript` instance
        start (int): position relative to the CDS
        end (int): optional, second position relative to the CDS

    Returns:
        tuple of int: position relative to the transcript
    """
    tstart = transcript.first_start_codon_spliced_offset + start - 1
    if not (1 <= tstart <= len(transcript.sequence)):
        raise TranscriptOutOfRange(transcript, start)
    if end:
        tend = _cpos_to_tpos(transcript, end)[1]
    else:
        tend = tstart

    return tstart, tend


def _epos_to_cpos(transcript, start, end):
    """Compute the equivalent CDS position for an exon.
    
    Args:
        transcript: `pyensembl.Transcript` instance
        start (int): genomic coordinate of one end of the exon
        end (int): genomic coordinate of the other end of the exon

    Returns:
        tuple of int: position relative to the CDS
    """
    for cds in transcript.coding_sequence_position_ranges:
        if start == cds[0] or end == cds[1]:
            gstart, gend = cds
            break
    else:
        raise CdsOutOfRange(transcript, start)

    tstart, tend = _gpos_to_tpos(transcript, gstart, gend)
    cstart, cend = _tpos_to_cpos(transcript, tstart, tend)
    return cstart, cend


def _gpos_to_tpos(transcript, start, end=None):
    """Compute the equivalent transcript position for a gene position.
    
    Args:
        transcript: `pyensembl.Transcript` instance
        start (int): genomic coordinate
        end (int): optional, second genomic coordinate

    Returns:
        tuple of int: position relative to the transcript
    """
    tstart = transcript.spliced_offset(start) + 1
    if end:
        tend = _gpos_to_tpos(transcript, end)[1]
    else:
        tend = tstart

    return tstart, tend


def _ppos_to_cpos(_, start, end=None):
    """Compute the equivalent protein position for a CDS position.

    Args:
        start (int): amino acid position
        end (int): optional, second amino acid position

    Returns:
        tuple of int: CDS position of the first base of the codon
    """
    cstart = (start - 1) * 3 + 1
    # add 2 to the end position to get the end of the codon
    if end:
        cend = _ppos_to_cpos(_, end)[1]
    else:
        cend = cstart + 2

    return cstart, cend


def _gpos_to_epos(transcript, start, end=None):
    """Return the genomic coordinates of the exon that contains a genomic coordinate.
    
    Args:
        transcript: `pyensembl.Transcript` instance
        start (int): genomic coordinate
        end (int): optional, second genomic coordinate

    Returns:
        tuple of int: genomic coordinates of the exon
    """
    for exon in sorted(transcript.exons, key=lambda x: x.start):
        if exon.start <= start <= exon.end:
            exon1 = (exon.start, exon.end, exon.exon_id)
            break
    else:
        raise ExonOutOfRange(transcript, start)

    if end:
        exon2 = _gpos_to_epos(transcript, end)
    else:
        exon2 = exon1

    if exon1 != exon2:
        raise ValueError(f"{start} and {end} are on different exons ({exon1[2]}, {exon2[2]}")
    else:
        return exon1


def _tpos_to_cpos(transcript, start, end=None):
    """Compute the equivalent transcript position for a CDS position.
    
    Args:
        transcript: `pyensembl.Transcript` instance
        start (int): position relative to the transcript
        end (int): optional, second position relative to the transcript

    Returns:
        tuple of int: position relative to the CDS
    """
    cstart = start - transcript.first_start_codon_spliced_offset
    if not (1 <= cstart <= len(transcript.coding_sequence)):
        raise CdsOutOfRange(transcript, cstart)

    if end:
        cend = _tpos_to_cpos(transcript, end)[1]
    else:
        cend = cstart

    return cstart, cend


def _tpos_to_gpos(transcript, start, end=None):
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

    remain = start - 1
    for i in ranges:
        length = i[1] - i[0] + 1
        if remain >= length:
            remain -= length
        else:
            if transcript.on_positive_strand:
                gstart = i[0] + remain
                break
            elif transcript.on_negative_strand:
                gstart = i[1] - remain
                break
    else:
        raise TranscriptOutOfRange(transcript, start)

    if end:
        gend = _tpos_to_gpos(transcript, end)[1]
    else:
        gend = gstart

    return gstart, gend
