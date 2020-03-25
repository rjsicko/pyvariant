from math import floor

from .exceptions import OutOfRangeError


def get_map_function(from_type, to_type):
    if from_type == to_type:
        raise ValueError("`from_type` and `to_type` must be different!")
    elif from_type == "cds" and to_type == "protein":
        return _cpos_to_ppos
    elif from_type == "cds" and to_type == "transcript":
        return _cpos_to_tpos
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


def _cpos_to_ppos(_, position):
    """Compute the equivalent protein position for a CDS position.

    Args:
        position (int): position relative to the CDS

    Returns:
        int: amino acid position
    """
    return floor((position - 1) / 3 + 1)


def _cpos_to_tpos(transcript, position):
    """Compute the equivalent CDS position for a transcript position.
    
    Args:
        transcript: `pyensembl.Transcript` instance
        position (int): position relative to the CDS

    Returns:
        int: position relative to the transcript
    """
    tpos = transcript.first_start_codon_spliced_offset + position - 1
    if not (1 <= tpos <= len(transcript.sequence)):
        raise OutOfRangeError(transcript, position, "transcript")
    return tpos


def _gpos_to_tpos(transcript, position):
    """Compute the equivalent transcript position for a gene position.
    
    Args:
        transcript: `pyensembl.Transcript` instance
        position (int): genomic coordinate

    Returns:
        int: position relative to the transcript
    """
    return transcript.spliced_offset(position) + 1


def _ppos_to_cpos(_, position):
    """Compute the equivalent protein position for a CDS position.

    Args:
        position (int): amino acid position

    Returns:
        int: CDS position of the first base of the codon
    """
    return (position - 1) * 3 + 1


def _gpos_to_epos(transcript, position):
    """Return the genomic coordinates of the exon that contains a genomic coordinate.
    
    Args:
        transcript: `pyensembl.Transcript` instance
        position (int): genomic coordinate

    Returns:
        tuple of int: genomic coordinates of the exon
    """
    for exon in sorted(transcript.exons, key=lambda x: x.start):
        if exon.start <= position <= exon.end:
            return exon.start, exon.end, exon.exon_id
    else:
        raise AssertionError("Iterated through all exons")


def _tpos_to_cpos(transcript, position):
    """Compute the equivalent transcript position for a CDS position.
    
    Args:
        transcript: `pyensembl.Transcript` instance
        position (int): position relative to the transcript

    Returns:
        cpos (int): position relative to the CDS
    """
    cpos = position - transcript.first_start_codon_spliced_offset
    if not (1 <= cpos <= len(transcript.coding_sequence)):
        raise OutOfRangeError(transcript, cpos, "CDS")
    return cpos


def _tpos_to_gpos(transcript, position):
    """Compute the equivalent gene position for a transcript position.
    
    Args:
        transcript: `pyensembl.Transcript` instance
        position (int): position relative to the transcript

    Returns:
        int: genomic coordinate
    """
    # make sure all the ranges are sorted from smallest to biggest
    ranges = sorted([sorted(i) for i in transcript.exon_intervals])

    # for transcripts on the negative strand, start counting from the last position
    if transcript.on_negative_strand:
        ranges = ranges[::-1]

    remain = position - 1
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
        raise OutOfRangeError(transcript, position, "transcript")
