from math import floor

from .exceptions import OutsideCdsError, OutsideTranscriptError


def cpos_to_ppos(_, position):
    """Compute the equivalent protein position for a CDS position.

    Args:
        position (int): position relative to the CDS

    Returns:
        int: amino acid position
    """
    return floor((position - 1) / 3 + 1)


def cpos_to_tpos(tobj, position):
    """Compute the equivalent CDS position for a transcript position.
    
    Args:
        tobj: `Transcript` instance
        position (int): position relative to the CDS

    Returns:
        int: position relative to the transcript
    """
    tpos = tobj.first_start_codon_spliced_offset + position - 1
    if not (1 <= tpos <= len(tobj.sequence)):
        raise OutsideTranscriptError(tobj, position)
    return tpos


def epos_to_gpos(tobj, position):
    """Return the genomic coordinates of the nth exon of the given transcript.
    
    Args:
        tobj: `Transcript` instance
        position (int): index of the exon relative to the gene (i.e. the nth exon)

    Returns:
        tuple of int: genomic coordinates of the exon
    """
    try:
        return sorted(tobj.exon_intervals)[position - 1]
    except IndexError:
        raise ValueError(f"{position} is greater than the number of exons")


def gpos_to_tpos(tobj, position):
    """Compute the equivalent transcript position for a gene position.
    
    Args:
        tobj: `Transcript` instance
        position (int): genomic coordinate

    Returns:
        int: position relative to the transcript
    """
    return tobj.spliced_offset(position) + 1


def ppos_to_cpos(_, position):
    """Compute the equivalent protein position for a CDS position.

    Args:
        position (int): amino acid position relative

    Returns:
        int: CDS position of the first base of the codon
    """
    return (position - 1) * 3 + 1


def gpos_to_epos(tobj, position):
    """Return the genomic coordinates of the exon that contains a genomic coordinate.
    
    Args:
        tobj: `Transcript` instance
        position (int): genomic coordinate

    Returns:
        tuple of int: genomic coordinates of the exon
    """
    for exon in sorted(tobj.exon_intervals):
        if exon[0] <= position <= exon[1]:
            return exon
    else:
        raise AssertionError("Iterated through all exons")


def tpos_to_cpos(tobj, position):
    """Compute the equivalent transcript position for a CDS position.
    
    Args:
        tobj: `Transcript` instance
        position (int): position relative to the transcript

    Returns:
        cpos (int): position relative to the CDS
    """
    cpos = position - tobj.first_start_codon_spliced_offset + 1
    if not (1 <= cpos <= len(tobj.coding_sequence)):
        raise OutsideCdsError(tobj, cpos)
    return cpos


def tpos_to_gpos(tobj, position):
    """Compute the equivalent gene position for a transcript position.
    
    Args:
        tobj: `Transcript` instance
        position (int): position relative to the transcript

    Returns:
        int: genomic coordinate
    """
    # make sure all the ranges are sorted from smallest to biggest
    ranges = sorted([sorted(i) for i in tobj.exon_intervals])

    # for transcripts on the negative strand, start counting from the last position
    if tobj.on_negative_strand:
        ranges = ranges[::-1]

    remain = position - 1
    for i in ranges:
        length = i[1] - i[0] + 1
        if remain > length:
            remain -= length
        else:
            if tobj.on_positive_strand:
                return i[0] + remain
            elif tobj.on_negative_strand:
                return i[1] - remain
    else:
        raise OutsideTranscriptError(tobj, position)
