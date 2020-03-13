from math import floor

from .exceptions import CdsRangeError, ExonRangeError
from .util import in_ranges


def cpos_to_ppos(_, position):
    """
    Compute the equivalent protein position for a CDS position.

    Args:
        position (int): position relative to the start of the CDS

    Returns:
        (int): amino acid position

    Examples:
        >>> cpos_to_ppos(_, 1)
        1
        >>> cpos_to_ppos(_, 4)
        2
        >>> cpos_to_ppos(_, 5)
        2
        >>> cpos_to_ppos(_, 7)
        3
    """
    return floor((position - 1) / 3 + 1)


def cpos_to_tpos(tobj, position):
    """Compute the equivalent CDS position for a transcript position."""
    return tobj.first_start_codon_spliced_offset + position - 1


def epos_to_tpos(tobj, position):
    """Compute the equivalent exon position for a transcript position."""
    return sorted(tobj.exon_intervals)[position - 1]


def gpos_to_tpos(tobj, position):
    """Compute the equivalent transcript position for a gene position."""
    return tobj.spliced_offset(position) + 1


def ppos_to_cpos(_, position):
    """
    Compute the equivalent protein position for a CDS position.

    Args:
        position (int): amino acid position relative to the peptide

    Returns:
        (int): CDS position of the first base of the codon

    Examples:
        >>> ppos_to_cpos(1)
        1
        >>> ppos_to_cpos(2)
        4
    """
    return (position - 1) * 3 + 1


def tpos_to_epos(tobj, position):
    """Compute the equivalent exon position for a transcript position."""
    for n, exon in enumerate(sorted(tobj.exons), 1):
        if exon.start <= position <= exon.end:
            return n
    else:
        raise ExonRangeError(tobj, position)


def tpos_to_cpos(tobj, position):
    """Compute the equivalent transcript position for a CDS position."""
    return tobj.first_start_codon_spliced_offset - position + 1


def tpos_to_gpos(tobj, position):
    """Compute the equivalent gene position for a transcript position."""
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
            if tobj.on_negative_strand:
                return i[0] + remain
            else:
                return i[1] - remain
    else:
        raise AssertionError(f"Iterated through all ranges")
