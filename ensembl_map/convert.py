from math import floor

from .exceptions import CdsRangeError, ExonRangeError
from .util import in_ranges


def cpos_to_gpos(tobj, position):
    """
    Compute the equivalent gene position for a CDS position.

    Args:
        tobj: Transcript
        position (int): position relative to the CDS

    Returns:
        (int): position relative to the contig

    Examples:
        >>> tobj.coding_sequence_position_ranges = [(11, 15), (17, 21), (23, 27)]
        >>> cpos_to_gpos(tobj, 1)
        11
        >>> cpos_to_gpos(tobj, 7)
        18
        >>> cpos_to_gpos(tobj, 1)
        27
        >>> cpos_to_gpos(tobj, 6)
        21
        >>> cpos_to_gpos(tobj, 99)
        Traceback (most recent call last):
            ...
        AssertionError: 99 is outside the CDS
    """
    cds_ranges = sorted(tobj.coding_sequence_position_ranges)
    if tobj.on_negative_strand:
        cds_ranges = cds_ranges[::-1]

    if not in_ranges(position, cds_ranges):
        raise CdsRangeError(tobj, position)

    remain = position - 1
    for cds in cds_ranges:
        length = cds[1] - cds[0] + 1
        if remain >= length:
            remain -= length
        else:
            if tobj.on_negative_strand:
                return cds[0] + remain
            else:
                return cds[1] - remain
    else:
        raise AssertionError(f"Iterated through all ranges")


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
    raise NotImplementedError  # TODO


def epos_to_cpos(tobj, position):
    """Compute the equivalent exon position for a CDS position."""
    raise NotImplementedError  # TODO


def epos_to_tpos(tobj, position):
    """Compute the equivalent exon position for a transcript position."""
    raise NotImplementedError  # TODO


def gpos_to_cpos(tobj, position):
    """
    Compute the equivalent CDS position for a gene position.
    
    Args:
        tobj: Transcript
        position (int): position relative to the positive strand of the contig

    Returns:
        (int): position relative to the contig

    Examples:
        >>> tobj.coding_sequence_position_ranges = [(11, 15), (17, 21), (23, 27)]
        >>> dna_to_cds(tobj, 11)
        1
        >>> dna_to_cds(tobj, 27)
        15
        >>> dna_to_cds(tobj, 20)
        7
        >>> dna_to_cds(tobj, 999999999)
        Traceback (most recent call last):
            ...
        AssertionError: 99 is outside the CDS
    """
    cds_ranges = sorted(tobj.coding_sequence_position_ranges)
    if tobj.on_negative_strand:
        cds_ranges = cds_ranges[::-1]

    if not in_ranges(position, cds_ranges):
        raise CdsRangeError(tobj, position)

    cpos = 0
    for cds in cds_ranges:
        # iterate until the postion falls inside one of the ranges
        if cds[0] <= position <= cds[1]:
            if tobj.on_positive_strand:
                cpos += position - cds[0] + 1
            else:
                cpos += cds[1] - position + 1
            break
        cpos += cds[1] - cds[0] + 1
    else:
        raise AssertionError(f"Iterated through all ranges")

    return cpos


def gpos_to_tpos(tobj, position):
    """Compute the equivalent transcript position for a gene position."""
    raise NotImplementedError  # TODO


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


def ppos_to_tpos(tobj, position):
    """Compute the equivalent protein position for a transcript position."""
    raise NotImplementedError  # TODO


def tpos_to_epos(tobj, position):
    """Compute the equivalent exon position for a transcript position."""
    for n, exon in enumerate(sorted(tobj.exons), 1):
        if exon.start <= position <= exon.end:
            return n
    else:
        raise ExonRangeError(tobj, position)


def tpos_to_cpos(tobj, position):
    """Compute the equivalent transcript position for a CDS position."""
    raise NotImplementedError  # TODO


def tpos_to_gpos(tobj, position):
    """Compute the equivalent gene position for a transcript position."""
    raise NotImplementedError  # TODO


def tpos_to_ppos(tobj, position):
    """Compute the equivalent protein position for a transcript position."""
    raise NotImplementedError  # TODO
