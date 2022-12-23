from itertools import product, zip_longest
from typing import Iterator, Tuple

from Bio.Seq import Seq

from .tables import DNA, DNA_CODON_TABLE, PROTEIN


def collapse_mutation(ref: str, alt: str) -> Tuple[str, str, int, int]:
    """Collapse a codon change to the shortest nucleotide change.

    Examples:
                  ref:  100 - GGT - 102
                  alt:  100 - GTG - 102
                               ^^
        --------------------------------
        collapsed alt:  101 -  TG - 102

        >>> collapse_mutation("GGT", "GTG")
        ("GT", "TG", 1, 0)
        >>> collapse_mutation("GGT", "GTT")
        ("G", "T", 1, 1)
        >>> collapse_mutation("GTT", "CTT")
        ("G", "C", 0, 2)
        >>> collapse_mutation("TTG", "TTA")
        ("G", "A", 2, 0)
        >>> collapse_mutation("ATG", "A")
        ("ATG", "A", 0, 0)
        >>> collapse_mutation("ATG", "ATG")
        ("ATG", "ATG", 0, 0)
        >>> collapse_mutation("ATG", "ATGATG")
        ("ATG", "ATGATG", 0, 0)
        >>> collapse_mutation("ATGATG", "ATG")
        ("ATGATG", "ATG", 0, 0)
    """
    # Special case when the ref and alt are the same
    if ref == alt:
        return ref, alt, 0, 0

    # Can't collapse if codons are different lengths
    if len(ref) != len(alt):
        return ref, alt, 0, 0

    alt_diff = [i != ref[n] for n, i in enumerate(alt)]  # "GCA" and "GTG" = [False, True, True]
    alt_index = [n for n, i in enumerate(alt_diff) if i]  # [False, True, True] = [1, 2]

    if alt_index:
        start_index = alt_index[0]  # [1, 2][0] = [1]
        end_index = alt_index[-1]  # [1, 2][-1] = [2]
    else:
        start_index = 0
        end_index = len(ref) - 1

    ref_collapsed = ref[start_index : end_index + 1]  # "GCA"[1:3] = CA
    alt_collapsed = alt[start_index : end_index + 1]  # "GTG"[1:3] = TG

    start_offset = start_index
    end_offset = len(ref) - 1 - end_index  # 3 - 1 - 2 = 0

    return ref_collapsed, alt_collapsed, start_offset, end_offset


def expand_nt(seq: str) -> Iterator[str]:
    """Expand any ambiguous nucleotides (e.g. 'ARA' -> 'AAA', 'AGA')."""
    yield from iter("".join(i) for i in product(*("".join(DNA[nt]) for nt in seq)))


def expand_pep(seq: str) -> Iterator[str]:
    """Expand any ambiguous peptides (e.g. 'B' -> 'D', 'N')."""
    yield from iter("".join(i) for i in product(*("".join(PROTEIN[aa]) for aa in seq)))


def format_hgvs_position(position: int, offset: int, is_3_prime_utr: bool = False) -> str:
    """Format a position as a string according to HGVS standards."""
    position_str = ""

    if is_3_prime_utr:
        position_str = "*"

    if position:
        if offset and position == 1 and not is_3_prime_utr:
            position_str += f"{offset:+}"
        elif offset:
            position_str += f"{position}{offset:+}"
        else:
            position_str += f"{position}"

    return position_str


def is_frameshift(cdna_start: int, cdna_end: int, cdna_altseq: str) -> bool:
    """Check if a cDNA nucleotide change would result in a frameshift mutation."""
    reflen = cdna_end - cdna_start + 1
    altlen = len(cdna_altseq)

    return abs(reflen - altlen) % 3 != 0


def reverse_complement(sequence: str) -> str:
    """Return the reverse complement of a given nucleotide sequence."""
    return str(Seq(sequence).reverse_complement())


def reverse_translate(peptide: str) -> Iterator[str]:
    """Return all reverse translations of a peptide (e.g. 'C' -> 'TGC', 'TGC', 'TGT', 'TGT')."""
    for pep in expand_pep(peptide):
        yield from iter("".join(i) for i in product(*[DNA_CODON_TABLE[aa] for aa in pep]))


def split_by_codon(iterable: str) -> Iterator[str]:
    """Collect a sequence into chunks of 3 (adapted from the 'grouper' recipe found
    here: https://docs.python.org/3.8/library/itertools.html#itertools-recipes).

    Examples:
        >>> split_by_codon('ABCDEF')
        ['ABC', 'DEF']
    """
    if len(iterable) % 3 != 0:
        raise ValueError(f"Iterable ({iterable}) is not divisible by 3")

    args = [iter(iterable)] * 3
    yield from iter("".join(i) for i in zip_longest(*args))


def strip_version(key: str) -> str:
    """Strip the version number from a transcript symbol.

    Examples:
        >>> _strip_version('NM_000546.5')
        'NM_000546'
        >>> _strip_version('ENST00000357191.1')
        'ENST00000357191'
    """
    return key.rsplit(".", 1)[0]
