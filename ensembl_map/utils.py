from itertools import product, zip_longest
from typing import Iterator, List, Optional, Tuple

from Bio.Seq import Seq

from .tables import DNA, DNA_CODON_TABLE, PROTEIN


def collapse_mutation(ref: str, alt: str) -> Tuple[str, str, int, int]:
    """Collapse a nucleotide change to the shortest change.

    Examples:
                  ref:  100 - GGT - 102
                  alt:  100 - GTG - 102
                               ^^
        --------------------------------
        collapsed alt:  101 -  TG - 102
    """
    ref_collapse = ref
    alt_collapse = alt
    start_offset = 0
    end_offset = 0

    # Special case when the ref and alt are the same
    if ref == alt:
        return ref, alt, 0, 0

    # If the event is a duplication, return the bases as is
    if is_duplication(ref, alt):
        return ref, alt, 0, 0

    # Trim bases from each end that are identical between the ref and alt
    ref_collapse, alt_collapse, same_5_prime, same_3_prime = collapse_mutation_ends(ref, alt)
    start_offset = len(same_5_prime)
    end_offset = len(same_3_prime)

    # If the event is an insertion or deletion, preserve the 5' and 3' ref bases
    if same_5_prime and same_3_prime:
        if len(ref_collapse) < len(alt_collapse):
            ref_collapse = same_5_prime[-1] + ref_collapse + same_3_prime[0]
            alt_collapse = same_5_prime[-1] + alt_collapse + same_3_prime[0]
            start_offset -= 1
            end_offset -= 1

    # If the ref is collapsed to an empty string it generally indicates an amino acid duplication
    if not ref_collapse:
        return ref, alt, 0, 0

    return ref_collapse, alt_collapse, start_offset, end_offset


def collapse_mutation_ends(refseq: str, altseq: str) -> Tuple[str, str, str, str]:
    """Compare two strings (ref and alt) and split each at the point where the two sequences are no
    longer the same when read from one end or the other (determined by `idx` being either 0 or -1).

    For example, 'AATTTC' and 'AAGC' both start with 'AA' and end with 'C'

    Examples:
        >>> _split_at_common_substring('AATTTC', 'AAGC') == ('TTT', 'G', 'AA', 'C')
    """
    common_left: List[str] = list()
    common_right: List[str] = list()
    refseq_ = list(refseq)
    altseq_ = list(altseq)

    def trim_one_side(idx: int, common: List):
        assert idx in (0, -1)

        while True:
            try:
                i = refseq_.pop(idx)
            except IndexError:
                i = ""

            try:
                j = altseq_.pop(idx)
            except IndexError:
                j = ""

            if i != j or not i or not j:
                # Put the 'popped' bases back then stop
                if i:
                    if idx == 0:
                        refseq_.insert(0, i)
                    elif idx == -1:
                        refseq_.append(i)
                if j:
                    if idx == 0:
                        altseq_.insert(0, j)
                    elif idx == -1:
                        altseq_.append(j)
                break
            else:
                if idx == 0:
                    common.append(i)
                elif idx == -1:
                    common.insert(0, i)

    trim_one_side(0, common_left)  # Trim 5' (left) end
    trim_one_side(-1, common_right)  # Trim 3' (right) end

    return "".join(refseq_), "".join(altseq_), "".join(common_left), "".join(common_right)


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


def is_deletion(refseq: str, altseq: str) -> bool:
    """Check if an allele change should be classified as a deletion mutation."""
    return len(refseq) > 0 and len(altseq) == 0


def is_duplication(refseq: str, altseq: str) -> bool:
    """Check if an allele change should be classified as a duplication mutation."""
    return altseq == refseq * 2


def is_frameshift(cdna_refseq: str, cdna_altseq: str) -> bool:
    """Check if a cDNA nucleotide change would result in a frameshift mutation."""
    return abs(len(cdna_refseq) - len(cdna_altseq)) % 3 != 0


def is_insertion(refseq: str, altseq: str) -> bool:
    """Check if an allele change should be classified as an insertion mutation."""
    return bool(split_insertion(refseq, altseq))


def is_substitution(refseq: str, altseq: str) -> bool:
    """Check if an allele change should be classified as a substitution mutation."""
    return len(refseq) == 1 and len(altseq) == 1


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


def split_insertion(refseq: str, altseq: str) -> Optional[Tuple[str, str, str]]:
    """Find an insertion (if any). This is done by finding a substring such that if the substring
    was removed from the alt, the alt would be the same as the ref.
    """
    # If the alt does not have more bases than the ref, it can't be an insertion
    if len(altseq) <= len(refseq):
        return None

    # Split the ref at different positions, if both halves match the ends of the alt, it's an insertion
    # e.g. 'GGGTTT' -> ['GGG', 'TTT']
    # Try in the middle first, since that's the most likely position
    mid = len(refseq) // 2
    idx = [mid] + [i for i in range(1, len(refseq)) if i != mid]
    for n in idx:
        l, r = refseq[:n], refseq[n:]
        if altseq[:n] == l and altseq[-n:] == r:
            return (altseq[:n], altseq[n:-n], altseq[-n:])

    return None


def strip_version(key: str) -> str:
    """Strip the version number from a transcript symbol.

    Examples:
        >>> _strip_version('NM_000546.5')
        'NM_000546'
        >>> _strip_version('ENST00000357191.1')
        'ENST00000357191'
    """
    return key.rsplit(".", 1)[0]
