"""Function definitions for manipulating reference sequences."""
from pyfaidx import Fasta


def get_sequence(fasta: Fasta, ref: str, start: int, end: int, window: int) -> str:
    """Return the reference sequence at the given position plus additional bases 5' and 3' such
    that the total returned sequence length adds up to `window`.

    Args:
        fasta (Fasta): FASTA containing all the reference sequences
        ref (str): Reference sequence ID
        start (int): First base of position of interest
        end (int): Last base of position of interest
        window (int): Total length of the returned sequence

    Returns:
        sequence (str): Reference sequence with flanking 5' and 3' bases
    """
    # `start` and `end` are 1-based positions but we need 0-based positions to index `fasta`
    start = start - 1

    # Determine the number of bases flanking the position site that need to be retrieved
    pos_len = end - start
    if not window or window < 1:
        window = pos_len
    elif window < pos_len:
        raise ValueError(f"sequence window is smaller than position length ({window} < {pos_len})")

    window -= pos_len

    # We want the mutation site centered in the returned sequence. Ideally, there would be an
    # equal number of bases from the reference sequence flanking the mutation site. If there can't
    # be, the 3' end will have 1 more base than the 5' end.
    pad_left = window // 2
    pad_right = window - pad_left

    # reference start = position start - 1/2 the desired length (rounded down)
    ref_start = start - pad_left

    # reference end = position end - 1/2 the desired length (rounded up)
    ref_end = end + pad_right

    # Get enough of the reference sequence that we can return a sequence of `length`
    sequence = fasta[ref][ref_start:ref_end]

    return sequence


def mutate_sequence(
    fasta: Fasta, ref: str, start: int, end: int, window: int, altseq: str, insertion: bool = False
) -> str:
    """Return a sequence that represents the alternate allele in the context of the reference
    sequence within a window of length `window`.

    Args:
        fasta (Fasta): FASTA containing all the reference sequences
        ref (str): Reference sequence ID
        start (int): First base of the mutation site. For insertions, this is the base 5' of the
            insertion site
        end (int): Last base of the mutation site. For insertions, this is the base 3' of the
            insertion site
        window (int): Total length of the returned sequence
        altseq (str): Alternate allele
        insertion (bool, optional): True if the mutation is an insertion. Defaults to False.

    Returns:
        sequence (str): Alternate allele with flanking 5' and 3' bases
    """
    # If a window is not give, we can just return the alternate sequence
    if not window or window < 0:
        return altseq

    # Catch cases where the given window size is less than the size of the alternate sequence
    if window < len(altseq):
        raise ValueError(
            f"sequence window is smaller than variant length ({window} < {len(altseq)})"
        )

    # If the window size is an odd number, pad the 3' more than the 5' end
    pad_left = window // 2
    pad_right = window - pad_left

    # Calculate the reference start and end positions required to return a sequence of size `window`
    del_len = end - start + 1
    ins_len_left = len(altseq) // 2
    ins_len_right = len(altseq) - ins_len_left
    ins_pad = 1 if insertion else 0

    ref_start = start - pad_left + ins_len_left + ins_pad - 1
    assert ref_start <= start, f"{ref_start} > {start}"

    ref_end = start + pad_right - ins_len_right + del_len - ins_pad - 1
    assert ref_end >= ref_start, f"{ref_end} < {ref_start}"

    refseq = fasta[ref][ref_start:ref_end]

    # Piece together the altseq and 5'/3' flanking sequences
    idx_left = start - ref_start + ins_pad - 1
    left = refseq[:idx_left]

    idx_right = ref_end - end + ins_pad
    right = refseq[-idx_right:]

    sequence = left + altseq + right

    return sequence
