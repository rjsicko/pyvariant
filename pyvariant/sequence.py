"""Function definitions for manipulating reference sequences."""
from __future__ import annotations

from functools import lru_cache
from typing import Optional

from pyfaidx import Fasta

from .constants import EMPTY_FASTA
from .utils import reverse_complement, strip_version


class PyfaidxFasta:
    read_ahead = 100000

    @classmethod
    def load(cls, path: str = "", strand: str = "") -> PyfaidxFasta:
        """Parse a FASTA/FASTQ file into a 'pyfaidx.Fasta' object.

        Args:
            path (str): Path to FASTA or FASTQ file

        Returns:
            PyfaidxFasta
        """
        if not path:
            path = EMPTY_FASTA

        fasta = Fasta(
            path,
            key_function=strip_version,
            as_raw=True,
            sequence_always_upper=True,
            build_index=False,
            rebuild=False,
            read_ahead=cls.read_ahead,
        )

        return cls(fasta, strand=strand)

    def __init__(self, fasta: Fasta, strand: str = ""):
        self.fasta = fasta
        self.strand = strand

    def __contains__(self, reference: str) -> bool:
        return reference in self.fasta

    def __getitem__(self, reference: str) -> str:
        return self.fasta[reference]

    def _try_buffer(self, reference: str, start: int, end: int):
        """This is a bit of a hack to buffer more of the sequence in memory. Normally, pyfaidx only
        buffers N bases *ahead* of the given position. This hack also buffers N bases *behind*. This
        is useful because we commonly need to fetch multiple times around a given position and that
        speeds the queries up.

        Args:
            reference (str): Reference sequence ID
            start (int): First base of position of interest
            end (int): Last base of position of interest
        """
        try:
            if reference != self.fasta.faidx.buffer["name"]:
                floor = max(start - self.read_ahead, 0)
                ceiling = min(end + self.read_ahead, self.length(reference) - 1)
                self.fasta[reference][floor:ceiling]
        except Exception:
            pass

    def fetch(self, reference: str, start: int, end: int) -> str:
        self._try_buffer(reference, start, end)
        return str(self.fasta[reference][start:end])

    def length(self, reference: str) -> int:
        return len(self.fasta[reference])


@lru_cache
def get_sequence(
    fasta: PyfaidxFasta,
    ref: str,
    start: int,
    end: int,
    strand: Optional[str],
    window: Optional[int],
    floor: Optional[int],
    ceiling: Optional[int],
) -> str:
    """Return the reference sequence at the given position plus additional bases 5' and 3' such
    that the total returned sequence length adds up to `window`.

    Args:
        fasta (PyfaidxFasta): FASTA containing all the reference sequences
        ref (str): Reference sequence ID
        start (int): First base of position of interest
        end (int): Last base of position of interest
        strand (Optional[str]): Strand relative to the given position (not the strand of `fasta`)
        window (Optional[int]): Total length of the returned sequence
        floor (Optional[int]): Absolute smallest start position (must be <= `start`)
        ceiling (Optional[int]): Absolute largest end position (must be >= `end`)

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

    # Adjust the reference start if it extends past the start of the molecule (i.e. can't be < 1)
    if ref_start < 0:
        ref_end -= ref_start
        ref_start = 0

    # If the floor is more than ref_start, make the new ref_start equal to floor and pad ref_end by
    # the difference
    floor = floor - 1 if (floor and floor > 0) else ref_start
    diff_left = floor - ref_start

    if ceiling and ceiling > 0:
        diff_right = ceiling - ref_end
        if diff_left > 0 and diff_left <= diff_right:
            ref_start -= diff_left
            ref_end -= diff_left
        elif diff_right < 0 and diff_right <= diff_left:
            ref_start += diff_right
            ref_end += diff_right
        else:
            raise ValueError("window constraints are smaller than the window length")
    elif diff_left > 0:
        ref_start += diff_left
        ref_end += diff_left

    # Adjust the reference start and end if it extends past the end of the molecule
    max_end = fasta.length(ref)
    if ref_end > max_end:
        diff = min(ref_end - max_end, ref_start)
        ref_end = max_end
        ref_start -= diff

    # Get enough of the reference sequence that we can return a sequence of `length`
    sequence = fasta.fetch(ref, ref_start, ref_end)

    # Reverse complement the sequence if the position strand isn't the same as the FASTA strand
    if (strand and fasta.strand) and (strand != fasta.strand):
        sequence = reverse_complement(sequence)

    return sequence


@lru_cache
def mutate_sequence(
    fasta: PyfaidxFasta,
    ref: str,
    start: int,
    end: int,
    strand: Optional[str],
    window: Optional[int],
    floor: Optional[int],
    ceiling: Optional[int],
    altseq: str,
    insertion: bool = False,
) -> str:
    """Return a sequence that represents the alternate allele in the context of the reference
    sequence within a window of length `window`.

    Args:
        fasta (PyfaidxFasta): FASTA containing all the reference sequences
        ref (str): Reference sequence ID
        start (int): First base of the mutation site. For insertions, this is the base 5' of the
            insertion site
        end (int): Last base of the mutation site. For insertions, this is the base 3' of the
            insertion site
        strand (Optional[str]): Strand relative to the given position (not the strand of `fasta`)
        window (Optional[int]): Total length of the returned sequence
        floor (Optional[int]): Absolute smallest start position (must be <= `start`)
        ceiling (Optional[int]): Absolute largest end position (must be >= `end`)
        altseq (str): Alternate allele. For insertions, requires the two bases immediately 5' and
            3' of the insertion site
        insertion (bool, optional): DEPRECATED

    Returns:
        sequence (str): Alternate allele with flanking 5' and 3' bases
    """
    # If a window is not given, we can just return the alternate sequence
    if not window or window < 0:
        return altseq

    # Validate that a reference name was given
    if not ref:
        raise ValueError(f"empty reference name {ref=}")

    # Validate the upper and lower bounds
    if not floor or floor < 1:
        floor = 1

    if not ceiling or ceiling < 1:
        ceiling = fasta.length(ref)

    if floor > ceiling:
        raise ValueError(f"{floor=} > {ceiling=}")

    # Validate that the given start/end are within the upper and lower bounds
    if not (floor <= start <= ceiling):
        raise ValueError(f"{start=} is out of range ({floor}, {ceiling})")

    if not (floor <= end <= ceiling):
        raise ValueError(f"{end=} is out of range ({floor}, {ceiling})")

    if start > end:
        raise ValueError(f"{start=} > {end=}")

    # Catch cases where the given window size is less than the size of the alternate sequence
    if window < len(altseq):
        raise ValueError(
            f"sequence window is smaller than variant length ({window} < {len(altseq)})"
        )

    # Decrease the effective window size by the number of bases inserted
    window -= len(altseq)

    # If the window size is an odd number, pad the 3' more than the 5' end
    pad_left = window // 2
    pad_right = window - pad_left

    # Get the positions 5' and 3' of the deletion/insertion site within the window
    left_start = start - pad_left
    left_end = start - 1
    right_start = end + 1
    right_end = end + pad_right

    # Shift the window 5' or 3' if it would exceed the upper or lower bound
    shift_right = max(floor - left_start, 0)
    shift_left = max(right_end - ceiling, 0)
    left_start = max(left_start - shift_left, floor)
    # left_end = max(left_end, floor)
    # right_start = min(right_start, ceiling)
    right_end = min(right_end + shift_right, ceiling)

    # Reverse complement the sequence if the position strand isn't the same as the FASTA strand
    opposite_strand = (strand and fasta.strand) and (strand != fasta.strand)
    if opposite_strand:
        altseq = reverse_complement(altseq)

    # Fetch the sequence and piece together the altseq and 5'/3' flanking sequences
    # Positions use 0-based numbering
    if left_start <= left_end:
        left = fasta.fetch(ref, left_start - 1, left_end)
    else:
        left = ""

    if right_start <= right_end:
        right = fasta.fetch(ref, right_start - 1, right_end)
    else:
        right = ""

    sequence = left + altseq + right

    if opposite_strand:
        sequence = reverse_complement(sequence)

    return sequence
