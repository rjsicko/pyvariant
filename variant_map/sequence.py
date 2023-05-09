"""Function definitions for manipulating reference sequences."""
from __future__ import annotations

import gzip
from abc import ABC
from functools import lru_cache
from typing import Dict, Optional

from Bio import SeqIO
from pyfaidx import Fasta

from .constants import EMPTY_FASTA
from .files import is_bgzipped
from .utils import strip_version


class _FastaABC(ABC):
    def __contains__(self, reference: str) -> bool:
        raise NotImplementedError()

    def __getitem__(self, reference: str) -> str:
        raise NotImplementedError()

    def fetch(self, reference: str, start: int, end: int) -> str:
        raise NotImplementedError()


class DictFasta(_FastaABC):
    @classmethod
    def load(cls, path: str = "") -> DictFasta:
        """Parse a FASTA/FASTQ file into a dictionary of sequences.

        Args:
            path (str): Path to FASTA or FASTQ file

        Returns:
            DictFasta
        """
        if not path:
            path = EMPTY_FASTA

        if is_bgzipped(path):
            with gzip.open(path, "rt") as fh:
                fasta = SeqIO.to_dict(SeqIO.parse(fh, "fasta"))
        else:
            fasta = SeqIO.to_dict(SeqIO.parse(path, "fasta"))

        return cls(fasta)

    def __init__(self, fasta: Dict[str, str]):
        self.fasta = fasta

    def __contains__(self, reference: str) -> bool:
        return reference in self.fasta

    def __getitem__(self, reference: str) -> str:
        return self.fasta[reference]

    def fetch(self, reference: str, start: int, end: int) -> str:
        return str(self.fasta[reference][start:end])


class PyfaidxFasta(_FastaABC):
    @classmethod
    def load(cls, path: str = "") -> PyfaidxFasta:
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
            read_ahead=100000,
        )

        return cls(fasta)

    def __init__(self, fasta: Fasta):
        self.fasta = fasta

    def __contains__(self, reference: str) -> bool:
        return reference in self.fasta

    def __getitem__(self, reference: str) -> str:
        return self.fasta[reference]

    def fetch(self, reference: str, start: int, end: int) -> str:
        return str(self.fasta[reference][start:end])


@lru_cache
def get_sequence(
    fasta: _FastaABC,
    ref: str,
    start: int,
    end: int,
    window: Optional[int],
    floor: Optional[int],
    ceiling: Optional[int],
) -> str:
    """Return the reference sequence at the given position plus additional bases 5' and 3' such
    that the total returned sequence length adds up to `window`.

    Args:
        fasta (_FastaABC): FASTA containing all the reference sequences
        ref (str): Reference sequence ID
        start (int): First base of position of interest
        end (int): Last base of position of interest
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

    # Get enough of the reference sequence that we can return a sequence of `length`
    sequence = fasta.fetch(ref, ref_start, ref_end)

    return sequence


@lru_cache
def mutate_sequence(
    fasta: _FastaABC,
    ref: str,
    start: int,
    end: int,
    window: Optional[int],
    floor: Optional[int],
    ceiling: Optional[int],
    altseq: str,
    insertion: bool = False,
) -> str:
    """Return a sequence that represents the alternate allele in the context of the reference
    sequence within a window of length `window`.

    Args:
        fasta (_FastaABC): FASTA containing all the reference sequences
        ref (str): Reference sequence ID
        start (int): First base of the mutation site. For insertions, this is the base 5' of the
            insertion site
        end (int): Last base of the mutation site. For insertions, this is the base 3' of the
            insertion site
        window (Optional[int]): Total length of the returned sequence
        floor (Optional[int]): Absolute smallest start position (must be <= `start`)
        ceiling (Optional[int]): Absolute largest end position (must be >= `end`)
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

    idx_left = start - ref_start + ins_pad - 1
    idx_right = ref_end - end + ins_pad

    # If the floor is more than ref_start, make the new ref_start equal to floor and pad ref_end by
    # the difference
    floor = floor - 1 if (floor and floor > 0) else ref_start
    diff_left = floor - ref_start

    if ceiling and ceiling > 0:
        diff_right = ceiling - ref_end
        if diff_left > 0 and diff_left <= diff_right:
            ref_start -= diff_left
            ref_end -= diff_left
            idx_left -= diff_left
            idx_right += diff_left
        elif diff_right < 0 and diff_right <= diff_left:
            ref_start += diff_right
            ref_end += diff_right
            idx_left -= diff_right
            idx_right += diff_right
        else:
            raise ValueError("window constraints are smaller than the window length")
    elif diff_left > 0:
        ref_start += diff_left
        ref_end += diff_left
        idx_left -= diff_left
        idx_right += diff_left

    refseq = fasta.fetch(ref, ref_start, ref_end)

    # Piece together the altseq and 5'/3' flanking sequences
    left = refseq[:idx_left]
    right = refseq[-idx_right:] if idx_right > 0 else ""
    sequence = left + altseq + right

    return sequence
