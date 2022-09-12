"""
3' rule: for all descriptions the most 3' position possible of the reference sequence is arbitrarily assigned to have been changed:
    * the 3' rule also applies for changes in single residue stretches and tandem repeats (nucleotide or amino acid)
    * the 3' rule applies to ALL descriptions (genome, gene, transcript and protein) of a given variant
    * exception: deletion/duplication around exon/exon junctions using c., r. or n. reference sequences

CACACTGTGGAT to 
CACACTG--GAT is c.8_9del

CACACTGTG--GAT to 
CACACTGTGTGGAT is c.8_9dup (not c.9_10insTG)
"""
from functools import lru_cache
from typing import Optional, Tuple

import pyensembl
from logzero import logger

from .features import Cdna, Protein, Transcript
from .map import (
    cdna_to_transcript,
    protein_to_transcript,
    transcript_to_cdna,
    transcript_to_protein,
    transcript_to_transcript,
)


@lru_cache()
def normalize_cdna(feature: str, start: int, end: Optional[int] = None) -> Cdna:
    """3' align a CDNA variant according to HGVS rules.

    Args:
        feature (str): transcript name or Ensembl ID
        start (int): first position on the transcript
        start (int): second position on the transcript

    Returns:
        CDNA: object with the re-aligned coordinates
    """
    transcript = cdna_to_transcript(feature, start, end)[0]
    start_shift, end_shift = _normalize(transcript)
    normalized = transcript_to_cdna(transcript.transcript_id, start_shift, end_shift)[0]

    in_coordinate = (feature, start, end)
    out_coordinate = (normalized.transcript_id, normalized.start, normalized.end)
    if in_coordinate != out_coordinate:
        logger.debug(f"Normalized {in_coordinate} to {out_coordinate}")
    else:
        logger.debug(f"No change in normalized position: {in_coordinate}")

    return normalized


@lru_cache()
def normalize_protein(feature: str, start: int, end: Optional[int] = None) -> Protein:
    """3' align a protein variant according to HGVS rules.

    Args:
        feature (str): protein name or Ensembl ID
        start (int): first position on the protein
        start (int): second position on the protein

    Returns:
        Protein: object with the re-aligned coordinates
    """
    transcript = protein_to_transcript(feature, start, end)[0]
    start_shift, end_shift = _normalize(transcript)
    normalized = transcript_to_protein(transcript.transcript_id, start_shift, end_shift)[0]

    in_coordinate = (feature, start, end)
    out_coordinate = (normalized.protein_id, normalized.start, normalized.end)
    if in_coordinate != out_coordinate:
        logger.debug(f"Normalized {in_coordinate} to {out_coordinate}")
    else:
        logger.debug(f"No change in normalized position: {in_coordinate}")

    return normalized


@lru_cache()
def normalize_transcript(feature: str, start: int, end: Optional[int] = None) -> Transcript:
    """3' align a transcript variant according to HGVS rules.

    Args:
        feature (str): transcript name or Ensembl ID
        start (int): first position on the transcript
        start (int): second position on the transcript

    Returns:
        Transcript: object with the re-aligned coordinates
    """
    transcript = _map_transcript(feature, start, end)
    start_shift, end_shift = _normalize(transcript)
    normalized = _map_transcript(feature, start_shift, end_shift)

    in_coordinate = (feature, start, end)
    out_coordinate = (normalized.transcript_id, normalized.start, normalized.end)
    if in_coordinate != out_coordinate:
        logger.debug(f"Normalized {in_coordinate} to {out_coordinate}")
    else:
        logger.debug(f"No change in normalized position: {in_coordinate}")

    return normalized


def _normalize(transcript: pyensembl.Transcript) -> Tuple[int, int]:
    """3' align a transcript variant according to HGVS rules."""
    # collapse the reference sequence to the shortest, non-overlapping substring (e.g. "GG" -> "G")
    collapsed = _collapse_repeats(transcript.sequence)
    end_new = transcript.start + len(collapsed) - 1
    if collapsed != transcript.sequence:
        logger.debug(
            f"Collapsed transcript sequence: {transcript.sequence} {transcript.start, transcript.end} -> {collapsed} {transcript.start, end_new}"
        )

    # shift the position
    start_shift, end_shift = _three_prime_rule(transcript.transcript_id, transcript.start, end_new)
    if (transcript.start, end_new) != (start_shift, end_shift):
        logger.debug(
            f"Shifted postion 3' by {start_shift - transcript.start}bp: {transcript.start, end_new} -> {start_shift, end_shift}"
        )

    return start_shift, end_shift


def _map_transcript(feature: str, start: int, end: Optional[int]) -> Transcript:
    """Convert a transcript coordinate to a ``Transcript`` object."""
    result = transcript_to_transcript(feature, start, end)
    if len(result) != 1:
        raise ValueError(f"Unable to uniquely map transcript {feature, start, end}: {result}")
    else:
        return result[0]


def _collapse_repeats(sequence: str) -> str:
    """Collapse a series of tandem repeats to the shortest substring."""
    # "GGGGG" -> "G"
    # "GTAGTA" -> "GTA"
    # "GTAGTAG" -> "GTAGTAG"

    for n in range(1, len(sequence)):
        # count the number of non-overlapping substrings of lenght `n`
        x = set([sequence[i : i + n] for i in range(0, len(sequence), n)])
        if len(x) == 1:
            return x.pop()
    else:
        return sequence


def _three_prime_rule(feature: str, start: int, end: int) -> Tuple[int, int]:
    """Shift the position as 3' as possible on the transcript."""
    window = end - start + 1
    while True:
        transcript = _map_transcript(feature, start, end + window)
        ref = transcript.sequence[:window]
        ref_lookahead = transcript.sequence[window:]
        if ref == ref_lookahead:
            start += window
            end += window
        else:
            return start, end
