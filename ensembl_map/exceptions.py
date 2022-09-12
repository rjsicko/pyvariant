from typing import Tuple

import pyensembl

from .constants import GENE, TRANSCRIPT


class OutOfRangeErrorBase(ValueError):
    """Base error for 'OutOfRange' errors."""

    def __init__(self):
        self.feature_id: str = ""
        self.feature_type: str = ""
        self.position: int = 0
        self.range: Tuple[int, int] = (0, 0)

    def __str__(self):
        return f"{self.position} is outside {self.feature_type} {self.feature_id} {self.range}"


class CdnaOutOfRange(OutOfRangeErrorBase):
    def __init__(self, transcript: pyensembl.Transcript, position: int):
        self.feature_id = transcript.transcript_id
        self.position = position
        self.feature_type = "cDNA"
        self.range = (1, len(transcript.coding_sequence))


class ExonOutOfRange(OutOfRangeErrorBase):
    def __init__(self, transcript: pyensembl.Transcript, position: int):
        self.feature_id = transcript.transcript_id
        self.position = position
        self.feature_type = "exons"
        self.range = transcript.exon_intervals


class GeneOutOfRange(OutOfRangeErrorBase):
    def __init__(self, gene: pyensembl.Gene, position: int):
        self.feature_id = gene.gene_id
        self.position = position
        self.feature_type = GENE
        self.range = (gene.start, gene.end)


class TranscriptOutOfRange(OutOfRangeErrorBase):
    def __init__(self, transcript: pyensembl.Transcript, position: int):
        self.feature_id = transcript.transcript_id
        self.position = position
        self.feature_type = TRANSCRIPT
        self.range = (1, len(transcript.sequence))
