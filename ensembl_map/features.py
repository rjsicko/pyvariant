from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from typing import Callable

import pyensembl

from .constants import CDS, CONTIG, EXON, GENE, PROTEIN, TRANSCRIPT
from .core import instance as CM


class Feature(ABC):
    """Every feature class requires these attributes."""

    @classmethod
    @abstractmethod
    def load(cls, transcript: pyensembl.Transcript, start: int, end: int) -> Feature:
        raise NotImplementedError


@dataclass(frozen=True, order=True)
class Cds(Feature):
    """CDS coordinate object.

    Attributes:
        biotype (str): biotype of the transcript
        contig (str): name of the contig the feature is mapped to
        end (int): end position, relative to the CDS
        sequence (str): CDS sequence from `start` to `end`, inclusive
        start (int): start position, relative to the CDS
        strand (str): orientation on the contig ("+" or "-")
        transcript (`pyensembl.Transcript`)
        transcript_id (str): Ensembl transcript ID
        transcript_name (str): Ensembl transcript name
    """

    transcript_id: str
    transcript_name: str
    start: int
    end: int
    strand: str
    biotype: str
    contig: str
    is_canonical: bool
    sequence: str

    @classmethod
    def load(cls, transcript: pyensembl.Transcript, start: int, end: int) -> Cds:
        if start > end:
            start, end = end, start

        sequence = CM().cdna[transcript.transcript_id]
        if sequence:
            sequence = sequence[start - 1 : end]

        return cls(
            biotype=transcript.transcript_biotype,
            contig=transcript.contig_id,
            end=end,
            is_canonical=CM().is_canonical_transcript(transcript.transcript_id),
            sequence=sequence,
            start=start,
            strand=transcript.strand,
            transcript_id=transcript.transcript_id,
            transcript_name=transcript.transcript_name,
        )


@dataclass(frozen=False, order=True)
class Contig(Feature):
    """Contig coordinate object.

    Attributes:
        contig (str): name of the contig the feature is mapped to
        end (int): end position, relative to the contig
        start (int): start position, relative to the contig
        strand (str): orientation on the contig ("+" or "-")
    """

    contig: str
    start: int
    end: int
    strand: str
    _sequence: str = field(default="", repr=False)

    def __hash__(self):
        return hash((self.contig, self.start, self.end, self.strand))

    @property
    def sequence(self) -> str:
        if not self._sequence:
            self._sequence = CM().contig_sequence(self.contig, self.start, self.end, self.strand)
        return self._sequence

    @classmethod
    def load(cls, transcript: pyensembl.Transcript, start: int, end: int) -> Contig:
        if start > end:
            start, end = end, start

        return cls(
            contig=transcript.contig_id,
            end=end,
            start=start,
            strand="+",  # describe relative to the positive strand to avoid confusion
        )


@dataclass(frozen=True, order=True)
class Exon(Feature):
    """Exon coordinate object.

    Attributes:
        biotype (str): biotype of the transcript the exon belongs to
        contig (str): name of the contig the exon is mapped to
        end (int): end position, relative to the contig
        exon (`pyensembl.Exon`)
        exon_id (str): Ensembl exon ID
        index (int): position of the exon relative to all exons in the transcript
        start (int): start position, relative to the contig
        strand (str): orientation on the contig ("+" or "-")
        transcript (`pyensembl.Transcript`)
        transcript_id (str): Ensembl transcript ID
        transcript_name (str): Ensembl transcript name
    """

    exon_id: str
    start: int
    end: int
    strand: str
    biotype: str
    contig: str
    transcript_id: str
    transcript_name: str
    index: int

    @classmethod
    def load(cls, transcript: pyensembl.Transcript, start: int, end: int) -> Exon:
        if start > end:
            start, end = end, start

        exon, index = cls.index_exon(transcript, start, end)

        return cls(
            biotype=transcript.transcript_biotype,
            contig=transcript.contig_id,
            end=end,
            exon_id=exon.exon_id,
            index=index,
            start=start,
            strand=transcript.strand,
            transcript_id=transcript.transcript_id,
            transcript_name=transcript.transcript_name,
        )

    @staticmethod
    def index_exon(transcript: pyensembl.Transcript, start: int, end: int) -> pyensembl.exon.Exon:
        exons = [
            (i, n) for n, i in enumerate(transcript.exons, 1) if start == i.start and end == i.end
        ]
        if len(exons) < 1:
            raise ValueError(
                f"Exon ({start}, {end}) not found in transcript {transcript.transcript_id}"
            )
        elif len(exons) > 1:
            raise ValueError(
                f"Multiple exons ({start}, {end}) found in transcript {transcript.transcript_id}"
            )
        else:
            return exons[0]


@dataclass(frozen=False, order=True)
class Gene(Feature):
    """Gene coordinate object.

    Attributes:
        biotype (str): biotype of the gene
        contig (str): name of the contig the feature is mapped to
        end (int): end position, relative to the contig
        gene (`pyensembl.Gene`)
        gene_id (str): Ensembl gene ID
        gene_name (str): Ensembl gene name
        start (int): start position, relative to the contig
        strand (str): orientation on the contig ("+" or "-")
    """

    gene_id: str
    gene_name: str
    start: int
    end: int
    strand: str
    biotype: str
    contig: str
    _sequence: str = field(default="", repr=False)

    def __hash__(self):
        return hash(
            (
                self.gene_id,
                self.gene_name,
                self.start,
                self.end,
                self.strand,
                self.biotype,
                self.contig,
            )
        )

    @property
    def sequence(self) -> str:
        if not self._sequence:
            self._sequence = CM().contig_sequence(self.contig, self.start, self.end, self.strand)
        return self._sequence

    @classmethod
    def load(cls, transcript: pyensembl.Transcript, start: int, end: int) -> Gene:
        if start > end:
            start, end = end, start

        gene = transcript.gene

        return cls(
            biotype=gene.biotype,
            contig=gene.contig,
            end=end,
            gene_id=gene.gene_id,
            gene_name=gene.gene_name,
            start=start,
            strand=gene.strand,
        )


@dataclass(frozen=True, order=True)
class Protein(Feature):
    """Protein coordinate object.

    Attributes:
        biotype (str): biotype of the transcript that encodes the protein (should be 'protein_coding')
        contig (str): name of the contig the protein is mapped to
        end (int): end position, relative to the protein
        protein_id (str): Ensembl protein ID
        start (int): start position, relative to the protein
        strand (str): orientation on the contig ("+" or "-")
        sequence (str): protein sequence from `start` to `end`, inclusive
        transcript (`pyensembl.Transcript`)
    """

    protein_id: str
    start: int
    end: int
    strand: None
    biotype: str
    contig: str
    sequence: str

    @classmethod
    def load(cls, transcript: pyensembl.Transcript, start: int, end: int) -> Protein:
        if transcript.transcript_biotype != "protein_coding":
            raise ValueError(
                f"{transcript.transcript_id} (biotype={transcript.transcript_biotype}) is not protein coding."
            )

        if start > end:
            start, end = end, start

        sequence = transcript.protein_sequence
        if sequence:
            sequence = sequence[start - 1 : end]

        return cls(
            biotype=transcript.transcript_biotype,
            contig=transcript.contig_id,
            end=end,
            protein_id=transcript.protein_id,
            sequence=sequence,
            start=start,
            strand=None,
        )


class Transcript(Cds):
    """Transcript coordinate object.

    Attributes:
        biotype (str): biotype of the transcript
        contig (str): name of the contig the feature is mapped to
        end (int): end position, relative to the transcript
        sequence (str): transcript sequence from `start` to `end`, inclusive
        start (int): start position, relative to the transcript
        strand (str): orientation on the contig ("+" or "-")
        transcript (`pyensembl.Transcript`)
        transcript_id (str): Ensembl transcript ID
        transcript_name (str): Ensembl transcript name
    """

    @classmethod
    def load(cls, transcript: pyensembl.Transcript, start: int, end: int) -> Transcript:
        if start > end:
            start, end = end, start

        transcript_id = transcript.transcript_id

        sequence = transcript.sequence
        if sequence:
            sequence = sequence[start - 1 : end]

        return cls(
            biotype=transcript.transcript_biotype,
            contig=transcript.contig_id,
            end=end,
            is_canonical=CM().is_best(transcript_id),
            sequence=sequence,
            start=start,
            strand=transcript.strand,
            transcript_id=transcript_id,
            transcript_name=transcript.transcript_name,
        )


def get_load_function(to_type: str) -> Callable:
    if to_type == CDS:
        return Cds.load
    elif to_type == CONTIG:
        return Contig.load
    elif to_type == EXON:
        return Exon.load
    elif to_type == GENE:
        return Gene.load
    elif to_type == PROTEIN:
        return Protein.load
    elif to_type == TRANSCRIPT:
        return Transcript.load
    else:
        raise TypeError(f"'{to_type}' is not a valid feature type")
