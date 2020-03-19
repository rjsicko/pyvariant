import logging

from .cache import Cache
from .util import is_ensembl_id


class CDS:
    """CDS coordinate object.

    Attributes:
        contig (str): name of the contig the feature is mapped to
        start (int): start position, relative to the CDS
        end (int): end position, relative to the CDS
        strand (str): orientation on the contig ("+" or "-")
        transcript_id (str): Ensembl transcript ID
        transcript_name (str): Ensembl transcript name
        seq (str): CDS sequence from `start` to `end`, inclusive
    """

    def __init__(self, contig, start, end, strand, transcript_id, transcript_name, seq):
        self.contig = contig
        self.start = start
        self.end = end
        self.strand = strand
        self.transcript_id = transcript_id
        self.transcript_name = transcript_name
        self.seq = seq

    @classmethod
    def format(cls, tobj, start, end):
        return cls(
            tobj.contig,
            start,
            end,
            tobj.strand,
            tobj.transcript_id,
            tobj.transcript_name,
            tobj.coding_sequence[start - 1 : end],
        )


class Exon:
    """Exon coordinate object.

    Attributes:
        contig (str): name of the contig the feature is mapped to
        start (int): start position, relative to the CDS
        end (int): end position, relative to the CDS
        strand (str): orientation on the contig ("+" or "-")
        exon_id (str): Ensembl exon ID
        transcript_id (str): Ensembl transcript ID
        transcript_name (str): Ensembl transcript name
        index (int): position of the exon relative to all exons in the transcript
        seq (str): CDS sequence from `start` to `end`, inclusive
    """

    def __init__(
        self,
        contig,
        start,
        end,
        strand,
        exon_id,
        transcript_id,
        transcript_name,
        index,
        seq,
    ):
        self.contig = contig
        self.start = start
        self.end = end
        self.strand = strand
        self.exon_id = exon_id
        self.transcript_id = transcript_id
        self.index = index
        self.seq = seq

    @classmethod
    def format(cls, tobj, start, end, exon_id, index):
        return cls(
            tobj.contig,
            start,
            end,
            tobj.strand,
            exon_id,
            tobj.transcript_id,
            tobj.transcript_name,
            index,
            tobj.sequence[start - 1 : end],
        )


class Gene:
    """Gene coordinate object.

    Attributes:
        contig (str): name of the contig the feature is mapped to
        start (int): start position, relative to the contig
        end (int): end position, relative to the contig
        strand (str): orientation on the contig ("+" or "-")
        gene_id (str): Ensembl gene ID
        gene_name (str): Ensembl gene name
    """

    def __init__(self, contig, start, end, strand, gene_id, gene_name):
        self.contig = contig
        self.start = start
        self.end = end
        self.strand = strand
        self.gene_id = gene_id
        self.gene_name = gene_name

    @classmethod
    def format(cls, tobj, start, end):
        return cls(tobj.gene_id, tobj.gene_name, start, end, tobj.contig, tobj.strand)


class Transcript:
    """Transcript coordinate object.

    Attributes:
        contig (str): name of the contig the feature is mapped to
        start (int): start position, relative to the transcript
        end (int): end position, relative to the transcript
        strand (str): orientation on the contig ("+" or "-")
        transcript_id (str): Ensembl transcript ID
        transcript_name (str): Ensembl transcript name
        seq (str): transcript sequence from `start` to `end`, inclusive
    """

    def __init__(self, contig, start, end, strand, transcript_id, transcript_name, seq):
        self.contig = contig
        self.start = start
        self.end = end
        self.strand = strand
        self.transcript_id = transcript_id
        self.transcript_name = transcript_name
        self.seq = seq

    @classmethod
    def format(cls, tobj, start, end):
        return cls(
            tobj.contig,
            start,
            end,
            tobj.transcript_id,
            tobj.transcript_name,
            tobj.strand,
            tobj.sequence[start - 1 : end],
        )

