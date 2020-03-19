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

    def __repr__(self):
        if len(self.seq) > 3:
            seq_str = self.seq[:3] + "..."
        else:
            seq_str = self.seq
        return (
            f"CDS("
            f"contig={self.contig}, "
            f"start={self.start}, "
            f"end={self.end}, "
            f"strand={self.strand}, "
            f"transcript_id={self.transcript_id}, "
            f"transcript_name={self.transcript_name}, "
            f"seq={seq_str})"
        )

    @classmethod
    def load(cls, tobj, start, end):
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
        self.transcript_name = transcript_name
        self.index = index
        self.seq = seq

    @classmethod
    def load(cls, tobj, start, end, exon_id, index):
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
    def load(cls, tobj, start, end):
        return cls(tobj.contig, start, end, tobj.strand, tobj.gene_id, tobj.gene_name)


class Protein:
    """Protein coordinate object.

    Attributes:
        contig (str): name of the contig the feature is mapped to
        start (int): start position, relative to the transcript
        end (int): end position, relative to the transcript
        strand (str): orientation on the contig ("+" or "-")
        protein_id_id (str): Ensembl protein ID
        seq (str): transcript sequence from `start` to `end`, inclusive
    """

    def __init__(self, contig, start, end, strand, protein_id, seq):
        self.contig = contig
        self.start = start
        self.end = end
        self.strand = strand
        self.protein_id = protein_id
        self.seq = seq

    @classmethod
    def load(cls, tobj, start, end):
        return cls(
            tobj.contig,
            start,
            end,
            tobj.strand,
            tobj.protein_id,
            tobj.protein_sequence[start - 1 : end],
        )


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
    def load(cls, tobj, start, end):
        return cls(
            tobj.contig,
            start,
            end,
            tobj.strand,
            tobj.transcript_id,
            tobj.transcript_name,
            tobj.sequence[start - 1 : end],
        )

