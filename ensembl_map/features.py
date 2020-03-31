from types import SimpleNamespace


class CDS(SimpleNamespace):
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

    @classmethod
    def load(cls, transcript, start, end):
        if start > end:
            start, end = end, start

        sequence = getattr(transcript, "coding_sequence", None)
        if sequence:
            sequence = sequence[start - 1 : end]

        return cls(
            biotype=getattr(transcript, "biotype", None),
            contig=getattr(transcript, "contig", None),
            end=end,
            sequence=sequence,
            start=start,
            strand=getattr(transcript, "strand", None),
            transcript_id=getattr(transcript, "transcript_id", None),
            transcript_name=getattr(transcript, "transcript_name", None),
        )

    def to_tuple(self):
        return self.transcript_id, self.start, self.end


class Exon(SimpleNamespace):
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

    @classmethod
    def load(cls, transcript, start, end, exon_id):
        if start > end:
            start, end = end, start

        exon, index = cls.index_exon(transcript, exon_id)

        return cls(
            biotype=getattr(transcript, "biotype", None),
            contig=getattr(transcript, "contig", None),
            end=end,
            exon_id=getattr(exon, "exon_id", None),
            index=index,
            start=start,
            strand=getattr(transcript, "strand", None),
            transcript_id=getattr(transcript, "transcript_id", None),
            transcript_name=getattr(transcript, "transcript_name", None),
        )

    @staticmethod
    def index_exon(transcript, exon_id):
        exons = [(i, n) for n, i in enumerate(transcript.exons, 1) if i.exon_id == exon_id]
        if len(exons) < 1:
            raise ValueError(f"Exon {exon_id} not found in transcript {transcript.transcript_id}")
        elif len(exons) > 1:
            raise ValueError(
                f"Multiple exons {exon_id} found in transcript {transcript.transcript_id}"
            )
        else:
            return exons[0]

    def to_tuple(self):
        return self.exon_id, self.start, self.end


class Gene(SimpleNamespace):
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

    @classmethod
    def load(cls, transcript, start, end):
        if start > end:
            start, end = end, start

        gene = transcript.gene

        return cls(
            biotype=getattr(gene, "biotype", None),
            contig=getattr(gene, "contig", None),
            end=end,
            gene_id=getattr(gene, "gene_id", None),
            gene_name=getattr(gene, "gene_name", None),
            start=start,
            strand=getattr(gene, "strand", None),
        )

    def to_tuple(self):
        return self.gene_id, self.start, self.end


class Protein(SimpleNamespace):
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

    @classmethod
    def load(cls, transcript, start, end):
        if start > end:
            start, end = end, start

        sequence = getattr(transcript, "protein_sequence", None)
        if sequence:
            sequence = sequence[start - 1 : end]

        return cls(
            biotype=getattr(transcript, "biotype", None),
            contig=getattr(transcript, "contig", None),
            end=end,
            protein_id=getattr(transcript, "protein_id", None),
            sequence=sequence,
            start=start,
            strand=getattr(transcript, "strand", None),
        )

    def to_tuple(self):
        return self.protein_id, self.start, self.end


class Transcript(CDS):
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
    def load(cls, transcript, start, end):
        if start > end:
            start, end = end, start

        sequence = getattr(transcript, "sequence", None)
        if sequence:
            sequence = sequence[start - 1 : end]

        return cls(
            biotype=getattr(transcript, "biotype", None),
            contig=getattr(transcript, "contig", None),
            end=end,
            sequence=sequence,
            start=start,
            strand=getattr(transcript, "strand", None),
            transcript_id=getattr(transcript, "transcript_id", None),
            transcript_name=getattr(transcript, "transcript_name", None),
        )


def get_load_function(to_type):
    if to_type == "cds":
        return CDS.load
    elif to_type == "exon":
        return Exon.load
    elif to_type == "gene":
        return Gene.load
    elif to_type == "protein":
        return Protein.load
    elif to_type == "transcript":
        return Transcript.load
    else:
        raise TypeError(f"Could not get parse function for {to_type}")
