class FeatureBase:
    _seqlim = 3  # only print the first N bases of a sequence

    def __init__(self, pyensembl_obj, start, end):
        self._pyensembl_obj = pyensembl_obj
        self.start = start
        self.end = end

    def __repr__(self):
        return f"{self.__class__.__name__}{self.to_tuple()}"

    @property
    def biotype(self):
        try:
            return self._pyensembl_obj.biotype
        except AttributeError:
            return None

    @property
    def contig(self):
        try:
            return self._pyensembl_obj.contig
        except AttributeError:
            return None

    @property
    def strand(self):
        try:
            return self._pyensembl_obj.strand
        except AttributeError:
            return None

    def to_tuple(self):
        return None, None, None


class CDS(FeatureBase):
    """CDS coordinate object.

    Attributes:
        biotype (str): biotype of the transcript
        contig (str): name of the contig the feature is mapped to
        start (int): start position, relative to the CDS
        end (int): end position, relative to the CDS
        strand (str): orientation on the contig ("+" or "-")
        transcript_id (str): Ensembl transcript ID
        transcript_name (str): Ensembl transcript name
        sequence (str): CDS sequence from `start` to `end`, inclusive
    """

    @classmethod
    def load(cls, transcript, start, end):
        if start > end:
            start, end = end, start
        return cls(transcript, start, end)

    @property
    def transcript(self):
        return self._pyensembl_obj

    @property
    def transcript_id(self):
        try:
            return self._pyensembl_obj.transcript_id
        except AttributeError:
            return None

    @property
    def transcript_name(self):
        try:
            return self._pyensembl_obj.transcript_name
        except AttributeError:
            return None

    @property
    def sequence(self):
        try:
            return self._pyensembl_obj.coding_sequence[self.start - 1 : self.end]
        except AttributeError:
            return None

    def to_tuple(self):
        return self.transcript_id, self.start, self.end


class Exon(FeatureBase):
    """Exon coordinate object.

    Attributes:
        biotype (str): biotype of the gene
        contig (str): name of the contig the feature is mapped to
        start (int): start position, relative to the CDS
        end (int): end position, relative to the CDS
        strand (str): orientation on the contig ("+" or "-")
        exon_id (str): Ensembl exon ID
        index (int): position of the exon relative to all exons in the transcript
        seq (str): CDS sequence from `start` to `end`, inclusive
    """

    def __init__(self, pyensembl_obj, exon, start, end):
        self._pyensembl_obj = pyensembl_obj
        self._exon = exon
        self.start = start
        self.end = end

    @classmethod
    def load(cls, transcript, start, end, exon_id):
        if start > end:
            start, end = end, start
        _, eobj = cls.get_exon_from_transcript(transcript, exon_id)
        return cls(transcript, eobj, start, end)

    @classmethod
    def get_exon_from_transcript(cls, transcript, exon_id):
        exons = [(n, i) for n, i in enumerate(transcript.exons, 1) if i.exon_id == exon_id]
        if len(exons) < 1:
            raise ValueError(f"Exon {exon_id} not found in transcript {transcript.transcript_id}")
        elif len(exons) > 1:
            raise ValueError(
                f"Multiple exons {exon_id} found in transcript {transcript.transcript_id}"
            )
        else:
            return exons[0]

    @property
    def exon(self):
        return self._exon

    @property
    def exon_id(self):
        try:
            return self._exon.exon_id
        except AttributeError:
            return None

    @property
    def exon_position(self):
        return self.get_exon_from_transcript(self.transcript, self.exon_id)[0]

    @property
    def transcript(self):
        return self._pyensembl_obj

    @property
    def transcript_id(self):
        try:
            return self._pyensembl_obj.transcript_id
        except AttributeError:
            return None

    @property
    def transcript_name(self):
        try:
            return self._pyensembl_obj.transcript_name
        except AttributeError:
            return None

    def to_tuple(self):
        return self.exon_id, self.start, self.end


class Gene(FeatureBase):
    """Gene coordinate object.

    Attributes:
        biotype (str): biotype of the gene
        contig (str): name of the contig the feature is mapped to
        start (int): start position, relative to the contig
        end (int): end position, relative to the contig
        strand (str): orientation on the contig ("+" or "-")
        gene_id (str): Ensembl gene ID
        gene_name (str): Ensembl gene name
    """

    @classmethod
    def load(cls, transcript, start, end):
        if start > end:
            start, end = end, start
        return cls(transcript.gene, start, end)

    @property
    def gene(self):
        return self._pyensembl_obj

    @property
    def gene_id(self):
        try:
            return self._pyensembl_obj.gene_id
        except AttributeError:
            return None

    @property
    def gene_name(self):
        try:
            return self._pyensembl_obj.gene_name
        except AttributeError:
            return None

    def to_tuple(self):
        return self.gene_id, self.start, self.end


class Protein(FeatureBase):
    """Protein coordinate object.

    Attributes:
        biotype (str): biotype of the gene 
        contig (str): name of the contig the feature is mapped to
        start (int): start position, relative to the transcript
        end (int): end position, relative to the transcript
        strand (str): orientation on the contig ("+" or "-")
        protein_id (str): Ensembl protein ID
        sequence (str): transcript sequence from `start` to `end`, inclusive
    """

    @classmethod
    def load(cls, transcript, start, end):
        if start > end:
            start, end = end, start
        return cls(transcript, start, end)

    @property
    def protein_id(self):
        try:
            return self._pyensembl_obj.protein_id
        except AttributeError:
            return None

    @property
    def sequence(self):
        return self._pyensembl_obj.protein_sequence[self.start - 1 : self.end]

    @property
    def transcript(self):
        return self._pyensembl_obj

    def to_tuple(self):
        return self.protein_id, self.start, self.end


class Transcript(CDS):
    """Transcript coordinate object.

    Attributes:
        biotype (str): biotype of the gene
        contig (str): name of the contig the feature is mapped to
        start (int): start position, relative to the transcript
        end (int): end position, relative to the transcript
        strand (str): orientation on the contig ("+" or "-")
        transcript_id (str): Ensembl transcript ID
        transcript_name (str): Ensembl transcript name
        sequence (str): transcript sequence from `start` to `end`, inclusive
    """

    @property
    def sequence(self):
        try:
            return self._pyensembl_obj.sequence[self.start - 1 : self.end]
        except AttributeError:
            return None


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
