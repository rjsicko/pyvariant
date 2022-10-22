from __future__ import annotations

from dataclasses import dataclass, fields
from itertools import product
from math import floor
from typing import Any, Callable, Dict, List, Optional, Tuple, Union

import pandas as pd
from logzero import logger
from pyfaidx import Fasta

from .constants import (
    CONTIG_ID,
    DEFAULT_RELEASE,
    DEFAULT_SPECIES,
    EXON_ID,
    GENE_ID,
    GENE_NAME,
    PROTEIN_ID,
    TRANSCRIPT_ID,
    TRANSCRIPT_NAME,
)
from .ensembl_cache import EnsemblCache
from .utils import reverse_complement, strip_version


# -------------------------------------------------------------------------------------------------
# position classes
# -------------------------------------------------------------------------------------------------
@dataclass(eq=True, frozen=True, order=True)
class _Position:
    _data: Core
    contig_id: str
    start: int
    end: int
    strand: str

    @classmethod
    def copy_from(cls, position: _Position, **kwargs):
        return cls(**{**position.asdict(), **kwargs})

    def __getitem__(self, item):
        return getattr(self, item)

    @property
    def length(self) -> int:
        return self.end - self.start + 1

    def asdict(self) -> Dict[str, Any]:
        return {f.name: self[f.name] for f in fields(self)}

    def to_str(self) -> str:
        raise NotImplementedError()

    def on_negative_strand(self) -> bool:
        return self.strand == "-"

    def on_positive_strand(self) -> bool:
        return self.strand == "+"

    def sequence(self) -> str:
        raise NotImplementedError()

    def to_cdna(self) -> List[CdnaPosition]:
        raise NotImplementedError()

    def to_dna(self) -> List[DnaPosition]:
        raise NotImplementedError()

    def to_exon(self) -> List[ExonPosition]:
        raise NotImplementedError()

    def to_protein(self) -> List[ProteinPosition]:
        raise NotImplementedError()

    def to_rna(self) -> List[RnaPosition]:
        raise NotImplementedError()


@dataclass(eq=True, frozen=True, order=True)
class CdnaPosition(_Position):
    gene_id: str
    gene_name: str
    transcript_id: str
    transcript_name: str
    protein_id: str

    def to_str(self) -> str:
        if self.start == self.end:
            return f"{self.transcript_id}:c.{self.start}"
        else:
            return f"{self.transcript_id}:c.{self.start}-{self.end}"

    def sequence(self) -> str:
        return self._data.cdna_sequence(self.transcript_id, self.start, self.end, self.strand)

    def to_cdna(self) -> List[CdnaPosition]:
        return self._data._cdna_to_cdna([self.transcript_id], self.start, self.end, [self.strand])

    def to_dna(self) -> List[DnaPosition]:
        return self._data._cdna_to_dna([self.transcript_id], self.start, self.end, [self.strand])

    def to_exon(self) -> List[ExonPosition]:
        return self._data._cdna_to_exon([self.transcript_id], self.start, self.end, [self.strand])

    def to_rna(self) -> List[RnaPosition]:
        return self._data._cdna_to_rna([self.transcript_id], self.start, self.end, [self.strand])

    def to_protein(self) -> List[ProteinPosition]:
        return self._data._cdna_to_protein(
            [self.transcript_id], self.start, self.end, [self.strand]
        )


@dataclass(eq=True, frozen=True, order=True)
class DnaPosition(_Position):
    def to_str(self) -> str:
        if self.start == self.end:
            return f"{self.contig_id}:g.{self.start}"
        else:
            return f"{self.contig_id}:g.{self.start}-{self.end}"

    @classmethod
    def load(cls, row: pd.Series = None, **kwargs):
        return cls(
            *[
                kwargs[f.name] if f.name in kwargs else row[f.name] if row is not None else None
                for f in fields(cls)
            ]
        )

    def sequence(self) -> str:
        return self._data.dna_sequence(self.contig_id, self.start, self.end, self.strand)

    def to_cdna(self) -> List[CdnaPosition]:
        return self._data._dna_to_cdna([self.contig_id], self.start, self.end, [self.strand])

    def to_dna(self) -> List[DnaPosition]:
        return self._data._dna_to_dna([self.contig_id], self.start, self.end, [self.strand])

    def to_exon(self) -> List[ExonPosition]:
        return self._data._dna_to_exon([self.contig_id], self.start, self.end, [self.strand])

    def to_protein(self) -> List[ProteinPosition]:
        return self._data._dna_to_protein([self.contig_id], self.start, self.end, [self.strand])

    def to_rna(self) -> List[RnaPosition]:
        return self._data._dna_to_rna([self.contig_id], self.start, self.end, [self.strand])


@dataclass(eq=True, frozen=True, order=True)
class ExonPosition(_Position):
    gene_id: str
    gene_name: str
    transcript_id: str
    transcript_name: str
    exon_id: str

    def to_str(self) -> str:
        if self.start == self.end:
            return f"{self.exon_id}:e.{self.start}"
        else:
            return f"{self.exon_id}:e.{self.start}-{self.end}"

    def sequence(self) -> str:
        raise NotImplementedError()  # TODO

    def to_cdna(self) -> List[CdnaPosition]:
        return self._data._exon_to_cdna([self.transcript_id], self.start, self.end, [self.strand])

    def to_dna(self) -> List[DnaPosition]:
        return self._data._exon_to_dna([self.transcript_id], self.start, self.end, [self.strand])

    def to_exon(self) -> List[ExonPosition]:
        return self._data._exon_to_exon([self.transcript_id], self.start, self.end, [self.strand])

    def to_protein(self) -> List[ProteinPosition]:
        return self._data._exon_to_protein(
            [self.transcript_id], self.start, self.end, [self.strand]
        )

    def to_rna(self) -> List[RnaPosition]:
        return self._data._exon_to_rna([self.transcript_id], self.start, self.end, [self.strand])


@dataclass(eq=True, frozen=True, order=True)
class ProteinPosition(_Position):
    gene_id: str
    gene_name: str
    transcript_id: str
    transcript_name: str
    protein_id: str

    def to_str(self) -> str:
        if self.start == self.end:
            return f"{self.protein_id}:p.{self.start}"
        else:
            return f"{self.protein_id}:p.{self.start}-{self.end}"

    def sequence(self) -> str:
        return self._data.peptide_sequence(self.protein_id, self.start, self.end, self.strand)

    def to_cdna(self) -> List[CdnaPosition]:
        return self._data._protein_to_cdna(
            [self.transcript_id], self.start, self.end, [self.strand]
        )

    def to_dna(self) -> List[DnaPosition]:
        return self._data._protein_to_dna([self.transcript_id], self.start, self.end, [self.strand])

    def to_exon(self) -> List[ExonPosition]:
        return self._data._protein_to_exon(
            [self.transcript_id], self.start, self.end, [self.strand]
        )

    def to_protein(self) -> List[ProteinPosition]:
        return self._data._protein_to_protein(
            [self.transcript_id], self.start, self.end, [self.strand]
        )

    def to_rna(self) -> List[RnaPosition]:
        return self._data._protein_to_rna([self.transcript_id], self.start, self.end, [self.strand])


@dataclass(eq=True, frozen=True, order=True)
class RnaPosition(_Position):
    gene_id: str
    gene_name: str
    transcript_id: str
    transcript_name: str

    def to_str(self) -> str:
        if self.start == self.end:
            return f"{self.transcript_id}:r.{self.start}"
        else:
            return f"{self.transcript_id}:r.{self.start}-{self.end}"

    def sequence(self) -> str:
        return self._data.rna_sequence(self.transcript_id, self.start, self.end, self.strand)

    def to_cdna(self) -> List[CdnaPosition]:
        return self._data._rna_to_cdna([self.transcript_id], self.start, self.end, [self.strand])

    def to_dna(self) -> List[DnaPosition]:
        return self._data._rna_to_dna([self.transcript_id], self.start, self.end, [self.strand])

    def to_exon(self) -> List[ExonPosition]:
        return self._data._rna_to_exon([self.transcript_id], self.start, self.end, [self.strand])

    def to_protein(self) -> List[ProteinPosition]:
        return self._data._rna_to_protein([self.transcript_id], self.start, self.end, [self.strand])

    def to_rna(self) -> List[RnaPosition]:
        return self._data._rna_to_rna([self.transcript_id], self.start, self.end, [self.strand])


# -------------------------------------------------------------------------------------------------
# cache logic
# -------------------------------------------------------------------------------------------------
class CachedCore(type):
    _instances: Dict[Any, Core] = {}
    _current: Optional[Core] = None

    def __call__(cls, *args, **kwargs):
        # convert the given arguments to a hashable tuple
        key = tuple(list(args) + [kwargs[k] for k in sorted(kwargs)])
        # check if an instance created with the same arguments is already cached
        if key in cls._instances:
            instance = cls._instances[key]
        else:
            # create a new instance and cache it
            instance = super().__call__(*args, **kwargs)
            cls._instances[key] = instance

        # track which instance was called last
        cls._current = instance

        return cls._instances[key]


class Core(metaclass=CachedCore):
    def __init__(
        self,
        df: pd.DataFrame,
        cdna: Fasta,
        dna: Fasta,
        pep: Fasta,
        ncrna: Fasta,
        canonical_transcript: List[str] = [],
        contig_alias: Dict[str, List[str]] = {},
        exon_alias: Dict[str, List[str]] = {},
        gene_alias: Dict[str, List[str]] = {},
        protein_alias: Dict[str, List[str]] = {},
        transcript_alias: Dict[str, List[str]] = {},
    ):
        """Load annotations for the given release."""
        self.df = df
        self.cdna = cdna
        self.dna = dna
        self.pep = pep
        self.ncrna = ncrna
        self.canonical_transcript = canonical_transcript
        self.contig_alias = contig_alias
        self.exon_alias = exon_alias
        self.gene_alias = gene_alias
        self.protein_alias = protein_alias
        self.transcript_alias = transcript_alias

    # ---------------------------------------------------------------------------------------------
    # all_<feature_symbol>s
    # ---------------------------------------------------------------------------------------------
    def all_contig_ids(self) -> List[str]:
        return self._uniquify_series(self.df[CONTIG_ID])

    def all_exon_ids(self) -> List[str]:
        return self._uniquify_series(self.df[EXON_ID])

    def all_gene_ids(self) -> List[str]:
        return self._uniquify_series(self.df[GENE_ID])

    def all_gene_names(self) -> List[str]:
        return self._uniquify_series(self.df[GENE_NAME])

    def all_protein_ids(self) -> List[str]:
        return self._uniquify_series(self.df[PROTEIN_ID])

    def all_transcript_ids(self) -> List[str]:
        return self._uniquify_series(self.df[TRANSCRIPT_ID])

    def all_transcript_names(self) -> List[str]:
        return self._uniquify_series(self.df[TRANSCRIPT_NAME])

    def _uniquify_series(self, series: pd.Series) -> List:
        return sorted(series.dropna().unique().tolist())

    # ---------------------------------------------------------------------------------------------
    # get_<feature_symbol>
    # ---------------------------------------------------------------------------------------------
    def contig_ids(self, feature: str, feature_type: str = "") -> List[str]:
        """Given a feature symbol, return the corresponding contig ID(s)."""
        return self._get_feature_attr(CONTIG_ID, feature, feature_type)

    def exon_ids(self, feature: str, feature_type: str = "") -> List[str]:
        """Given a feature symbol, return the corresponding exon ID(s)."""
        return self._get_feature_attr(EXON_ID, feature, feature_type)

    def gene_ids(self, feature: str, feature_type: str = "") -> List[str]:
        """Given a feature symbol, return the corresponding gene ID(s)."""
        return self._get_feature_attr(GENE_ID, feature, feature_type)

    def gene_names(self, feature: str, feature_type: str = "") -> List[str]:
        """Given a feature symbol, return the corresponding gene names(s)."""
        return self._get_feature_attr(GENE_NAME, feature, feature_type)

    def protein_ids(self, feature: str, feature_type: str = "") -> List[str]:
        """Given a feature symbol, return the corresponding protein ID(s)."""
        return self._get_feature_attr(PROTEIN_ID, feature, feature_type)

    def transcript_ids(self, feature: str, feature_type: str = "") -> List[str]:
        """Given a feature symbol, return the corresponding transcript ID(s)."""
        return self._get_feature_attr(TRANSCRIPT_ID, feature, feature_type)

    def transcript_names(self, feature: str, feature_type: str = "") -> List[str]:
        """Given a feature symbol, return the corresponding transcript names(s)."""
        return self._get_feature_attr(TRANSCRIPT_NAME, feature, feature_type)

    def _get_feature_attr(self, key: str, feature: str, feature_type: str = "") -> List[str]:
        """Given a feature symbol, return the corresponding ID(s)."""
        parts = []

        for feature, feature_type in self.normalize_feature(feature, feature_type):
            if feature_type == CONTIG_ID:
                func = self._query_by_contig_id
            elif feature_type == EXON_ID:
                func = self._query_by_exon_id
            elif feature_type == GENE_ID:
                func = self._query_by_gene_id
            elif feature_type == GENE_NAME:
                func = self._query_by_gene_name
            elif feature_type == PROTEIN_ID:
                func = self._query_by_protein_id
            elif feature_type == TRANSCRIPT_ID:
                func = self._query_by_transcript_id
            elif feature_type == TRANSCRIPT_NAME:
                func = self._query_by_transcript_name
            else:
                raise ValueError(f"Unable to get {key} for {feature} ({feature_type})")

            parts.append(func(feature)[key])

        if parts:
            result = pd.concat(parts)
        else:
            result = pd.Series(dtype="object")

        return self._uniquify_series(result)

    def _query_by_contig_id(self, feature: Union[List[str], str]) -> pd.DataFrame:
        return self._query(feature, CONTIG_ID)

    def _query_by_exon_id(self, feature: Union[List[str], str]) -> pd.DataFrame:
        return self._query(feature, EXON_ID)

    def _query_by_gene_id(self, feature: Union[List[str], str]) -> pd.DataFrame:
        return self._query(feature, GENE_ID)

    def _query_by_gene_name(self, feature: Union[List[str], str]) -> pd.DataFrame:
        return self._query(feature, GENE_NAME)

    def _query_by_protein_id(self, feature: Union[List[str], str]) -> pd.DataFrame:
        return self._query(feature, PROTEIN_ID)

    def _query_by_transcript_id(self, feature: Union[List[str], str]) -> pd.DataFrame:
        return self._query(feature, TRANSCRIPT_ID)

    def _query_by_transcript_name(self, feature: Union[List[str], str]) -> pd.DataFrame:
        return self._query(feature, TRANSCRIPT_NAME)

    def _query(self, feature: Union[List[str], str], col: str) -> pd.DataFrame:
        """Generic function for querying the data cache."""
        feature = [feature] if isinstance(feature, str) else feature
        sudbf = self.df.loc[self.df[col].isin(feature)]

        return sudbf

    # ---------------------------------------------------------------------------------------------
    # normalize_feature
    # ---------------------------------------------------------------------------------------------
    def normalize_feature(self, feature: str, feature_type: str = "") -> List[Tuple[str, str]]:
        """Normalize a feature ID to the representation used by Ensembl."""
        normalized = []

        feature_type, result = self._normalize_feature(feature, feature_type=feature_type)
        if feature_type:
            normalized = self._uniquify_series(result[feature_type])

        return [(i, feature_type) for i in normalized]

    def _normalize_feature(self, feature: str, feature_type: str = "") -> Tuple[str, pd.DataFrame]:
        """Normalize a feature to the representation used by Ensembl."""
        result = pd.DataFrame()

        for key, func in [
            (CONTIG_ID, self._normalize_contig_id),
            (EXON_ID, self._normalize_exon_id),
            (GENE_ID, self._normalize_gene_id),
            (GENE_NAME, self._normalize_gene_name),
            (PROTEIN_ID, self._normalize_protein_id),
            (TRANSCRIPT_ID, self._normalize_transcript_id),
            (TRANSCRIPT_NAME, self._normalize_transcript_name),
        ]:
            if feature_type == key or not feature_type:
                result = func(feature)
                if not result.empty:
                    feature_type = key
                    break

        return feature_type, result

    def _normalize_contig_id(self, feature: str) -> pd.DataFrame:
        """Normalize a contig ID or it's alias to one or more matching Ensembl contig ID."""
        featurel = [feature] + self.get_contig_alias(feature)

        return self._query_by_contig_id(featurel)

    def _normalize_exon_id(self, feature: str) -> pd.DataFrame:
        """Normalize an exon ID or it's alias to one or more matching Ensembl exon ID."""
        featurel = [feature] + self.get_exon_alias(feature)

        return self._query_by_exon_id(featurel)

    def _normalize_gene_id(self, feature: str) -> pd.DataFrame:
        """Normalize a gene ID or it's alias to one or more matching Ensembl gene ID."""
        featurel = [feature] + self.get_gene_alias(feature)

        return self._query_by_gene_id(featurel)

    def _normalize_gene_name(self, feature: str) -> pd.DataFrame:
        """Normalize a gene ID or it's alias to one or more matching Ensembl gene ID."""
        featurel = [feature] + self.get_gene_alias(feature)

        return self._query_by_gene_name(featurel)

    def _normalize_protein_id(self, feature: str) -> pd.DataFrame:
        """Normalize a protein ID or it's alias to one or more matching Ensembl protein ID."""
        featurel = [feature] + self.get_protein_alias(feature)

        return self._query_by_protein_id(featurel)

    def _normalize_transcript_id(self, feature: str) -> pd.DataFrame:
        """Normalize a transcript ID or it's alias to one or more matching Ensembl transcript ID."""
        featurel = [feature] + self.get_transcript_alias(feature)

        return self._query_by_transcript_id(featurel)

    def _normalize_transcript_name(self, feature: str) -> pd.DataFrame:
        """Normalize a transcript name or it's alias to one or more matching Ensembl transcript ID."""
        featurel = [feature] + self.get_transcript_alias(feature)

        return self._query_by_transcript_name(featurel)

    # ---------------------------------------------------------------------------------------------
    # is_<feature>(feature)
    # ---------------------------------------------------------------------------------------------
    def is_contig(self, feature: str) -> bool:
        """Return True if the given feature is a contig."""
        return any((i[1] == CONTIG_ID for i in self.normalize_feature(feature, CONTIG_ID)))

    def is_exon(self, feature: str) -> bool:
        """Return True if the given feature is an exon."""
        return any((i[1] == EXON_ID for i in self.normalize_feature(feature, EXON_ID)))

    def is_gene(self, feature: str) -> bool:
        """Return True if the given feature is a gene."""
        return any((i[1] == GENE_ID for i in self.normalize_feature(feature, GENE_ID))) or any(
            (i[1] == GENE_NAME for i in self.normalize_feature(feature, GENE_NAME))
        )

    def is_protein(self, feature: str) -> bool:
        """Return True if the given feature is a protein."""
        return any((i[1] == PROTEIN_ID for i in self.normalize_feature(feature, PROTEIN_ID)))

    def is_transcript(self, feature: str) -> bool:
        """Return True if the given feature is a transcript."""
        return any(
            (i[1] == TRANSCRIPT_ID for i in self.normalize_feature(feature, TRANSCRIPT_ID))
        ) or any(
            (i[1] == TRANSCRIPT_NAME for i in self.normalize_feature(feature, TRANSCRIPT_NAME))
        )

    # ---------------------------------------------------------------------------------------------
    # <feature>_sequence
    # ---------------------------------------------------------------------------------------------
    def cdna_sequence(
        self,
        transcript_id: str,
        start: Optional[int] = None,
        end: Optional[int] = None,
        strand: Optional[str] = None,
    ) -> str:
        """Return the nucleotide sequence at the given cDNA coordinates."""
        return self._get_sequence(self.cdna, transcript_id, start, end, strand)

    def dna_sequence(
        self,
        contig_id: str,
        start: Optional[int] = None,
        end: Optional[int] = None,
        strand: str = "+",
    ) -> str:
        """Return the nucleotide sequence at the given contig coordinates."""
        return self._get_sequence(self.dna, contig_id, start, end, strand)

    def ncrna_sequence(
        self,
        transcript_id: str,
        start: Optional[int] = None,
        end: Optional[int] = None,
        strand: Optional[str] = None,
    ) -> str:
        """Return the nucleotide sequence at the given ncRNA coordinates."""
        return self._get_sequence(self.ncrna, transcript_id, start, end, strand)

    def peptide_sequence(
        self,
        protein_id: str,
        start: Optional[int] = None,
        end: Optional[int] = None,
        strand: Optional[str] = None,
    ) -> str:
        """Return the amino acid sequence at the given peptide coordinates."""
        return self._get_sequence(self.pep, protein_id, start, end, strand)

    def rna_sequence(
        self,
        transcript_id: str,
        start: Optional[int] = None,
        end: Optional[int] = None,
        strand: Optional[str] = None,
    ) -> str:
        """Return the nucleotide sequence at the given cDNA or ncRNA coordinates."""
        try:
            return self._get_sequence(self.cdna, transcript_id, start, end, strand)
        except KeyError:
            return self._get_sequence(self.ncrna, transcript_id, start, end, strand)

    def _get_sequence(
        self,
        fasta: Fasta,
        ref: str,
        start: Optional[int] = None,
        end: Optional[int] = None,
        strand: Optional[str] = None,
    ):
        """Return the sequence between the given positions (inclusive)."""
        seq = fasta[ref]
        seqlen = len(seq)

        # if no positions are given, return the whole sequence
        start = start if start is not None else 1
        end = end if end is not None else seqlen

        # validate that the given positions fall within the sequence
        if not (0 <= start < seqlen):
            raise ValueError(f"Start must be from 1 to {seqlen} ({start})")
        if not (1 < end <= seqlen):
            raise ValueError(f"End must be from 2 to {seqlen + 1} ({end})")

        # sanity check that the end position is after the start
        if end <= start:
            raise ValueError(f"End must be > start ({end} <= {start})")

        subseq = seq[start - 1 : end - 1]
        if strand == "-":
            subseq = reverse_complement(subseq)

        return subseq

    # ---------------------------------------------------------------------------------------------
    # get_<feature>_alias
    # ---------------------------------------------------------------------------------------------
    def get_contig_alias(self, feature: str) -> List[str]:
        """Return the aliases of the given contig (if any)."""
        return self._get_alias(feature, self.contig_alias)

    def get_exon_alias(self, feature: str) -> List[str]:
        """Return the aliases of the given exon (if any)."""
        return self._get_alias(feature, self.exon_alias)

    def get_gene_alias(self, feature: str) -> List[str]:
        """Return the aliases of the given gene (if any)."""
        return self._get_alias(feature, self.gene_alias)

    def get_protein_alias(self, feature: str) -> List[str]:
        """Return the aliases of the given protein (if any)."""
        return self._get_alias(feature, self.protein_alias)

    def get_transcript_alias(self, feature: str) -> List[str]:
        """Return the aliases of the given transcript (if any)."""
        return self._get_alias(feature, self.transcript_alias)

    def _get_alias(self, feature: str, source: Dict[str, List[str]]) -> List[str]:
        for key in [feature, strip_version(feature)]:
            if alias := source.get(key, []):
                return alias
        else:
            return []

    # ---------------------------------------------------------------------------------------------
    # canonical transcript
    # ---------------------------------------------------------------------------------------------
    def is_canonical_transcript(self, feature: str) -> bool:
        """Return True if the transcript is in the list of canonical transcripts."""
        return feature in self.canonical_transcript

    # ---------------------------------------------------------------------------------------------
    # get_<feature>
    # ---------------------------------------------------------------------------------------------
    def get_cdna(self, feature: str) -> List[CdnaPosition]:
        """Return the cDNA position(s) of the given feature."""
        result = []

        transcript_ids = self.transcript_ids(feature)
        mask = (self.df[TRANSCRIPT_ID].isin(transcript_ids)) & (self.df["feature"] == "cdna")
        for _, cdna in self.df[mask].iterrows():
            result.append(
                CdnaPosition(
                    _data=self,
                    contig_id=cdna.contig_id,
                    start=cdna.cdna_start,
                    end=cdna.cdna_end,
                    strand=cdna.strand,
                    gene_id=cdna.gene_id,
                    gene_name=cdna.gene_name,
                    transcript_id=cdna.transcript_id,
                    transcript_name=cdna.transcript_name,
                    protein_id=cdna.protein_id,
                )
            )

        return uniquify(result)

    def get_dna(self, feature: str) -> List[CdnaPosition]:
        """Return the DNA position(s) of the given feature."""
        result = []

        # get the strand of the original feature
        _, df = self._normalize_feature(feature)
        strand_list = self._uniquify_series(df["strand"])

        for contig_id in self.contig_ids(feature):
            contig_seq = self.dna[contig_id]
            start = 1
            end = len(contig_seq)
            for strand in strand_list:
                result.append(
                    DnaPosition(
                        _data=self, contig_id=contig_id, start=start, end=end, strand=strand
                    )
                )

        return uniquify(result)

    def get_exons(self, feature: str) -> List[ExonPosition]:
        """Return the gene position(s) of the given feature."""
        result = []

        exon_ids = self.exon_ids(feature)
        mask = (self.df[EXON_ID].isin(exon_ids)) & (self.df["feature"] == "exon")
        for _, exon in self.df[mask].iterrows():
            result.append(
                ExonPosition(
                    _data=self,
                    contig_id=exon.contig_id,
                    start=int(exon.exon_number),
                    end=int(exon.exon_number),
                    strand=exon.strand,
                    gene_id=exon.gene_id,
                    gene_name=exon.gene_name,
                    transcript_id=exon.transcript_id,
                    transcript_name=exon.transcript_name,
                    exon_id=exon.exon_id,
                )
            )

        return uniquify(result)

    def get_genes(self, feature: str) -> List[DnaPosition]:
        """Return the gene position(s) of the given feature."""
        result = []

        gene_ids = self.gene_ids(feature)
        mask = (self.df[GENE_ID].isin(gene_ids)) & (self.df["feature"] == "gene")
        for _, gene in self.df[mask].iterrows():
            result.append(
                DnaPosition(
                    _data=self,
                    contig_id=gene.contig_id,
                    start=gene.start,
                    end=gene.end,
                    strand=gene.strand,
                )
            )

        return uniquify(result)

    def get_transcripts(self, feature: str) -> List[RnaPosition]:
        """Return the transcript position(s) of the given feature."""
        result = []

        transcript_ids = self.transcript_ids(feature)
        mask = (self.df[TRANSCRIPT_ID].isin(transcript_ids)) & (self.df["feature"] == "transcript")
        for _, transcript in self.df[mask].iterrows():
            result.append(
                RnaPosition(
                    _data=self,
                    contig_id=transcript.contig_id,
                    start=transcript.transcript_start,
                    end=transcript.transcript_end,
                    strand=transcript.strand,
                    gene_id=transcript.gene_id,
                    gene_name=transcript.gene_name,
                    transcript_id=transcript.transcript_id,
                    transcript_name=transcript.transcript_name,
                )
            )

        return uniquify(result)

    # ---------------------------------------------------------------------------------------------
    # <feature>_to_<feature>
    # ---------------------------------------------------------------------------------------------
    def cdna_to_cdna(
        self, feature: str, start: int, end: Optional[int] = None, strand: Optional[str] = None
    ) -> List[CdnaPosition]:
        """Map a cDNA position to a cDNA position."""
        return self._map(
            self.transcript_ids, self._cdna_to_cdna, feature, start, end=end, strand=strand
        )

    def cdna_to_dna(
        self, feature: str, start: int, end: Optional[int] = None, strand: Optional[str] = None
    ) -> List[DnaPosition]:
        """Map a cDNA position to a DNA position."""
        return self._map(
            self.transcript_ids, self._cdna_to_dna, feature, start, end=end, strand=strand
        )

    def cdna_to_exon(
        self, feature: str, start: int, end: Optional[int] = None, strand: Optional[str] = None
    ) -> List[ExonPosition]:
        """Map a cDNA position to an exon position."""
        return self._map(
            self.transcript_ids, self._cdna_to_exon, feature, start, end=end, strand=strand
        )

    def cdna_to_protein(
        self, feature: str, start: int, end: Optional[int] = None, strand: Optional[str] = None
    ) -> List[ProteinPosition]:
        """Map a cDNA position to a protein position."""
        return self._map(
            self.transcript_ids, self._cdna_to_protein, feature, start, end=end, strand=strand
        )

    def cdna_to_rna(
        self, feature: str, start: int, end: Optional[int] = None, strand: Optional[str] = None
    ) -> List[RnaPosition]:
        """Map a cDNA position to a RNA position."""
        return self._map(
            self.transcript_ids, self._cdna_to_rna, feature, start, end=end, strand=strand
        )

    def dna_to_cdna(
        self, feature: str, start: int, end: Optional[int] = None, strand: Optional[str] = None
    ) -> List[CdnaPosition]:
        """Map a DNA position to a cDNA position."""
        return self._map(self.contig_ids, self._dna_to_cdna, feature, start, end=end, strand=strand)

    def dna_to_dna(
        self, feature: str, start: int, end: Optional[int] = None, strand: Optional[str] = None
    ) -> List[DnaPosition]:
        """Map a DNA position to a DNA position."""
        return self._map(
            self.contig_ids, self._dna_to_dna, feature, start=start, end=end, strand=strand
        )

    def dna_to_exon(
        self, feature: str, start: int, end: Optional[int] = None, strand: Optional[str] = None
    ) -> List[ExonPosition]:
        """Map a DNA position to an exon position."""
        return self._map(
            self.contig_ids, self._dna_to_exon, feature, start=start, end=end, strand=strand
        )

    def dna_to_protein(
        self, feature: str, start: int, end: Optional[int] = None, strand: Optional[str] = None
    ) -> List[ProteinPosition]:
        """Map a DNA position to a protein position."""
        return self._map(
            self.contig_ids, self._dna_to_protein, feature, start=start, end=end, strand=strand
        )

    def dna_to_rna(
        self, feature: str, start: int, end: Optional[int] = None, strand: Optional[str] = None
    ) -> List[RnaPosition]:
        """Map a DNA position to a RNA position."""
        return self._map(
            self.contig_ids, self._dna_to_rna, feature, start=start, end=end, strand=strand
        )

    def exon_to_cdna(
        self,
        feature: str,
        start: Optional[int] = None,
        end: Optional[int] = None,
        strand: Optional[str] = None,
    ) -> List[CdnaPosition]:
        """Map an exon position to a cDNA position."""
        return self._map_exon(self._exon_to_cdna, feature, start, end, strand)

    def exon_to_dna(
        self,
        feature: str,
        start: Optional[int] = None,
        end: Optional[int] = None,
        strand: Optional[str] = None,
    ) -> List[DnaPosition]:
        """Map an exon position to a DNA position."""
        return self._map_exon(self._exon_to_dna, feature, start, end, strand)

    def exon_to_exon(
        self,
        feature: str,
        start: Optional[int] = None,
        end: Optional[int] = None,
        strand: Optional[str] = None,
    ) -> List[ExonPosition]:
        """Map an exon position to an exon position."""
        return self._map_exon(self._exon_to_exon, feature, start, end, strand)

    def exon_to_protein(
        self,
        feature: str,
        start: Optional[int] = None,
        end: Optional[int] = None,
        strand: Optional[str] = None,
    ) -> List[ProteinPosition]:
        """Map an exon position to a protein position."""
        return self._map_exon(self._exon_to_protein, feature, start, end, strand)

    def exon_to_rna(
        self,
        feature: str,
        start: Optional[int] = None,
        end: Optional[int] = None,
        strand: Optional[str] = None,
    ) -> List[RnaPosition]:
        """Map an exon position to a RNA position."""
        return self._map_exon(self._exon_to_rna, feature, start, end, strand)

    def protein_to_cdna(
        self, feature: str, start: int, end: Optional[int] = None, strand: Optional[str] = None
    ) -> List[CdnaPosition]:
        """Map a protein position to a cDNA position."""
        return self._map(
            self.transcript_ids, self._protein_to_cdna, feature, start, end=end, strand=strand
        )

    def protein_to_dna(
        self, feature: str, start: int, end: Optional[int] = None, strand: Optional[str] = None
    ) -> List[DnaPosition]:
        """Map a protein position to a DNA position."""
        return self._map(
            self.transcript_ids, self._protein_to_dna, feature, start, end=end, strand=strand
        )

    def protein_to_exon(
        self, feature: str, start: int, end: Optional[int] = None, strand: Optional[str] = None
    ) -> List[ExonPosition]:
        """Map a protein position to an exon position."""
        return self._map(
            self.transcript_ids, self._protein_to_exon, feature, start, end=end, strand=strand
        )

    def protein_to_protein(
        self, feature: str, start: int, end: Optional[int] = None, strand: Optional[str] = None
    ) -> List[ProteinPosition]:
        """Map a protein position to a protein position."""
        return self._map(
            self.transcript_ids, self._protein_to_protein, feature, start, end=end, strand=strand
        )

    def protein_to_rna(
        self, feature: str, start: int, end: Optional[int] = None, strand: Optional[str] = None
    ) -> List[RnaPosition]:
        """Map a protein position to a RNA position."""
        return self._map(
            self.transcript_ids, self._protein_to_rna, feature, start, end=end, strand=strand
        )

    def rna_to_cdna(
        self, feature: str, start: int, end: Optional[int] = None, strand: Optional[str] = None
    ) -> List[CdnaPosition]:
        """Map a RNA position to a cDNA position."""
        return self._map(
            self.transcript_ids, self._rna_to_cdna, feature, start, end=end, strand=strand
        )

    def rna_to_dna(
        self, feature: str, start: int, end: Optional[int] = None, strand: Optional[str] = None
    ) -> List[DnaPosition]:
        """Map a RNA position to a DNA position."""
        return self._map(
            self.transcript_ids, self._rna_to_dna, feature, start, end=end, strand=strand
        )

    def rna_to_exon(
        self, feature: str, start: int, end: Optional[int] = None, strand: Optional[str] = None
    ) -> List[ExonPosition]:
        """Map a RNA position to an exon position."""
        return self._map(
            self.transcript_ids, self._rna_to_exon, feature, start, end=end, strand=strand
        )

    def rna_to_protein(
        self, feature: str, start: int, end: Optional[int] = None, strand: Optional[str] = None
    ) -> List[ProteinPosition]:
        """Map a RNA position to a protein position."""
        return self._map(
            self.transcript_ids, self._rna_to_protein, feature, start, end=end, strand=strand
        )

    def rna_to_rna(
        self, feature: str, start: int, end: Optional[int] = None, strand: Optional[str] = None
    ) -> List[RnaPosition]:
        """Map a RNA position to a RNA position."""
        return self._map(
            self.transcript_ids, self._rna_to_rna, feature, start, end=end, strand=strand
        )

    def _cdna_to_cdna(
        self,
        transcript_ids: List[str],
        start: int,
        end: int,
        strand: List[str],
        include_stop: bool = True,
    ) -> List[CdnaPosition]:
        feature = ["CDS", "stop_codon"] if include_stop else ["CDS"]

        def convert(position: int):
            result = []

            mask = (
                (self.df[TRANSCRIPT_ID].isin(transcript_ids))
                & (self.df["cdna_start"] <= position)
                & (self.df["cdna_end"] >= position)
                & (self.df["strand"].isin(strand))
                & (self.df["feature"].isin(feature))
            )
            for _, cds in self.df[mask].iterrows():
                result.append(
                    CdnaPosition(
                        _data=self,
                        contig_id=cds.contig_id,
                        start=start,
                        end=end,
                        strand=cds.strand,
                        gene_id=cds.gene_id,
                        gene_name=cds.gene_name,
                        transcript_id=cds.transcript_id,
                        transcript_name=cds.transcript_name,
                        protein_id=cds.protein_id,
                    )
                )

            return result

        result_start = convert(start)
        if start == end:
            return result_start
        else:
            result_end = convert(end)
            return merge_positions(result_start, result_end, TRANSCRIPT_ID)

    def _cdna_to_dna(
        self,
        transcript_ids: List[str],
        start: int,
        end: int,
        strand: List[str],
        include_stop: bool = True,
    ) -> List[DnaPosition]:
        feature = ["CDS", "stop_codon"] if include_stop else ["CDS"]

        def convert(position: int):
            result = []

            mask = (
                (self.df[TRANSCRIPT_ID].isin(transcript_ids))
                & (self.df["cdna_start"] <= position)
                & (self.df["cdna_end"] >= position)
                & (self.df["strand"].isin(strand))
                & (self.df["feature"].isin(feature))
            )
            for _, cds in self.df[mask].iterrows():
                offset = position - cds.cdna_start
                if cds.strand == "-":
                    new_start = new_end = cds.end - offset
                else:
                    new_start = new_end = cds.start + offset

                result.append(
                    DnaPosition(
                        _data=self,
                        contig_id=cds.contig_id,
                        start=new_start,
                        end=new_end,
                        strand=cds.strand,
                    )
                )

            return result

        result_start = convert(start)
        if start == end:
            return result_start
        else:
            result_end = convert(end)
            return merge_positions(result_start, result_end, CONTIG_ID)

    def _cdna_to_exon(
        self,
        transcript_ids: List[str],
        start: int,
        end: int,
        strand: List[str],
        include_stop: bool = True,
    ) -> List[ExonPosition]:
        feature = ["CDS", "stop_codon"] if include_stop else ["CDS"]

        def convert(position: int):
            result = []

            mask_cds = (
                (self.df[TRANSCRIPT_ID].isin(transcript_ids))
                & (self.df["cdna_start"] <= position)
                & (self.df["cdna_end"] >= position)
                & (self.df["strand"].isin(strand))
                & (self.df["feature"].isin(feature))
            )
            for _, cds in self.df[mask_cds].iterrows():
                mask_exon = (
                    (self.df[TRANSCRIPT_ID].isin(transcript_ids))
                    & (self.df["exon_number"] == cds.exon_number)
                    & (self.df["feature"] == "exon")
                )
                for _, exon in self.df[mask_exon].iterrows():
                    result.append(
                        ExonPosition(
                            _data=self,
                            contig_id=exon.contig_id,
                            start=int(exon.exon_number),
                            end=int(exon.exon_number),
                            strand=exon.strand,
                            gene_id=exon.gene_id,
                            gene_name=exon.gene_name,
                            transcript_id=exon.transcript_id,
                            transcript_name=exon.transcript_name,
                            exon_id=exon.exon_id,
                        )
                    )

            return result

        result_start = convert(start)
        result_end = convert(end)
        assert result_start == result_end  # TODO: mapping across introns?

        return result_start

    def _cdna_to_protein(
        self, transcript_ids: List[str], start: int, end: int, strand: List[str]
    ) -> List[ProteinPosition]:
        def convert(x):
            return floor((x - 1) / 3 + 1)

        protein = []
        for cdna in self._cdna_to_cdna(transcript_ids, start, end, strand, include_stop=False):
            pstart = convert(start)
            pend = convert(end)
            protein.append(ProteinPosition.copy_from(cdna, start=pstart, end=pend))

        return protein

    def _cdna_to_rna(
        self,
        transcript_ids: List[str],
        start: int,
        end: int,
        strand: List[str],
        include_stop: bool = True,
    ) -> List[RnaPosition]:
        feature = ["CDS", "stop_codon"] if include_stop else ["CDS"]

        def convert(position: int):
            result = []

            mask = (
                (self.df[TRANSCRIPT_ID].isin(transcript_ids))
                & (self.df["cdna_start"] <= position)
                & (self.df["cdna_end"] >= position)
                & (self.df["strand"].isin(strand))
                & (self.df["feature"].isin(feature))
            )
            for _, cds in self.df[mask].iterrows():
                offset = position - cds.cdna_start
                new_start = new_end = cds.transcript_start + offset
                result.append(
                    RnaPosition(
                        _data=self,
                        contig_id=cds.contig_id,
                        start=new_start,
                        end=new_end,
                        strand=cds.strand,
                        gene_id=cds.gene_id,
                        gene_name=cds.gene_name,
                        transcript_id=cds.transcript_id,
                        transcript_name=cds.transcript_name,
                    )
                )

            return result

        result_start = convert(start)
        if start == end:
            return result_start
        else:
            result_end = convert(end)
            return merge_positions(result_start, result_end, TRANSCRIPT_ID)

    def _dna_to_cdna(
        self,
        contig_ids: List[str],
        start: int,
        end: int,
        strand: List[str],
        include_stop: bool = True,
    ) -> List[CdnaPosition]:
        feature = ["CDS", "stop_codon"] if include_stop else ["CDS"]

        def convert(position: int):
            result = []

            mask = (
                (self.df[CONTIG_ID].isin(contig_ids))
                & (self.df["start"] <= position)
                & (self.df["end"] >= position)
                & (self.df["strand"].isin(strand))
                & (self.df["feature"].isin(feature))
            )
            for _, cds in self.df[mask].iterrows():
                if cds.strand == "-":
                    new_start = new_end = cds.end - position + cds.cdna_start
                else:
                    new_start = new_end = position - cds.start + cds.cdna_start

                result.append(
                    CdnaPosition(
                        _data=self,
                        contig_id=cds.contig_id,
                        start=new_start,
                        end=new_end,
                        strand=cds.strand,
                        gene_id=cds.gene_id,
                        gene_name=cds.gene_name,
                        transcript_id=cds.transcript_id,
                        transcript_name=cds.transcript_name,
                        protein_id=cds.protein_id,
                    )
                )

            return result

        result_start = convert(start)
        if start == end:
            return result_start
        else:
            result_end = convert(end)
            return merge_positions(result_start, result_end, TRANSCRIPT_ID)

    def _dna_to_dna(
        self, contig_ids: List[str], start: int, end: int, strand: List[str]
    ) -> List[DnaPosition]:
        result = []
        for c, s in product(contig_ids, strand):
            result.append(DnaPosition(self, c, start, end, s))

        return result

    def _dna_to_exon(
        self, contig_ids: List[str], start: int, end: int, strand: List[str]
    ) -> List[ExonPosition]:
        def convert(position: int):
            result = []

            mask = (
                (self.df[CONTIG_ID].isin(contig_ids))
                & (self.df["start"] <= position)
                & (self.df["end"] >= position)
                & (self.df["strand"].isin(strand))
                & (self.df["feature"].isin(["exon"]))
            )
            for _, exon in self.df[mask].iterrows():
                result.append(
                    ExonPosition(
                        _data=self,
                        contig_id=exon.contig_id,
                        start=int(exon.exon_number),
                        end=int(exon.exon_number),
                        strand=exon.strand,
                        gene_id=exon.gene_id,
                        gene_name=exon.gene_name,
                        transcript_id=exon.transcript_id,
                        transcript_name=exon.transcript_name,
                        exon_id=exon.exon_id,
                    )
                )

            return result

        result_start = convert(start)
        result_end = convert(end)
        assert result_start == result_end  # TODO: mapping across introns?

        return result_start

    def _dna_to_protein(
        self, contig_ids: List[str], start: int, end: int, strand: List[str]
    ) -> List[ProteinPosition]:
        def convert(x):
            return floor((x - 1) / 3 + 1)

        protein = []
        for cdna in self._dna_to_cdna(contig_ids, start, end, strand, include_stop=False):
            pstart = convert(cdna.start)
            pend = convert(cdna.end)
            protein.append(ProteinPosition.copy_from(cdna, start=pstart, end=pend, _data=self))

        return protein

    def _dna_to_rna(
        self, contig_ids: List[str], start: int, end: int, strand: List[str]
    ) -> List[RnaPosition]:
        def convert(position: int):
            result = []

            mask = (
                (self.df[CONTIG_ID].isin(contig_ids))
                & (self.df["start"] <= position)
                & (self.df["end"] >= position)
                & (self.df["strand"].isin(strand))
                & (self.df["feature"] == "exon")
            )
            for _, exon in self.df[mask].iterrows():
                if exon.strand == "-":
                    new_start = new_end = exon.end - position + exon.transcript_start
                else:
                    new_start = new_end = position - exon.start + exon.transcript_start

                result.append(
                    RnaPosition(
                        _data=self,
                        contig_id=exon.contig_id,
                        start=new_start,
                        end=new_end,
                        strand=exon.strand,
                        gene_id=exon.gene_id,
                        gene_name=exon.gene_name,
                        transcript_id=exon.transcript_id,
                        transcript_name=exon.transcript_name,
                    )
                )

            return result

        result_start = convert(start)
        if start == end:
            return result_start
        else:
            result_end = convert(end)
            return merge_positions(result_start, result_end, TRANSCRIPT_ID)

    def _exon_to_cdna(
        self,
        transcript_ids: List[str],
        start: int,
        end: int,
        strand: List[str],
        include_stop: bool = True,
    ) -> List[CdnaPosition]:
        feature = ["CDS", "stop_codon"] if include_stop else ["CDS"]

        def convert(position: int):
            result = []

            mask = (
                (self.df[TRANSCRIPT_ID].isin(transcript_ids))
                & (self.df["exon_number"] == str(position))  # TODO: exon number should be int
                & (self.df["strand"].isin(strand))
                & (self.df["feature"].isin(feature))
            )
            for _, cds in self.df[mask].iterrows():
                result.append(
                    CdnaPosition(
                        _data=self,
                        contig_id=cds.contig_id,
                        start=cds.cdna_start,
                        end=cds.cdna_end,
                        strand=cds.strand,
                        gene_id=cds.gene_id,
                        gene_name=cds.gene_name,
                        transcript_id=cds.transcript_id,
                        transcript_name=cds.transcript_name,
                        protein_id=cds.protein_id,
                    )
                )

            return result

        result_start = convert(start)
        if start == end:
            return result_start
        else:
            result_end = convert(end)
            return merge_positions(result_start, result_end, TRANSCRIPT_ID)

    def _exon_to_dna(
        self, transcript_ids: List[str], start: int, end: int, strand: List[str]
    ) -> List[DnaPosition]:
        def convert(position: int):
            result = []

            mask = (
                (self.df[TRANSCRIPT_ID].isin(transcript_ids))
                # TODO: exon number should be int
                & (self.df["exon_number"] == str(position))
                & (self.df["strand"].isin(strand))
                & (self.df["feature"] == "exon")
            )
            for _, exon in self.df[mask].iterrows():
                result.append(
                    DnaPosition(
                        _data=self,
                        contig_id=exon.contig_id,
                        start=exon.start,
                        end=exon.end,
                        strand=exon.strand,
                    )
                )

            return result

        result_start = convert(start)
        if start == end:
            return result_start
        else:
            result_end = convert(end)
            return merge_positions(result_start, result_end, CONTIG_ID)

    def _exon_to_exon(
        self, transcript_ids: List[str], start: int, end: int, strand: List[str]
    ) -> List[ExonPosition]:
        def convert(position: int):
            result = []

            mask = (
                (self.df[TRANSCRIPT_ID].isin(transcript_ids))
                # TODO: exon number should be int
                & (self.df["exon_number"] == str(position))
                & (self.df["strand"].isin(strand))
                & (self.df["feature"] == "exon")
            )
            for _, exon in self.df[mask].iterrows():
                result.append(
                    ExonPosition(
                        _data=self,
                        contig_id=exon.contig_id,
                        start=int(exon.exon_number),
                        end=int(exon.exon_number),
                        strand=exon.strand,
                        gene_id=exon.gene_id,
                        gene_name=exon.gene_name,
                        transcript_id=exon.transcript_id,
                        transcript_name=exon.transcript_name,
                        exon_id=exon.exon_id,
                    )
                )

            return result

        result_start = convert(start)
        if start == end:
            return result_start
        else:
            result_end = convert(end)
            return merge_positions(result_start, result_end, TRANSCRIPT_ID)

    def _exon_to_protein(
        self, transcript_ids: List[str], start: int, end: int, strand: List[str]
    ) -> List[ProteinPosition]:
        result = []
        for cdna in self._exon_to_cdna(transcript_ids, start, end, strand, include_stop=False):
            result.extend(cdna.to_protein())

        return result

    def _exon_to_rna(
        self, transcript_ids: List[str], start: int, end: int, strand: List[str]
    ) -> List[RnaPosition]:
        def convert(position: int):
            result = []

            mask = (
                (self.df[TRANSCRIPT_ID].isin(transcript_ids))
                & (self.df["exon_number"] == str(position))  # TODO: exon number should be int
                & (self.df["strand"].isin(strand))
                & (self.df["feature"] == "exon")
            )
            for _, exon in self.df[mask].iterrows():
                result.append(
                    RnaPosition(
                        _data=self,
                        contig_id=exon.contig_id,
                        start=exon.transcript_start,
                        end=exon.transcript_end,
                        strand=exon.strand,
                        gene_id=exon.gene_id,
                        gene_name=exon.gene_name,
                        transcript_id=exon.transcript_id,
                        transcript_name=exon.transcript_name,
                    )
                )

            return result

        result_start = convert(start)
        if start == end:
            return result_start
        else:
            result_end = convert(end)
            return merge_positions(result_start, result_end, TRANSCRIPT_ID)

    def _protein_to_cdna(
        self, transcript_ids: List[str], start: int, end: int, strand: List[str]
    ) -> List[CdnaPosition]:
        def convert(position: int):
            return ((position - 1) * 3) + 1

        cdna_start = convert(start)
        cdna_end = convert(end) + 2

        return self._cdna_to_cdna(transcript_ids, cdna_start, cdna_end, strand, include_stop=False)

    def _protein_to_dna(
        self, transcript_ids: List[str], start: int, end: int, strand: List[str]
    ) -> List[DnaPosition]:
        result = []
        for cdna in self._protein_to_cdna(transcript_ids, start, end, strand):
            result.extend(cdna.to_dna())

        return result

    def _protein_to_exon(
        self, transcript_ids: List[str], start: int, end: int, strand: List[str]
    ) -> List[ExonPosition]:
        result = []
        for cdna in self._protein_to_cdna(transcript_ids, start, end, strand):
            result.extend(cdna.to_exon())

        return result

    def _protein_to_protein(
        self, transcript_ids: List[str], start: int, end: int, strand: List[str]
    ) -> List[ProteinPosition]:
        result = []
        for cdna in self._protein_to_cdna(transcript_ids, start, end, strand):
            result.extend(cdna.to_protein())

        return result

    def _protein_to_rna(
        self, transcript_ids: List[str], start: int, end: int, strand: List[str]
    ) -> List[RnaPosition]:
        result = []
        for cdna in self._protein_to_cdna(transcript_ids, start, end, strand):
            result.extend(cdna.to_rna())

        return result

    def _rna_to_cdna(
        self,
        transcript_ids: List[str],
        start: int,
        end: int,
        strand: List[str],
        include_stop: bool = True,
    ) -> List[CdnaPosition]:
        feature = ["CDS", "stop_codon"] if include_stop else ["CDS"]

        def convert(position: int):
            result = []

            mask = (
                (self.df[TRANSCRIPT_ID].isin(transcript_ids))
                & (self.df["transcript_start"] <= position)
                & (self.df["transcript_end"] >= position)
                & (self.df["strand"].isin(strand))
                & (self.df["feature"].isin(feature))
            )
            for _, cds in self.df[mask].iterrows():
                offset = position - cds.transcript_start
                new_start = new_end = cds.cdna_start + offset
                result.append(
                    CdnaPosition(
                        _data=self,
                        contig_id=cds.contig_id,
                        start=new_start,
                        end=new_end,
                        strand=cds.strand,
                        gene_id=cds.gene_id,
                        gene_name=cds.gene_name,
                        transcript_id=cds.transcript_id,
                        transcript_name=cds.transcript_name,
                        protein_id=cds.protein_id,
                    )
                )

            return result

        result_start = convert(start)
        if start == end:
            return result_start
        else:
            result_end = convert(end)
            return merge_positions(result_start, result_end, TRANSCRIPT_ID)

    def _rna_to_dna(
        self, transcript_ids: List[str], start: int, end: int, strand: List[str]
    ) -> List[DnaPosition]:
        def convert(position: int):
            result = []

            mask = (
                (self.df[TRANSCRIPT_ID].isin(transcript_ids))
                & (self.df["transcript_start"] <= position)
                & (self.df["transcript_end"] >= position)
                & (self.df["strand"].isin(strand))
                & (self.df["feature"].isin(["exon"]))
            )
            exon_df = self.df[mask]
            for _, exon in exon_df.iterrows():
                offset = position - exon.transcript_start
                if exon.strand == "-":
                    new_start = new_end = exon.end - offset
                else:
                    new_start = new_end = exon.start + offset

                result.append(
                    DnaPosition(
                        _data=self,
                        contig_id=exon.contig_id,
                        start=new_start,
                        end=new_end,
                        strand=exon.strand,
                    )
                )

            return result

        result_start = convert(start)
        if start == end:
            return result_start
        else:
            result_end = convert(end)
            print(result_end)
            return merge_positions(result_start, result_end, CONTIG_ID)

    def _rna_to_exon(
        self, transcript_ids: List[str], start: int, end: int, strand: List[str]
    ) -> List[ExonPosition]:
        def convert(position: int):
            result = []

            mask = (
                (self.df[TRANSCRIPT_ID].isin(transcript_ids))
                & (self.df["transcript_start"] <= position)
                & (self.df["transcript_end"] >= position)
                & (self.df["strand"].isin(strand))
                & (self.df["feature"].isin(["exon"]))
            )
            for _, exon in self.df[mask].iterrows():
                result.append(
                    ExonPosition(
                        _data=self,
                        contig_id=exon.contig_id,
                        start=int(exon.exon_number),
                        end=int(exon.exon_number),
                        strand=exon.strand,
                        gene_id=exon.gene_id,
                        gene_name=exon.gene_name,
                        transcript_id=exon.transcript_id,
                        transcript_name=exon.transcript_name,
                        exon_id=exon.exon_id,
                    )
                )

            return result

        result_start = convert(start)
        if start == end:
            return result_start
        else:
            result_end = convert(end)
            return merge_positions(result_start, result_end, TRANSCRIPT_ID)

    def _rna_to_protein(
        self, transcript_ids: List[str], start: int, end: int, strand: List[str]
    ) -> List[ProteinPosition]:
        result = []
        for cdna in self._rna_to_cdna(transcript_ids, start, end, strand, include_stop=False):
            result.extend(cdna.to_protein())

        return result

    def _rna_to_rna(
        self, transcript_ids: List[str], start: int, end: int, strand: List[str]
    ) -> List[RnaPosition]:
        def convert(position: int):
            result = []

            mask = (
                (self.df[TRANSCRIPT_ID].isin(transcript_ids))
                & (self.df["transcript_start"] <= position)
                & (self.df["transcript_end"] >= position)
                & (self.df["strand"].isin(strand))
                & (self.df["feature"] == "exon")
            )
            for _, exon in self.df[mask].iterrows():
                result.append(
                    RnaPosition(
                        _data=self,
                        contig_id=exon.contig_id,
                        start=start,
                        end=end,
                        strand=exon.strand,
                        gene_id=exon.gene_id,
                        gene_name=exon.gene_name,
                        transcript_id=exon.transcript_id,
                        transcript_name=exon.transcript_name,
                    )
                )

            return result

        result_start = convert(start)
        if start == end:
            return result_start
        else:
            result_end = convert(end)
            return merge_positions(result_start, result_end, TRANSCRIPT_ID)

    def _map(
        self,
        idfunc: Callable,
        mapfunc: Callable,
        feature: str,
        start: int,
        end: Optional[int] = None,
        strand: Optional[str] = None,
    ) -> List:
        end = end if end is not None else start
        strandl = [strand] if strand is not None else ["+", "-"]
        featurel = idfunc(feature)
        result = mapfunc(featurel, start, end, strandl)

        return uniquify(result)

    def _map_exon(
        self,
        function: Callable,
        feature: str,
        start: Optional[int],
        end: Optional[int],
        strand: Optional[str],
    ) -> List:
        if start is None and end is not None:
            start = end
        if end is None and start is not None:
            end = start

        if start is None:
            positions = []
            for exon_id in self.exon_ids(feature):
                positions.extend(self._get_exon_number(exon_id))
        else:
            strandl = [strand] if strand is not None else ["+", "-"]
            positions = product(self.transcript_ids(feature), [start], [end], strandl)  # type: ignore

        result = []
        for transcript_id, start2, end2, strand2 in positions:
            result.extend(function([transcript_id], start2, end2, [strand2]))

        return uniquify(result)

    def _get_exon_number(self, exon_id: str) -> List[Tuple[str, int, int, str]]:
        result = []

        mask = (self.df[EXON_ID] == exon_id) & (self.df["feature"] == "exon")
        for _, exon in self.df[mask].iterrows():
            result.append((exon.transcript_id, exon.exon_number, exon.exon_number, exon.strand))

        return result


class EnsemblRelease(Core):
    def __init__(
        self,
        species: str = DEFAULT_SPECIES,
        release: int = DEFAULT_RELEASE,
        cache_dir: str = "",
        canonical_transcript: str = "",
        contig_alias: str = "",
        exon_alias: str = "",
        gene_alias: str = "",
        protein_alias: str = "",
        transcript_alias: str = "",
    ):
        """Load annotations for the given release."""
        self.ensembl_cache = EnsemblCache(species, release, cache_dir=cache_dir)

        self.cache_dir = self.ensembl_cache.release_cache_dir
        self.reference = self.ensembl_cache.reference
        self.release = self.ensembl_cache.release
        self.species = self.ensembl_cache.species

        self.df = self.ensembl_cache.load_df()
        self.cdna = self.ensembl_cache.load_cdna_fasta()
        self.dna = self.ensembl_cache.load_dna_fasta()
        self.pep = self.ensembl_cache.load_pep_fasta()
        self.ncrna = self.ensembl_cache.load_ncrna_fasta()

        self.canonical_transcript = _parse_txt_to_list(canonical_transcript, "canonical transcript")
        self.contig_alias = _parse_tsv_to_dict(contig_alias, "contig aliases")
        self.exon_alias = _parse_tsv_to_dict(exon_alias, "exon aliases")
        self.gene_alias = _parse_tsv_to_dict(gene_alias, "gene aliases")
        self.protein_alias = _parse_tsv_to_dict(protein_alias, "protein aliases")
        self.transcript_alias = _parse_tsv_to_dict(transcript_alias, "transcript aliases")

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}(species={self.species}, release={self.release})"


def instance():
    """Return the last called instance of `EnsemblRelease` otherwise initialize an instance with
    the default parameters and return that.
    """
    if EnsemblRelease._current is None:
        EnsemblRelease._current = EnsemblRelease()

    return EnsemblRelease._current


def _parse_txt_to_list(path: str, message: str) -> List[str]:
    """Parse a text file into a list of unique values."""
    result = set()

    if path:
        logger.debug(f"Loading {message} from {path}")
        with open(path, "r") as fh:
            for line in fh:
                if line:
                    result.add(line.strip())
    else:
        logger.debug(f"No {message} to load")

    return sorted(result)


def _parse_tsv_to_dict(path: str, message: str) -> Dict[str, List[str]]:
    """Parse a TSV of one-to-one alias-to-name mappings."""
    result: Dict = {}

    if path:
        # load data from the TSV file, keep only unique values
        logger.debug(f"Loading {message} from {path}")
        with open(path, "r") as fh:
            for line in fh:
                if line:
                    alias, name = line.strip().split("\t")[:2]
                    result.setdefault(alias, set())
                    result[alias].add(name)

        # convert values to sorted lists
        for alias in result:
            result[alias] = sorted(result[alias])
    else:
        logger.debug(f"No {message} to load")

    return result


def merge_positions(
    start_positions: List[_Position], end_positions: List[_Position], key: str
) -> List:
    result = set()

    for start_pos, end_pos in product(start_positions, end_positions):
        start_key = getattr(start_pos, key)
        end_key = getattr(end_pos, key)
        if start_key == end_key:
            assert start_pos.__class__ == end_pos.__class__
            kwargs = start_pos.asdict()
            start = min((start_pos.start, end_pos.start))
            end = max((start_pos.end, end_pos.end))
            kwargs.update({"start": start, "end": end})
            new_position = start_pos.__class__(**kwargs)
            result.add(new_position)

    return sorted(result)


def uniquify(positions: List) -> List:
    """Return a list of only unique positions."""
    return sorted(set(positions))
