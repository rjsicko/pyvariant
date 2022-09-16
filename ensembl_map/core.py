from __future__ import annotations

from dataclasses import dataclass, fields
from itertools import product
from math import floor
from typing import Any, Callable, Dict, List, Optional, Tuple, Union

import pandas
from logzero import logger
from pyfaidx import Fasta

from .cache import Cache
from .constants import (
    CDNA,
    CONTIG,
    DEFAULT_REFERENCE,
    DEFAULT_RELEASE,
    DEFAULT_SPECIES,
    EXON,
    GENE,
    PROTEIN,
    TRANSCRIPT,
)
from .utils import is_ensembl_id, reverse_complement, strip_version


# -------------------------------------------------------------------------------------------------
# position classes
# -------------------------------------------------------------------------------------------------
@dataclass(eq=True, frozen=True, order=True)
class _Position:
    _data: EnsemblRelease
    contig_id: str
    start: int
    end: int
    strand: str

    @classmethod
    def copy_from(cls, position: _Position, **kwargs):
        return cls(**{**position.asdict(), **kwargs})

    def __getitem__(self, item):
        return getattr(self, item)

    # def __repr__(self) -> str:
    #     return str(self)

    @property
    def length(self) -> int:
        return self.end - self.start + 1

    def asdict(self) -> Dict[str, Any]:
        return {f.name: self[f.name] for f in fields(self)}

    def on_negative_strand(self) -> bool:
        return self.strand == "-"

    def on_positive_strand(self) -> bool:
        return self.strand == "+"

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

    # def __str__(self) -> str:
    #     if self.start == self.end:
    #         return f"{self.transcript_id}:c.{self.start}"
    #     else:
    #         return f"{self.transcript_id}:c.{self.start}-{self.end}"

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
    # def __str__(self) -> str:
    #     if self.start == self.end:
    #         return f"{self.contig_id}:g.{self.start}"
    #     else:
    #         return f"{self.contig_id}:g.{self.start}-{self.end}"

    @classmethod
    def load(cls, row: pandas.Series = None, **kwargs):
        return cls(
            *[
                kwargs[f.name] if f.name in kwargs else row[f.name] if row is not None else None
                for f in fields(cls)
            ]
        )

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

    # def __str__(self) -> str:
    #     if self.start == self.end:
    #         return f"{self.exon_id}:e.{self.start}"
    #     else:
    #         return f"{self.exon_id}:e.{self.start}-{self.end}"

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

    # def __str__(self) -> str:
    #     if self.start == self.end:
    #         return f"{self.protein_id}:p.{self.start}"
    #     else:
    #         return f"{self.protein_id}:p.{self.start}-{self.end}"

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

    # def __str__(self) -> str:
    #     if self.start == self.end:
    #         return f"{self.transcript_id}:r.{self.start}"
    #     else:
    #         return f"{self.transcript_id}:r.{self.start}-{self.end}"

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
class CachedEnsemblRelease(type):
    _instances: Dict[Any, EnsemblRelease] = {}
    _current: Optional[EnsemblRelease] = None

    def __call__(cls, *args, **kwargs):
        # convert the given arguments to a hashable tuple
        key = tuple(list(args) + [kwargs[k] for k in sorted(kwargs)])
        # check if an instance created with the same arguments is already cached
        if key in cls._instances:
            instance = cls._instances[key]
        else:
            # create a new instance and cache it
            instance = super(CachedEnsemblRelease, cls).__call__(*args, **kwargs)
            cls._instances[key] = instance

        # track which instance was called last
        cls._current = instance

        return cls._instances[key]


class EnsemblRelease(metaclass=CachedEnsemblRelease):
    def __init__(
        self,
        species: str = DEFAULT_SPECIES,
        reference: str = DEFAULT_REFERENCE,
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
        self.species = species
        self.reference = reference
        self.release = release
        self.cache_dir = cache_dir

        self.cache = Cache(species, reference, release, cache_dir=cache_dir)
        self.ensembl = self.cache.load_gtf()
        self.cdna = self.cache.load_cdna_fasta()
        self.dna = self.cache.load_dna_fasta()
        self.pep = self.cache.load_pep_fasta()
        self.ncrna = self.cache.load_ncrna_fasta()

        self.canonical_transcript = _parse_txt_to_list(canonical_transcript, "canonical transcript")
        self.contig_alias = _parse_tsv_to_dict(contig_alias, "contig aliases")
        self.exon_alias = _parse_tsv_to_dict(exon_alias, "exon aliases")
        self.gene_alias = _parse_tsv_to_dict(gene_alias, "gene aliases")
        self.protein_alias = _parse_tsv_to_dict(protein_alias, "protein aliases")
        self.transcript_alias = _parse_tsv_to_dict(transcript_alias, "transcript aliases")

    # ---------------------------------------------------------------------------------------------
    # <feature_symbol>()
    # ---------------------------------------------------------------------------------------------
    def all_contig_ids(self) -> List[str]:
        return self._uniquify_series(self.ensembl["contig_id"])

    def all_exon_ids(self) -> List[str]:
        return self._uniquify_series(self.ensembl["exon_id"])

    def all_gene_ids(self) -> List[str]:
        return self._uniquify_series(self.ensembl["gene_id"])

    def all_gene_names(self) -> List[str]:
        return self._uniquify_series(self.ensembl["gene_name"])

    def all_protein_ids(self) -> List[str]:
        return self._uniquify_series(self.ensembl["protein_id"])

    def all_transcript_ids(self) -> List[str]:
        return self._uniquify_series(self.ensembl["transcript_id"])

    def all_transcript_names(self) -> List[str]:
        return self._uniquify_series(self.ensembl["transcript_name"])

    # ---------------------------------------------------------------------------------------------
    # get_<feature_symbol>(feature)
    # ---------------------------------------------------------------------------------------------
    def contig_ids(self, feature: str, feature_type: str = "") -> List[str]:
        """Given a feature symbol, return the corresponding contig ID(s)."""
        result = []

        for feature, feature_type in self.normalize_feature(feature, feature_type):
            if is_ensembl_id(feature):
                if feature_type == CONTIG:
                    result.extend(self._contig_ids_of_contig_id(feature))
                elif feature_type == EXON:
                    result.extend(self._contig_ids_of_exon_id(feature))
                elif feature_type == GENE:
                    result.extend(self._contig_ids_of_gene_id(feature))
                elif feature_type == PROTEIN:
                    result.extend(self._contig_ids_of_protein_id(feature))
                elif feature_type in (CDNA, TRANSCRIPT):
                    result.extend(self._contig_ids_of_transcript_id(feature))
                else:
                    raise ValueError(f"Unable to get contig IDs for {feature} ({feature_type})")
            else:
                if feature_type == CONTIG:
                    result.extend(self._contig_ids_of_contig_id(feature))
                elif feature_type == GENE:
                    result.extend(self._contig_ids_of_gene_name(feature))
                elif feature_type in (CDNA, TRANSCRIPT):
                    result.extend(self._contig_ids_of_transcript_name(feature))
                else:
                    raise ValueError(f"Unable to get contig IDs for {feature} ({feature_type})")

        return result

    def exon_ids(self, feature: str, feature_type: str = "") -> List[str]:
        """Given a feature symbol, return the corresponding exon ID(s)."""
        result = []

        for feature, feature_type in self.normalize_feature(feature, feature_type):
            if is_ensembl_id(feature):
                if feature_type == CONTIG:
                    result.extend(self._exon_ids_of_contig_id(feature))
                elif feature_type == EXON:
                    result.extend(self._exon_ids_of_exon_id(feature))
                elif feature_type == GENE:
                    result.extend(self._exon_ids_of_gene_id(feature))
                elif feature_type == PROTEIN:
                    result.extend(self._exon_ids_of_protein_id(feature))
                elif feature_type in (CDNA, TRANSCRIPT):
                    result.extend(self._exon_ids_of_transcript_id(feature))
                else:
                    raise ValueError(f"Unable to get exon IDs for {feature} ({feature_type})")
            else:
                if feature_type == CONTIG:
                    result.extend(self._exon_ids_of_contig_id(feature))
                elif feature_type == GENE:
                    result.extend(self._exon_ids_of_gene_name(feature))
                elif feature_type in (CDNA, TRANSCRIPT):
                    result.extend(self._exon_ids_of_transcript_name(feature))
                else:
                    raise ValueError(f"Unable to get exon IDs for {feature} ({feature_type})")

        return result

    def gene_ids(self, feature: str, feature_type: str = "") -> List[str]:
        """Given a feature symbol, return the corresponding gene ID(s)."""
        result = []

        for feature, feature_type in self.normalize_feature(feature, feature_type):
            if is_ensembl_id(feature):
                if feature_type == CONTIG:
                    result.extend(self._gene_ids_of_contig_id(feature))
                elif feature_type == EXON:
                    result.extend(self._gene_ids_of_exon_id(feature))
                elif feature_type == GENE:
                    result.extend(self._gene_ids_of_gene_id(feature))
                elif feature_type == PROTEIN:
                    result.extend(self._gene_ids_of_protein_id(feature))
                elif feature_type in (CDNA, TRANSCRIPT):
                    result.extend(self._gene_ids_of_transcript_id(feature))
                else:
                    raise ValueError(f"Unable to get gene IDs for {feature} ({feature_type})")
            else:
                if feature_type == CONTIG:
                    result.extend(self._gene_ids_of_contig_id(feature))
                elif feature_type == GENE:
                    result.extend(self._gene_ids_of_gene_name(feature))
                elif feature_type in (CDNA, TRANSCRIPT):
                    result.extend(self._gene_ids_of_transcript_name(feature))
                else:
                    raise ValueError(f"Unable to get gene IDs for {feature} ({feature_type})")

        return result

    def gene_names(self, feature: str, feature_type: str = "") -> List[str]:
        """Given a feature symbol, return the corresponding gene names(s)."""
        result = []

        for feature, feature_type in self.normalize_feature(feature, feature_type):
            if is_ensembl_id(feature):
                if feature_type == CONTIG:
                    result.extend(self._gene_names_of_contig_id(feature))
                elif feature_type == EXON:
                    result.extend(self._gene_names_of_exon_id(feature))
                elif feature_type == GENE:
                    result.extend(self._gene_names_of_gene_id(feature))
                elif feature_type == PROTEIN:
                    result.extend(self._gene_names_of_protein_id(feature))
                elif feature_type in (CDNA, TRANSCRIPT):
                    result.extend(self._gene_names_of_transcript_id(feature))
                else:
                    raise ValueError(f"Unable to get gene names for {feature} ({feature_type})")
            else:
                if feature_type == CONTIG:
                    result.extend(self._gene_names_of_contig_id(feature))
                elif feature_type == GENE:
                    result.extend(self._gene_names_of_gene_name(feature))
                elif feature_type in (CDNA, TRANSCRIPT):
                    result.extend(self._gene_names_of_transcript_name(feature))
                else:
                    raise ValueError(f"Unable to get gene names for {feature} ({feature_type})")

        return result

    def protein_ids(self, feature: str, feature_type: str = "") -> List[str]:
        """Given a feature symbol, return the corresponding protein ID(s)."""
        result = []

        for feature, feature_type in self.normalize_feature(feature, feature_type):
            if is_ensembl_id(feature):
                if feature_type == CONTIG:
                    result.extend(self._protein_ids_of_contig_id(feature))
                elif feature_type == EXON:
                    result.extend(self._protein_ids_of_exon_id(feature))
                elif feature_type == GENE:
                    result.extend(self._protein_ids_of_gene_id(feature))
                elif feature_type == PROTEIN:
                    result.extend(self._protein_ids_of_protein_id(feature))
                elif feature_type in (CDNA, TRANSCRIPT):
                    result.extend(self._protein_ids_of_transcript_id(feature))
                else:
                    raise ValueError(f"Unable to get protein IDs for {feature} ({feature_type})")
            else:
                if feature_type == CONTIG:
                    result.extend(self._protein_ids_of_contig_id(feature))
                elif feature_type == GENE:
                    result.extend(self._protein_ids_of_gene_name(feature))
                elif feature_type in (CDNA, TRANSCRIPT):
                    result.extend(self._protein_ids_of_transcript_name(feature))
                else:
                    raise ValueError(f"Unable to get protein IDs for {feature} ({feature_type})")

        return result

    def transcript_ids(self, feature: str, feature_type: str = "") -> List[str]:
        """Given a feature symbol, return the corresponding transcript ID(s)."""
        result = []

        for feature, feature_type in self.normalize_feature(feature, feature_type):
            if is_ensembl_id(feature):
                if feature_type == CONTIG:
                    result.extend(self._transcript_ids_of_contig_id(feature))
                elif feature_type == EXON:
                    result.extend(self._transcript_ids_of_exon_id(feature))
                elif feature_type == GENE:
                    result.extend(self._transcript_ids_of_gene_id(feature))
                elif feature_type == PROTEIN:
                    result.extend(self._transcript_ids_of_protein_id(feature))
                elif feature_type in (CDNA, TRANSCRIPT):
                    result.extend(self._transcript_ids_of_transcript_id(feature))
                else:
                    raise ValueError(f"Unable to get transcript IDs for {feature} ({feature_type})")
            else:
                if feature_type == CONTIG:
                    result.extend(self._transcript_ids_of_contig_id(feature))
                elif feature_type == GENE:
                    result.extend(self._transcript_ids_of_gene_name(feature))
                elif feature_type in (CDNA, TRANSCRIPT):
                    result.extend(self._transcript_ids_of_transcript_name(feature))
                else:
                    raise ValueError(f"Unable to get transcript IDs for {feature} ({feature_type})")

        return result

    def transcript_names(self, feature: str, feature_type: str = "") -> List[str]:
        """Given a feature symbol, return the corresponding transcript names(s)."""
        result = []

        for feature, feature_type in self.normalize_feature(feature, feature_type):
            if is_ensembl_id(feature):
                if feature_type == CONTIG:
                    result.extend(self._transcript_names_of_contig_id(feature))
                elif feature_type == EXON:
                    result.extend(self._transcript_names_of_exon_id(feature))
                elif feature_type == GENE:
                    result.extend(self._transcript_names_of_gene_id(feature))
                elif feature_type == PROTEIN:
                    result.extend(self._transcript_names_of_protein_id(feature))
                elif feature_type in (CDNA, TRANSCRIPT):
                    result.extend(self._transcript_names_of_transcript_id(feature))
                else:
                    raise ValueError(
                        f"Unable to get transcript names for {feature} ({feature_type})"
                    )
            else:
                if feature_type == CONTIG:
                    result.extend(self._transcript_names_of_contig_id(feature))
                elif feature_type == GENE:
                    result.extend(self._transcript_names_of_gene_name(feature))
                elif feature_type in (CDNA, TRANSCRIPT):
                    result.extend(self._transcript_names_of_transcript_name(feature))
                else:
                    raise ValueError(
                        f"Unable to get transcript names for {feature} ({feature_type})"
                    )

        return result

    def _contig_ids_of_contig_id(self, feature: Union[List[str], str]) -> List[str]:
        return self._query(feature, "contig_id", "contig_id")

    def _contig_ids_of_exon_id(self, feature: Union[List[str], str]) -> List[str]:
        return self._query(feature, "exon_id", "contig_id")

    def _contig_ids_of_gene_id(self, feature: Union[List[str], str]) -> List[str]:
        return self._query(feature, "gene_id", "contig_id")

    def _contig_ids_of_gene_name(self, feature: Union[List[str], str]) -> List[str]:
        return self._query(feature, "gene_name", "contig_id")

    def _contig_ids_of_protein_id(self, feature: Union[List[str], str]) -> List[str]:
        return self._query(feature, "protein_id", "contig_id")

    def _contig_ids_of_transcript_id(self, feature: Union[List[str], str]) -> List[str]:
        return self._query(feature, "transcript_id", "contig_id")

    def _contig_ids_of_transcript_name(self, feature: Union[List[str], str]) -> List[str]:
        return self._query(feature, "transcript_name", "contig_id")

    def _exon_ids_of_contig_id(self, feature: Union[List[str], str]) -> List[str]:
        return self._query(feature, "contig_id", "exon_id")

    def _exon_ids_of_exon_id(self, feature: Union[List[str], str]) -> List[str]:
        return self._query(feature, "exon_id", "exon_id")

    def _exon_ids_of_gene_id(self, feature: Union[List[str], str]) -> List[str]:
        return self._query(feature, "gene_id", "exon_id")

    def _exon_ids_of_gene_name(self, feature: Union[List[str], str]) -> List[str]:
        return self._query(feature, "gene_name", "exon_id")

    def _exon_ids_of_protein_id(self, feature: Union[List[str], str]) -> List[str]:
        return self._query(self._transcript_ids_of_protein_id(feature), "transcript_id", "exon_id")

    def _exon_ids_of_transcript_id(self, feature: Union[List[str], str]) -> List[str]:
        return self._query(feature, "transcript_id", "exon_id")

    def _exon_ids_of_transcript_name(self, feature: Union[List[str], str]) -> List[str]:
        return self._query(feature, "transcript_name", "exon_id")

    def _gene_ids_of_contig_id(self, feature: Union[List[str], str]) -> List[str]:
        return self._query(feature, "contig_id", "gene_id")

    def _gene_ids_of_exon_id(self, feature: Union[List[str], str]) -> List[str]:
        return self._query(feature, "exon_id", "gene_id")

    def _gene_ids_of_gene_id(self, feature: Union[List[str], str]) -> List[str]:
        return self._query(feature, "gene_id", "gene_id")

    def _gene_ids_of_gene_name(self, feature: Union[List[str], str]) -> List[str]:
        return self._query(feature, "gene_name", "gene_id")

    def _gene_ids_of_protein_id(self, feature: Union[List[str], str]) -> List[str]:
        return self._query(feature, "protein_id", "gene_id")

    def _gene_ids_of_transcript_id(self, feature: Union[List[str], str]) -> List[str]:
        return self._query(feature, "transcript_id", "gene_id")

    def _gene_ids_of_transcript_name(self, feature: Union[List[str], str]) -> List[str]:
        return self._query(feature, "transcript_name", "gene_id")

    def _gene_names_of_contig_id(self, feature: Union[List[str], str]) -> List[str]:
        return self._query(feature, "contig_id", "gene_name")

    def _gene_names_of_exon_id(self, feature: Union[List[str], str]) -> List[str]:
        return self._query(feature, "exon_id", "gene_name")

    def _gene_names_of_gene_id(self, feature: Union[List[str], str]) -> List[str]:
        return self._query(feature, "gene_id", "gene_name")

    def _gene_names_of_gene_name(self, feature: Union[List[str], str]) -> List[str]:
        return self._query(feature, "gene_name", "gene_name")

    def _gene_names_of_protein_id(self, feature: Union[List[str], str]) -> List[str]:
        return self._query(feature, "protein_id", "gene_name")

    def _gene_names_of_transcript_id(self, feature: Union[List[str], str]) -> List[str]:
        return self._query(feature, "transcript_id", "gene_name")

    def _gene_names_of_transcript_name(self, feature: Union[List[str], str]) -> List[str]:
        return self._query(feature, "transcript_name", "gene_name")

    def _protein_ids_of_contig_id(self, feature: Union[List[str], str]) -> List[str]:
        return self._query(feature, "contig_id", "protein_id")

    def _protein_ids_of_exon_id(self, feature: Union[List[str], str]) -> List[str]:
        return self._query(self._transcript_ids_of_exon_id(feature), "transcript_id", "protein_id")

    def _protein_ids_of_gene_id(self, feature: Union[List[str], str]) -> List[str]:
        return self._query(feature, "gene_id", "protein_id")

    def _protein_ids_of_gene_name(self, feature: Union[List[str], str]) -> List[str]:
        return self._query(feature, "gene_name", "protein_id")

    def _protein_ids_of_protein_id(self, feature: Union[List[str], str]) -> List[str]:
        return self._query(feature, "protein_id", "protein_id")

    def _protein_ids_of_transcript_id(self, feature: Union[List[str], str]) -> List[str]:
        return self._query(feature, "transcript_id", "protein_id")

    def _protein_ids_of_transcript_name(self, feature: Union[List[str], str]) -> List[str]:
        return self._query(feature, "transcript_name", "protein_id")

    def _transcript_ids_of_contig_id(self, feature: Union[List[str], str]) -> List[str]:
        return self._query(feature, "contig_id", "transcript_id")

    def _transcript_ids_of_exon_id(self, feature: Union[List[str], str]) -> List[str]:
        return self._query(feature, "exon_id", "transcript_id")

    def _transcript_ids_of_gene_id(self, feature: Union[List[str], str]) -> List[str]:
        return self._query(feature, "gene_id", "transcript_id")

    def _transcript_ids_of_gene_name(self, feature: Union[List[str], str]) -> List[str]:
        return self._query(feature, "gene_name", "transcript_id")

    def _transcript_ids_of_protein_id(self, feature: Union[List[str], str]) -> List[str]:
        return self._query(feature, "protein_id", "transcript_id")

    def _transcript_ids_of_transcript_id(self, feature: Union[List[str], str]) -> List[str]:
        return self._query(feature, "transcript_id", "transcript_id")

    def _transcript_ids_of_transcript_name(self, feature: Union[List[str], str]) -> List[str]:
        return self._query(feature, "transcript_name", "transcript_id")

    def _transcript_names_of_contig_id(self, feature: Union[List[str], str]) -> List[str]:
        return self._query(feature, "contig_id", "transcript_name")

    def _transcript_names_of_exon_id(self, feature: Union[List[str], str]) -> List[str]:
        return self._query(feature, "exon_id", "transcript_name")

    def _transcript_names_of_gene_id(self, feature: Union[List[str], str]) -> List[str]:
        return self._query(feature, "gene_id", "transcript_name")

    def _transcript_names_of_gene_name(self, feature: Union[List[str], str]) -> List[str]:
        return self._query(feature, "gene_name", "transcript_name")

    def _transcript_names_of_protein_id(self, feature: Union[List[str], str]) -> List[str]:
        return self._query(feature, "protein_id", "transcript_name")

    def _transcript_names_of_transcript_id(self, feature: Union[List[str], str]) -> List[str]:
        return self._query(feature, "transcript_id", "transcript_name")

    def _transcript_names_of_transcript_name(self, feature: Union[List[str], str]) -> List[str]:
        return self._query(feature, "transcript_name", "transcript_name")

    def _query(self, feature: Union[List[str], str], col: str, key: str) -> List[str]:
        """Generic function for querying the data cache."""
        feature = [feature] if isinstance(feature, str) else feature
        return self._uniquify_series(self.ensembl.loc[self.ensembl[col].isin(feature)][key])

    @staticmethod
    def _uniquify_series(values: pandas.Series) -> List:
        """Convert a `pandas.Series` object into a sorted list of unique values."""
        return sorted((i for i in values.unique().tolist() if i))

    # ---------------------------------------------------------------------------------------------
    # normalize_feature(feature, feature_type)
    # ---------------------------------------------------------------------------------------------
    def normalize_feature(self, feature: str, feature_type: str = "") -> List[Tuple[str, str]]:
        """Normalize a feature to the representation used by Ensembl."""
        normalized = []

        # first, check if the feature is found in the database
        if feature_type == CONTIG or not feature_type:
            if normalized := self._normalize_contig_id(feature):
                feature_type = CONTIG

        if feature_type == EXON or not feature_type:
            if normalized := self._normalize_exon_id(feature):
                feature_type = EXON

        if feature_type == GENE or not feature_type:
            if normalized := self._normalize_gene_id(feature):
                feature_type = GENE
            elif normalized := self._normalize_gene_name(feature):
                feature_type = GENE

        if feature_type == PROTEIN or not feature_type:
            if normalized := self._normalize_protein_id(feature):
                feature_type = PROTEIN

        if feature_type in (CDNA, TRANSCRIPT) or not feature_type:
            if normalized := self._normalize_transcript_id(feature):
                feature_type = TRANSCRIPT
            elif normalized := self._normalize_transcript_name(feature):
                feature_type = TRANSCRIPT

        return [(i, feature_type) for i in normalized]

    def _normalize_contig_id(self, feature: str) -> List[str]:
        """Normalize a contig ID or it's alias to one or more matching Ensembl contig ID."""
        featurel = [feature] + self.get_contig_alias(feature)
        if normalized := self._contig_ids_of_contig_id(featurel):
            return normalized
        else:
            return []

    def _normalize_exon_id(self, feature: str) -> List[str]:
        """Normalize an exon ID or it's alias to one or more matching Ensembl exon ID."""
        featurel = [feature] + self.get_exon_alias(feature)
        if normalized := self._exon_ids_of_exon_id(featurel):
            return normalized
        else:
            return []

    def _normalize_gene_id(self, feature: str) -> List[str]:
        """Normalize a gene ID or it's alias to one or more matching Ensembl gene ID."""
        featurel = [feature] + self.get_gene_alias(feature)
        if normalized := self._gene_ids_of_gene_id(featurel):
            return normalized
        else:
            return []

    def _normalize_gene_name(self, feature: str) -> List[str]:
        """Normalize a gene ID or it's alias to one or more matching Ensembl gene ID."""
        featurel = [feature] + self.get_gene_alias(feature)
        if normalized := self._gene_names_of_gene_name(featurel):
            return normalized
        else:
            return []

    def _normalize_protein_id(self, feature: str) -> List[str]:
        """Normalize a protein ID or it's alias to one or more matching Ensembl protein ID."""
        featurel = [feature] + self.get_protein_alias(feature)
        if normalized := self._protein_ids_of_protein_id(featurel):
            return normalized
        else:
            return []

    def _normalize_transcript_id(self, feature: str) -> List[str]:
        """Normalize a transcript ID or it's alias to one or more matching Ensembl transcript ID."""
        featurel = [feature] + self.get_transcript_alias(feature)
        if normalized := self._transcript_ids_of_transcript_id(featurel):
            return normalized
        else:
            return []

    def _normalize_transcript_name(self, feature: str) -> List[str]:
        """Normalize a transcript name or it's alias to one or more matching Ensembl transcript ID."""
        featurel = [feature] + self.get_transcript_alias(feature)
        if normalized := self._transcript_names_of_transcript_name(featurel):
            return normalized
        else:
            return []

    # ---------------------------------------------------------------------------------------------
    # is_<feature>(feature)
    # ---------------------------------------------------------------------------------------------
    def is_contig(self, feature: str) -> bool:
        """Return True if the given feature is a contig."""
        return any((i[1] == CONTIG for i in self.normalize_feature(feature, CONTIG)))

    def is_exon(self, feature: str) -> bool:
        """Return True if the given feature is an exon."""
        return any((i[1] == EXON for i in self.normalize_feature(feature, EXON)))

    def is_gene(self, feature: str) -> bool:
        """Return True if the given feature is a gene."""
        return any((i[1] == GENE for i in self.normalize_feature(feature, GENE)))

    def is_protein(self, feature: str) -> bool:
        """Return True if the given feature is a protein."""
        return any((i[1] == PROTEIN for i in self.normalize_feature(feature, PROTEIN)))

    def is_transcript(self, feature: str) -> bool:
        """Return True if the given feature is a transcript."""
        return any(
            (i[1] in (CDNA, TRANSCRIPT) for i in self.normalize_feature(feature, TRANSCRIPT))
        )

    # ---------------------------------------------------------------------------------------------
    # sequence
    # ---------------------------------------------------------------------------------------------
    def cdna_sequence(
        self,
        transcript: str,
        start: Optional[int] = None,
        end: Optional[int] = None,
        strand: str = "+",
    ) -> str:
        """Return the nucleotide sequence at the given cDNA coordinates."""
        return self._get_sequence(self.cdna, transcript, start, end, strand)

    def dna_sequence(
        self, contig: str, start: Optional[int] = None, end: Optional[int] = None, strand: str = "+"
    ) -> str:
        """Return the nucleotide sequence at the given contig coordinates."""
        return self._get_sequence(self.dna, contig, start, end, strand)

    def peptide_sequence(
        self,
        protein: str,
        start: Optional[int] = None,
        end: Optional[int] = None,
        strand: str = "+",
    ) -> str:
        """Return the amino acid sequence at the given peptide coordinates."""
        return self._get_sequence(self.pep, protein, start, end, strand)

    def ncrna_sequence(
        self,
        transcript: str,
        start: Optional[int] = None,
        end: Optional[int] = None,
        strand: str = "+",
    ) -> str:
        """Return the nucleotide sequence at the given ncRNA coordinates."""
        return self._get_sequence(self.dna, transcript, start, end, strand)

    def _get_sequence(
        self,
        fasta: Fasta,
        ref: str,
        start: Optional[int] = None,
        end: Optional[int] = None,
        strand: str = "+",
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
    # aliases
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
    # <feature>_to_<feature>
    # ---------------------------------------------------------------------------------------------
    def cdna_to_cdna(
        self, feature: str, start: int, end: Optional[int] = None, strand: Optional[str] = None
    ) -> List[CdnaPosition]:
        """Map a cDNA position to a cDNA position."""
        end = end if end is not None else start
        strandl = [strand] if strand is not None else ["+", "-"]
        transcript_ids = self.transcript_ids(feature)
        result = self._cdna_to_cdna(transcript_ids, start, end, strandl)

        return uniquify(result)

    def cdna_to_dna(
        self, feature: str, start: int, end: Optional[int] = None, strand: Optional[str] = None
    ) -> List[DnaPosition]:
        """Map a cDNA position to a DNA position."""
        end = end if end is not None else start
        strandl = [strand] if strand is not None else ["+", "-"]
        transcript_ids = self.transcript_ids(feature)
        result = self._cdna_to_dna(transcript_ids, start, end, strandl)

        return uniquify(result)

    def cdna_to_exon(
        self, feature: str, start: int, end: Optional[int] = None, strand: Optional[str] = None
    ) -> List[ExonPosition]:
        """Map a cDNA position to an exon position."""
        end = end if end is not None else start
        strandl = [strand] if strand is not None else ["+", "-"]
        transcript_ids = self.transcript_ids(feature)
        result = self._cdna_to_exon(transcript_ids, start, end, strandl)

        return uniquify(result)

    def cdna_to_protein(
        self, feature: str, start: int, end: Optional[int] = None, strand: Optional[str] = None
    ) -> List[ProteinPosition]:
        """Map a cDNA position to a protein position."""
        end = end if end is not None else start
        strandl = [strand] if strand is not None else ["+", "-"]
        transcript_ids = self.transcript_ids(feature)
        result = self._cdna_to_protein(transcript_ids, start, end, strandl)

        return uniquify(result)

    def cdna_to_rna(
        self, feature: str, start: int, end: Optional[int] = None, strand: Optional[str] = None
    ) -> List[RnaPosition]:
        """Map a cDNA position to a RNA position."""
        end = end if end is not None else start
        strandl = [strand] if strand is not None else ["+", "-"]
        transcript_ids = self.transcript_ids(feature)
        result = self._cdna_to_rna(transcript_ids, start, end, strandl)

        return uniquify(result)

    def dna_to_cdna(
        self, feature: str, start: int, end: Optional[int] = None, strand: Optional[str] = None
    ) -> List[CdnaPosition]:
        """Map a DNA position to a cDNA position."""
        end = end if end is not None else start
        strandl = [strand] if strand is not None else ["+", "-"]
        contig_ids = self.contig_ids(feature)
        result = self._dna_to_cdna(contig_ids, start, end, strandl)

        return uniquify(result)

    def dna_to_dna(
        self, feature: str, start: int, end: Optional[int] = None, strand: Optional[str] = None
    ) -> List[DnaPosition]:
        """Map a DNA position to a DNA position."""
        end = end if end is not None else start
        strandl = [strand] if strand is not None else ["+", "-"]
        contig_ids = self.contig_ids(feature)
        result = self._dna_to_dna(contig_ids, start, end, strandl)

        return uniquify(result)

    def dna_to_exon(
        self, feature: str, start: int, end: Optional[int] = None, strand: Optional[str] = None
    ) -> List[ExonPosition]:
        """Map a DNA position to an exon position."""
        end = end if end is not None else start
        strandl = [strand] if strand is not None else ["+", "-"]
        contig_ids = self.contig_ids(feature)
        result = self._dna_to_exon(contig_ids, start, end, strandl)

        return uniquify(result)

    def dna_to_protein(
        self, feature: str, start: int, end: Optional[int] = None, strand: Optional[str] = None
    ) -> List[ProteinPosition]:
        """Map a DNA position to a protein position."""
        end = end if end is not None else start
        strandl = [strand] if strand is not None else ["+", "-"]
        contig_ids = self.contig_ids(feature)
        result = self._dna_to_protein(contig_ids, start, end, strandl)

        return uniquify(result)

    def dna_to_rna(
        self, feature: str, start: int, end: Optional[int] = None, strand: Optional[str] = None
    ) -> List[RnaPosition]:
        """Map a DNA position to a RNA position."""
        end = end if end is not None else start
        strandl = [strand] if strand is not None else ["+", "-"]
        contig_ids = self.contig_ids(feature)
        result = self._dna_to_rna(contig_ids, start, end, strandl)

        return uniquify(result)

    def exon_to_cdna(
        self,
        feature: str,
        start: Optional[int] = None,
        end: Optional[int] = None,
        strand: Optional[str] = None,
    ) -> List[CdnaPosition]:
        """Map an exon position to a cDNA position."""
        return self._to_exon(self._exon_to_cdna, feature, start, end, strand)

    def exon_to_dna(
        self,
        feature: str,
        start: Optional[int] = None,
        end: Optional[int] = None,
        strand: Optional[str] = None,
    ) -> List[DnaPosition]:
        """Map an exon position to a DNA position."""
        return self._to_exon(self._exon_to_dna, feature, start, end, strand)

    def exon_to_exon(
        self,
        feature: str,
        start: Optional[int] = None,
        end: Optional[int] = None,
        strand: Optional[str] = None,
    ) -> List[ExonPosition]:
        """Map an exon position to an exon position."""
        return self._to_exon(self._exon_to_exon, feature, start, end, strand)

    def exon_to_protein(
        self,
        feature: str,
        start: Optional[int] = None,
        end: Optional[int] = None,
        strand: Optional[str] = None,
    ) -> List[ProteinPosition]:
        """Map an exon position to a protein position."""
        return self._to_exon(self._exon_to_protein, feature, start, end, strand)

    def exon_to_rna(
        self,
        feature: str,
        start: Optional[int] = None,
        end: Optional[int] = None,
        strand: Optional[str] = None,
    ) -> List[RnaPosition]:
        """Map an exon position to a RNA position."""
        return self._to_exon(self._exon_to_rna, feature, start, end, strand)

    def _to_exon(
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

        mask = (self.ensembl["exon_id"] == exon_id) & (self.ensembl["feature"] == "exon")
        for _, exon in self.ensembl[mask].iterrows():
            result.append((exon.transcript_id, exon.exon_number, exon.exon_number, exon.strand))

        return result

    def protein_to_cdna(
        self, feature: str, start: int, end: Optional[int] = None, strand: Optional[str] = None
    ) -> List[CdnaPosition]:
        """Map a protein position to a cDNA position."""
        end = end if end is not None else start
        strandl = [strand] if strand is not None else ["+", "-"]
        transcript_ids = self.transcript_ids(feature)
        result = self._protein_to_cdna(transcript_ids, start, end, strandl)

        return uniquify(result)

    def protein_to_dna(
        self, feature: str, start: int, end: Optional[int] = None, strand: Optional[str] = None
    ) -> List[DnaPosition]:
        """Map a protein position to a DNA position."""
        end = end if end is not None else start
        strandl = [strand] if strand is not None else ["+", "-"]
        transcript_ids = self.transcript_ids(feature)
        result = self._protein_to_dna(transcript_ids, start, end, strandl)

        return uniquify(result)

    def protein_to_exon(
        self, feature: str, start: int, end: Optional[int] = None, strand: Optional[str] = None
    ) -> List[ExonPosition]:
        """Map a protein position to an exon position."""
        end = end if end is not None else start
        strandl = [strand] if strand is not None else ["+", "-"]
        transcript_ids = self.transcript_ids(feature)
        result = self._protein_to_exon(transcript_ids, start, end, strandl)

        return uniquify(result)

    def protein_to_protein(
        self, feature: str, start: int, end: Optional[int] = None, strand: Optional[str] = None
    ) -> List[ProteinPosition]:
        """Map a protein position to a protein position."""
        end = end if end is not None else start
        strandl = [strand] if strand is not None else ["+", "-"]
        transcript_ids = self.transcript_ids(feature)
        result = self._protein_to_protein(transcript_ids, start, end, strandl)

        return uniquify(result)

    def protein_to_rna(
        self, feature: str, start: int, end: Optional[int] = None, strand: Optional[str] = None
    ) -> List[RnaPosition]:
        """Map a protein position to a RNA position."""
        end = end if end is not None else start
        strandl = [strand] if strand is not None else ["+", "-"]
        transcript_ids = self.transcript_ids(feature)
        result = self._protein_to_rna(transcript_ids, start, end, strandl)

        return uniquify(result)

    def rna_to_cdna(
        self, feature: str, start: int, end: Optional[int] = None, strand: Optional[str] = None
    ) -> List[CdnaPosition]:
        """Map a RNA position to a cDNA position."""
        end = end if end is not None else start
        strandl = [strand] if strand is not None else ["+", "-"]
        transcript_ids = self.transcript_ids(feature)
        result = self._rna_to_cdna(transcript_ids, start, end, strandl)

        return uniquify(result)

    def rna_to_dna(
        self, feature: str, start: int, end: Optional[int] = None, strand: Optional[str] = None
    ) -> List[DnaPosition]:
        """Map a RNA position to a DNA position."""
        end = end if end is not None else start
        strandl = [strand] if strand is not None else ["+", "-"]
        transcript_ids = self.transcript_ids(feature)
        result = self._rna_to_dna(transcript_ids, start, end, strandl)

        return uniquify(result)

    def rna_to_exon(
        self, feature: str, start: int, end: Optional[int] = None, strand: Optional[str] = None
    ) -> List[ExonPosition]:
        """Map a RNA position to an exon position."""
        end = end if end is not None else start
        strandl = [strand] if strand is not None else ["+", "-"]
        transcript_ids = self.transcript_ids(feature)
        result = self._rna_to_exon(transcript_ids, start, end, strandl)

        return uniquify(result)

    def rna_to_protein(
        self, feature: str, start: int, end: Optional[int] = None, strand: Optional[str] = None
    ) -> List[ProteinPosition]:
        """Map a RNA position to a protein position."""
        end = end if end is not None else start
        strandl = [strand] if strand is not None else ["+", "-"]
        transcript_ids = self.transcript_ids(feature)
        result = self._rna_to_protein(transcript_ids, start, end, strandl)

        return uniquify(result)

    def rna_to_rna(
        self, feature: str, start: int, end: Optional[int] = None, strand: Optional[str] = None
    ) -> List[RnaPosition]:
        """Map a RNA position to a RNA position."""
        end = end if end is not None else start
        strandl = [strand] if strand is not None else ["+", "-"]
        transcript_ids = self.transcript_ids(feature)
        result = self._rna_to_rna(transcript_ids, start, end, strandl)

        return uniquify(result)

    def _cdna_to_cdna(
        self, transcript_ids: List[str], start: int, end: int, strand: List[str]
    ) -> List[CdnaPosition]:
        def convert(position: int):
            result = []

            mask = (
                (self.ensembl["transcript_id"].isin(transcript_ids))
                & (self.ensembl["cdna_start"] <= position)
                & (self.ensembl["cdna_end"] >= position)
                & (self.ensembl["strand"].isin(strand))
                & (self.ensembl["feature"].isin(["CDS", "stop_codon"]))
            )
            for _, cds in self.ensembl[mask].iterrows():
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
            return merge_positions(result_start, result_end, "transcript_id")

    def _cdna_to_dna(
        self, transcript_ids: List[str], start: int, end: int, strand: List[str]
    ) -> List[DnaPosition]:
        def convert(position: int):
            result = []

            mask = (
                (self.ensembl["transcript_id"].isin(transcript_ids))
                & (self.ensembl["cdna_start"] <= position)
                & (self.ensembl["cdna_end"] >= position)
                & (self.ensembl["strand"].isin(strand))
                & (self.ensembl["feature"].isin(["CDS", "stop_codon"]))
            )
            for _, cds in self.ensembl[mask].iterrows():
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
            return merge_positions(result_start, result_end, "contig_id")

    def _cdna_to_exon(
        self, transcript_ids: List[str], start: int, end: int, strand: List[str]
    ) -> List[ExonPosition]:
        def convert(position: int):
            result = []

            mask_cds = (
                (self.ensembl["transcript_id"].isin(transcript_ids))
                & (self.ensembl["cdna_start"] <= position)
                & (self.ensembl["cdna_end"] >= position)
                & (self.ensembl["strand"].isin(strand))
                & (self.ensembl["feature"].isin(["CDS", "stop_codon"]))
            )
            for _, cds in self.ensembl[mask_cds].iterrows():
                mask_exon = (
                    (self.ensembl["transcript_id"].isin(transcript_ids))
                    & (self.ensembl["exon_number"] == cds.exon_number)
                    & (self.ensembl["feature"] == "exon")
                )
                for _, exon in self.ensembl[mask_exon].iterrows():
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
        for cdna in self._cdna_to_cdna(transcript_ids, start, end, strand):
            pstart = convert(start)
            pend = convert(end)
            protein.append(ProteinPosition.copy_from(cdna, start=pstart, end=pend))

        return protein

    def _cdna_to_rna(
        self, transcript_ids: List[str], start: int, end: int, strand: List[str]
    ) -> List[RnaPosition]:
        def convert(position: int):
            result = []

            mask = (
                (self.ensembl["transcript_id"].isin(transcript_ids))
                & (self.ensembl["cdna_start"] <= position)
                & (self.ensembl["cdna_end"] >= position)
                & (self.ensembl["strand"].isin(strand))
                & (self.ensembl["feature"].isin(["CDS", "stop_codon"]))
            )
            for _, cds in self.ensembl[mask].iterrows():
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
            return merge_positions(result_start, result_end, "transcript_id")

    def _dna_to_cdna(
        self, contig_ids: List[str], start: int, end: int, strand: List[str]
    ) -> List[CdnaPosition]:
        def convert(position: int):
            result = []

            mask = (
                (self.ensembl["contig_id"].isin(contig_ids))
                & (self.ensembl["start"] <= position)
                & (self.ensembl["end"] >= position)
                & (self.ensembl["strand"].isin(strand))
                & (self.ensembl["feature"].isin(["CDS", "stop_codon"]))
            )
            for _, cds in self.ensembl[mask].iterrows():
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
            return merge_positions(result_start, result_end, "transcript_id")

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
                (self.ensembl["contig_id"].isin(contig_ids))
                & (self.ensembl["start"] <= position)
                & (self.ensembl["end"] >= position)
                & (self.ensembl["strand"].isin(strand))
                & (self.ensembl["feature"].isin(["exon"]))
            )
            for _, exon in self.ensembl[mask].iterrows():
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
        for cdna in self._dna_to_cdna(contig_ids, start, end, strand):
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
                (self.ensembl["contig_id"].isin(contig_ids))
                & (self.ensembl["start"] <= position)
                & (self.ensembl["end"] >= position)
                & (self.ensembl["strand"].isin(strand))
                & (self.ensembl["feature"] == "exon")
            )
            for _, exon in self.ensembl[mask].iterrows():
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
            return merge_positions(result_start, result_end, "transcript_id")

    def _exon_to_cdna(
        self, transcript_ids: List[str], start: int, end: int, strand: List[str]
    ) -> List[CdnaPosition]:
        def convert(position: int):
            result = []

            mask = (
                (self.ensembl["transcript_id"].isin(transcript_ids))
                & (self.ensembl["exon_number"] == str(position))  # TODO: exon number should be int
                & (self.ensembl["strand"].isin(strand))
                & (self.ensembl["feature"].isin(["CDS", "stop_codon"]))
            )
            for _, cds in self.ensembl[mask].iterrows():
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
            return merge_positions(result_start, result_end, "transcript_id")

    def _exon_to_dna(
        self, transcript_ids: List[str], start: int, end: int, strand: List[str]
    ) -> List[DnaPosition]:
        def convert(position: int):
            result = []

            mask = (
                (self.ensembl["transcript_id"].isin(transcript_ids))
                # TODO: exon number should be int
                & (self.ensembl["exon_number"] == str(position))
                & (self.ensembl["strand"].isin(strand))
                & (self.ensembl["feature"] == "exon")
            )
            for _, exon in self.ensembl[mask].iterrows():
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
            return merge_positions(result_start, result_end, "contig_id")

    def _exon_to_exon(
        self, transcript_ids: List[str], start: int, end: int, strand: List[str]
    ) -> List[ExonPosition]:
        def convert(position: int):
            result = []

            mask = (
                (self.ensembl["transcript_id"].isin(transcript_ids))
                # TODO: exon number should be int
                & (self.ensembl["exon_number"] == str(position))
                & (self.ensembl["strand"].isin(strand))
                & (self.ensembl["feature"] == "exon")
            )
            for _, exon in self.ensembl[mask].iterrows():
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
            return merge_positions(result_start, result_end, "transcript_id")

    def _exon_to_protein(
        self, transcript_ids: List[str], start: int, end: int, strand: List[str]
    ) -> List[ProteinPosition]:
        result = []
        for cdna in self._exon_to_cdna(transcript_ids, start, end, strand):
            result.extend(cdna.to_protein())

        return result

    def _exon_to_rna(
        self, transcript_ids: List[str], start: int, end: int, strand: List[str]
    ) -> List[RnaPosition]:
        def convert(position: int):
            result = []

            mask = (
                (self.ensembl["transcript_id"].isin(transcript_ids))
                & (self.ensembl["exon_number"] == str(position))  # TODO: exon number should be int
                & (self.ensembl["strand"].isin(strand))
                & (self.ensembl["feature"] == "exon")
            )
            for _, exon in self.ensembl[mask].iterrows():
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
            return merge_positions(result_start, result_end, "transcript_id")

    def _protein_to_cdna(
        self, transcript_ids: List[str], start: int, end: int, strand: List[str]
    ) -> List[CdnaPosition]:
        def convert(position: int):
            return (position - 1) * 3 + 1

        cdna_start = convert(start)
        cdna_end = convert(end) + 2

        return self._cdna_to_cdna(transcript_ids, cdna_start, cdna_end, strand)

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
        self, transcript_ids: List[str], start: int, end: int, strand: List[str]
    ) -> List[CdnaPosition]:
        def convert(position: int):
            result = []

            mask = (
                (self.ensembl["transcript_id"].isin(transcript_ids))
                & (self.ensembl["transcript_start"] <= position)
                & (self.ensembl["transcript_end"] >= position)
                & (self.ensembl["strand"].isin(strand))
                & (self.ensembl["feature"].isin(["CDS", "stop_codon"]))
            )
            for _, cds in self.ensembl[mask].iterrows():
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
            return merge_positions(result_start, result_end, "transcript_id")

    def _rna_to_dna(
        self, transcript_ids: List[str], start: int, end: int, strand: List[str]
    ) -> List[DnaPosition]:
        def convert(position: int):
            result = []

            mask = (
                (self.ensembl["transcript_id"].isin(transcript_ids))
                & (self.ensembl["transcript_start"] <= position)
                & (self.ensembl["transcript_end"] >= position)
                & (self.ensembl["strand"].isin(strand))
                & (self.ensembl["feature"].isin(["exon"]))
            )
            exon_df = self.ensembl[mask]
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
            return merge_positions(result_start, result_end, "contig_id")

    def _rna_to_exon(
        self, transcript_ids: List[str], start: int, end: int, strand: List[str]
    ) -> List[ExonPosition]:
        def convert(position: int):
            result = []

            mask = (
                (self.ensembl["transcript_id"].isin(transcript_ids))
                & (self.ensembl["transcript_start"] <= position)
                & (self.ensembl["transcript_end"] >= position)
                & (self.ensembl["strand"].isin(strand))
                & (self.ensembl["feature"].isin(["exon"]))
            )
            for _, exon in self.ensembl[mask].iterrows():
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

    def _rna_to_protein(
        self, transcript_ids: List[str], start: int, end: int, strand: List[str]
    ) -> List[ProteinPosition]:
        result = []
        for cdna in self._rna_to_cdna(transcript_ids, start, end, strand):
            result.extend(cdna.to_protein())

        return result

    def _rna_to_rna(
        self, transcript_ids: List[str], start: int, end: int, strand: List[str]
    ) -> List[RnaPosition]:
        def convert(position: int):
            result = []

            mask = (
                (self.ensembl["transcript_id"].isin(transcript_ids))
                & (self.ensembl["transcript_start"] <= position)
                & (self.ensembl["transcript_end"] >= position)
                & (self.ensembl["strand"].isin(strand))
                & (self.ensembl["feature"] == "exon")
            )
            for _, exon in self.ensembl[mask].iterrows():
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
            return merge_positions(result_start, result_end, "transcript_id")


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

    start_dict = {getattr(pos, key): pos for pos in start_positions}
    end_dict = {getattr(pos, key): pos for pos in end_positions}
    for key, start_pos in start_dict.items():
        if end_pos := end_dict.get(key):
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
