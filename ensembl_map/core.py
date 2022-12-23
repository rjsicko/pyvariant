from __future__ import annotations

from dataclasses import dataclass, fields
from itertools import product
from math import floor
from typing import Any, Callable, Dict, List, Optional, Tuple, TypeVar, Union

import pandas as pd
from gtfparse import read_gtf
from pyfaidx import Fasta

from .cache import EnsemblCache
from .constants import (
    CONTIG_ID,
    DEFAULT_ENSEMBL_RELEASE,
    DEFAULT_SPECIES,
    DELETION,
    DELINS,
    DUPLICATION,
    EXON_ID,
    FRAMESHIFT,
    FUSION,
    GENE_ID,
    GENE_NAME,
    INSERTION,
    PROTEIN_ID,
    SUBTITUTION,
    TRANSCRIPT_ID,
    TRANSCRIPT_NAME,
)
from .files import read_fasta, tsv_to_dict, txt_to_list
from .tables import AMINO_ACID_TABLE
from .utils import (
    collapse_mutation,
    expand_nt,
    format_hgvs_position,
    is_frameshift,
    reverse_complement,
    reverse_translate,
    split_by_codon,
    strip_version,
)


# -------------------------------------------------------------------------------------------------
# Position classes
# -------------------------------------------------------------------------------------------------
@dataclass(eq=True, frozen=True)
class Position:
    """Base class for position data objects."""

    _data: Core
    contig_id: str
    start: int
    start_offset: int
    end: int
    end_offset: int
    strand: str

    @classmethod
    def copy_from(cls, position: Position, **kwargs):
        return cls(
            **{
                **{
                    k: v for k, v in position.asdict().items() if k in [i.name for i in fields(cls)]
                },
                **kwargs,
            }
        )

    def __getitem__(self, item: Any) -> Any:
        return getattr(self, item)

    def __lt__(self, other: Position) -> bool:
        return str(self) < str(other)

    def __str__(self) -> str:
        raise NotImplementedError()

    @property
    def is_cdna(self) -> bool:
        return False

    @property
    def is_dna(self) -> bool:
        return False

    @property
    def is_exon(self) -> bool:
        return False

    @property
    def is_protein(self) -> bool:
        return False

    @property
    def is_rna(self) -> bool:
        return False

    @property
    def on_negative_strand(self) -> bool:
        return self.strand == "-"

    @property
    def on_positive_strand(self) -> bool:
        return self.strand == "+"

    def asdict(self) -> Dict[str, Any]:
        return {f.name: self[f.name] for f in fields(self)}

    def sequence(self) -> str:
        raise NotImplementedError()


@dataclass(eq=True, frozen=True)
class CdnaPosition(Position):
    """Stores information on cDNA data objects."""

    gene_id: str
    gene_name: str
    transcript_id: str
    transcript_name: str
    protein_id: str

    def __str__(self) -> str:
        if self.start == self.end:
            return f"{self.transcript_id}:c.{self.start}"
        else:
            return f"{self.transcript_id}:c.{self.start}_{self.end}"

    @property
    def is_cdna(self) -> bool:
        return True

    def sequence(self) -> str:
        if self.start_offset or self.end_offset:
            raise ValueError(
                f"Unable to get cDNA sequence of offset position ({self.start_offset}, {self.end_offset})"
            )

        return self._data.cds_sequence(self.transcript_id, self.start, self.end)


@dataclass(eq=True, frozen=True)
class DnaPosition(Position):
    """Stores information on DNA data objects."""

    def __str__(self) -> str:
        if self.start == self.end:
            return f"{self.contig_id}:g.{self.start}"
        else:
            return f"{self.contig_id}:g.{self.start}_{self.end}"

    @property
    def is_dna(self) -> bool:
        return True

    def sequence(self) -> str:
        if self.start_offset or self.end_offset:
            raise ValueError(
                f"Unable to get DNA sequence of offset position ({self.start_offset}, {self.end_offset})"
            )

        return self._data.dna_sequence(self.contig_id, self.start, self.end, self.strand)


@dataclass(eq=True, frozen=True)
class ExonPosition(Position):
    """Stores information on exon data objects."""

    gene_id: str
    gene_name: str
    transcript_id: str
    transcript_name: str
    exon_id: str

    def __str__(self) -> str:
        if self.start == self.end:
            return f"{self.transcript_id}:e.{self.start}"
        else:
            return f"{self.transcript_id}:e.{self.start}_{self.end}"

    @property
    def is_exon(self) -> bool:
        return True

    def sequence(self) -> str:
        raise NotImplementedError()  # TODO


@dataclass(eq=True, frozen=True)
class ProteinPosition(Position):
    """Stores information on protein data objects."""

    gene_id: str
    gene_name: str
    transcript_id: str
    transcript_name: str
    protein_id: str

    def __str__(self) -> str:
        if self.start == self.end:
            return f"{self.protein_id}:p.{self.start}"
        else:
            return f"{self.protein_id}:p.{self.start}_{self.end}"

    @property
    def is_protein(self) -> bool:
        return True

    def sequence(self) -> str:
        if self.start_offset or self.end_offset:
            raise ValueError(
                f"Unable to get peptide sequence of offset position ({self.start_offset}, {self.end_offset})"
            )

        return self._data.peptide_sequence(self.protein_id, self.start, self.end)


@dataclass(eq=True, frozen=True)
class RnaPosition(Position):
    """Stores information on RNA data objects."""

    gene_id: str
    gene_name: str
    transcript_id: str
    transcript_name: str

    def __str__(self) -> str:
        if self.start == self.end:
            return f"{self.transcript_id}:r.{self.start}"
        else:
            return f"{self.transcript_id}:r.{self.start}_{self.end}"

    @property
    def is_rna(self) -> bool:
        return True

    def sequence(self) -> str:
        if self.start_offset or self.end_offset:
            raise ValueError(
                f"Unable to get RNA sequence of offset position ({self.start_offset}, {self.end_offset})"
            )

        return self._data.rna_sequence(self.transcript_id, self.start, self.end)


# -------------------------------------------------------------------------------------------------
# MappablePosition classes
# -------------------------------------------------------------------------------------------------
@dataclass(eq=True, frozen=True)
class MappablePosition:
    """Base class for mappable position data objects."""

    def to_cdna(self) -> List[CdnaMappablePosition]:
        raise NotImplementedError()

    def to_dna(self) -> List[DnaMappablePosition]:
        raise NotImplementedError()

    def to_exon(self) -> List[ExonMappablePosition]:
        raise NotImplementedError()

    def to_protein(self) -> List[ProteinMappablePosition]:
        raise NotImplementedError()

    def to_rna(self) -> List[RnaMappablePosition]:
        raise NotImplementedError()


@dataclass(eq=True, frozen=True)
class CdnaMappablePosition(CdnaPosition, MappablePosition):
    """Stores information on a cDNA position and converts to other position types."""

    def to_cdna(self) -> List[CdnaMappablePosition]:
        return self._data._cdna_to_cdna(
            [self.transcript_id],
            self.start,
            self.start_offset,
            self.end,
            self.end_offset,
            [self.strand],
        )

    def to_dna(self) -> List[DnaMappablePosition]:
        return self._data._cdna_to_dna(
            [self.transcript_id],
            self.start,
            self.start_offset,
            self.end,
            self.end_offset,
            [self.strand],
        )

    def to_exon(self) -> List[ExonMappablePosition]:
        return self._data._cdna_to_exon(
            [self.transcript_id],
            self.start,
            self.start_offset,
            self.end,
            self.end_offset,
            [self.strand],
        )

    def to_rna(self) -> List[RnaMappablePosition]:
        return self._data._cdna_to_rna(
            [self.transcript_id],
            self.start,
            self.start_offset,
            self.end,
            self.end_offset,
            [self.strand],
        )

    def to_protein(self) -> List[ProteinMappablePosition]:
        return self._data._cdna_to_protein(
            [self.transcript_id],
            self.start,
            self.start_offset,
            self.end,
            self.end_offset,
            [self.strand],
        )


@dataclass(eq=True, frozen=True)
class DnaMappablePosition(DnaPosition, MappablePosition):
    """Stores information on a DNA position and converts to other position types."""

    def to_cdna(self) -> List:
        return self._data._dna_to_cdna(
            [self.contig_id],
            self.start,
            self.start_offset,
            self.end,
            self.end_offset,
            [self.strand],
        )

    def to_dna(self) -> List:
        return self._data._dna_to_dna(
            [self.contig_id],
            self.start,
            self.start_offset,
            self.end,
            self.end_offset,
            [self.strand],
        )

    def to_exon(self) -> List:
        return self._data._dna_to_exon(
            [self.contig_id],
            self.start,
            self.start_offset,
            self.end,
            self.end_offset,
            [self.strand],
        )

    def to_protein(self) -> List:
        return self._data._dna_to_protein(
            [self.contig_id],
            self.start,
            self.start_offset,
            self.end,
            self.end_offset,
            [self.strand],
        )

    def to_rna(self) -> List:
        return self._data._dna_to_rna(
            [self.contig_id],
            self.start,
            self.start_offset,
            self.end,
            self.end_offset,
            [self.strand],
        )


@dataclass(eq=True, frozen=True)
class ExonMappablePosition(ExonPosition, MappablePosition):
    """Stores information on an exon position and converts to other position types."""

    def to_cdna(self) -> List:
        return self._data._exon_to_cdna(
            [self.transcript_id],
            self.start,
            self.start_offset,
            self.end,
            self.end_offset,
            [self.strand],
        )

    def to_dna(self) -> List:
        return self._data._exon_to_dna(
            [self.transcript_id],
            self.start,
            self.start_offset,
            self.end,
            self.end_offset,
            [self.strand],
        )

    def to_exon(self) -> List:
        return self._data._exon_to_exon(
            [self.transcript_id],
            self.start,
            self.start_offset,
            self.end,
            self.end_offset,
            [self.strand],
        )

    def to_protein(self) -> List:
        return self._data._exon_to_protein(
            [self.transcript_id],
            self.start,
            self.start_offset,
            self.end,
            self.end_offset,
            [self.strand],
        )

    def to_rna(self) -> List:
        return self._data._exon_to_rna(
            [self.transcript_id],
            self.start,
            self.start_offset,
            self.end,
            self.end_offset,
            [self.strand],
        )


@dataclass(eq=True, frozen=True)
class ProteinMappablePosition(ProteinPosition, MappablePosition):
    """Stores information on a protein position and converts to other position types."""

    def to_cdna(self) -> List:
        return self._data._protein_to_cdna(
            [self.transcript_id],
            self.start,
            self.start_offset,
            self.end,
            self.end_offset,
            [self.strand],
        )

    def to_dna(self) -> List:
        return self._data._protein_to_dna(
            [self.transcript_id],
            self.start,
            self.start_offset,
            self.end,
            self.end_offset,
            [self.strand],
        )

    def to_exon(self) -> List:
        return self._data._protein_to_exon(
            [self.transcript_id],
            self.start,
            self.start_offset,
            self.end,
            self.end_offset,
            [self.strand],
        )

    def to_protein(self) -> List:
        return self._data._protein_to_protein(
            [self.transcript_id],
            self.start,
            self.start_offset,
            self.end,
            self.end_offset,
            [self.strand],
        )

    def to_rna(self) -> List:
        return self._data._protein_to_rna(
            [self.transcript_id],
            self.start,
            self.start_offset,
            self.end,
            self.end_offset,
            [self.strand],
        )


@dataclass(eq=True, frozen=True)
class RnaMappablePosition(RnaPosition, MappablePosition):
    """Stores information on an RNA position and converts to other position types."""

    def to_cdna(self) -> List:
        return self._data._rna_to_cdna(
            [self.transcript_id],
            self.start,
            self.start_offset,
            self.end,
            self.end_offset,
            [self.strand],
        )

    def to_dna(self) -> List:
        return self._data._rna_to_dna(
            [self.transcript_id],
            self.start,
            self.start_offset,
            self.end,
            self.end_offset,
            [self.strand],
        )

    def to_exon(self) -> List:
        return self._data._rna_to_exon(
            [self.transcript_id],
            self.start,
            self.start_offset,
            self.end,
            self.end_offset,
            [self.strand],
        )

    def to_protein(self) -> List:
        return self._data._rna_to_protein(
            [self.transcript_id],
            self.start,
            self.start_offset,
            self.end,
            self.end_offset,
            [self.strand],
        )

    def to_rna(self) -> List:
        return self._data._rna_to_rna(
            [self.transcript_id],
            self.start,
            self.start_offset,
            self.end,
            self.end_offset,
            [self.strand],
        )


# -------------------------------------------------------------------------------------------------
# SmallVariant classes
# -------------------------------------------------------------------------------------------------
@dataclass(eq=True, frozen=True)
class SmallVariant:
    """Base class for small variant data objects."""

    refseq: str
    altseq: str

    def __str__(self) -> str:
        raise NotImplementedError()

    @property
    def is_deletion(self) -> bool:
        return self.type == DELETION

    @property
    def is_delins(self) -> bool:
        return self.type == DELINS

    @property
    def is_duplication(self) -> bool:
        return self.type == DUPLICATION

    @property
    def is_frameshift(self) -> bool:
        return self.type == FRAMESHIFT

    @property
    def is_fusion(self) -> bool:
        return self.type == FUSION

    @property
    def is_insertion(self) -> bool:
        return self.type == INSERTION

    @property
    def is_substitution(self) -> bool:
        return self.type == SUBTITUTION

    @property
    def type(self) -> str:
        raise NotImplementedError()

    def to_cdna(self) -> List:
        raise NotImplementedError()

    def to_dna(self) -> List:
        raise NotImplementedError()

    def to_exon(self) -> List:
        raise NotImplementedError()

    def to_protein(self) -> List:
        raise NotImplementedError()

    def to_rna(self) -> List:
        raise NotImplementedError()


@dataclass(eq=True, frozen=True)
class CdnaSmallVariant(CdnaPosition, SmallVariant):
    """Base class for cDNA variant data objects."""

    @classmethod
    def from_cdna(
        cls, cdna: CdnaMappablePosition, refseq: str, altseq: str
    ) -> List[CdnaSmallVariant]:
        """Convert a cDNA position plus ref/alt sequences into a cDNA variant data object."""
        variant_list = []

        ref_annotated = cdna.sequence()
        for ref, alt in product(expand_nt(refseq), expand_nt(altseq)):
            # Assert that the given ref matches the annotated one
            if ref != ref_annotated:
                raise ValueError(
                    f"Given ref allele '{ref}' does not match annotated ref allele '{ref_annotated}'"
                )

            # Remove bases that are unchanged between the ref and alt alleles
            new_ref, new_alt, new_start, new_end = collapse_mutation(ref, alt)
            start = cdna.start + new_start
            end = cdna.end - new_end

            # Determine the type of variant
            if new_alt == new_ref * 2:
                # If the bases are a copy of the bases immediately 5' it is an duplication
                variant = CdnaDuplication.copy_from(
                    cdna, start=start, end=end, refseq=new_ref, altseq=new_alt
                )
            elif len(new_ref) == 1 and len(new_alt) == 1:
                # If exactly one base changes it is a substitution
                variant = CdnaSubstitution.copy_from(
                    cdna, start=start, end=end, refseq=new_ref, altseq=new_alt
                )
            elif len(new_ref) > 0 and len(new_alt) == 0:
                # If only the reference sequence is removed it is a deletion
                variant = CdnaDeletion.copy_from(
                    cdna, start=start, end=end, refseq=new_ref, altseq=new_alt
                )
            elif len(new_ref) == 0 and len(new_alt) > 0:
                # If only new bases are added it is an insertion
                variant = CdnaInsertion.copy_from(
                    cdna, start=start, end=end, refseq=new_ref, altseq=new_alt
                )
            elif new_alt == new_ref * 2:
                # If the bases are a copy of the bases immediately 5' it is an duplication
                variant = CdnaDuplication.copy_from(
                    cdna, start=start, end=end, refseq=new_ref, altseq=new_alt
                )
            elif len(new_ref) * len(new_alt) > 1:
                # If more than one ref and alt base is changed it is a delins
                variant = CdnaDelins.copy_from(
                    cdna, start=start, end=end, refseq=new_ref, altseq=new_alt
                )
            else:
                raise NotImplementedError(f"Unknown ref/alt '{new_ref}/{new_alt}' for {cdna}")

            variant_list.append(variant)

        return variant_list

    @classmethod
    def from_protein(
        cls, cdna: CdnaMappablePosition, refseq: str, altseq: str
    ) -> List[CdnaSmallVariant]:
        """Convert a cDNA position plus ref/alt amino acids into a cDNA variant data object."""
        variant_list = []

        ref_annotated = cdna.sequence()
        for ref, alt in product(reverse_translate(refseq), reverse_translate(altseq)):
            # Assert that the given ref matches the annotated one
            if ref != ref_annotated:
                continue

            # Remove bases that are unchanged between the ref and alt alleles
            new_ref, new_alt, new_start, new_end = collapse_mutation(ref, alt)
            start = cdna.start + new_start
            end = cdna.end - new_end

            # Determine the type of variant
            if len(new_ref) == 1 and len(new_alt) == 1:
                # If exactly one base changes it is a substitution
                variant = CdnaSubstitution.copy_from(
                    cdna, start=start, end=end, refseq=new_ref, altseq=new_alt
                )
            elif len(new_ref) > 0 and len(new_alt) == 0:
                # If only the reference sequence is removed it is a deletion
                variant = CdnaDeletion.copy_from(
                    cdna, start=start, end=end, refseq=new_ref, altseq=new_alt
                )
            elif len(new_ref) == 0 and len(new_alt) > 0:
                # If only new bases are added it is an insertion
                variant = CdnaInsertion.copy_from(
                    cdna, start=start, end=end, refseq=new_ref, altseq=new_alt
                )
            elif new_alt == new_ref * 2:
                # If the bases are a copy of the bases immediately 5' it is an duplication
                variant = CdnaDuplication.copy_from(
                    cdna, start=start, end=end, refseq=new_ref, altseq=new_alt
                )
            elif len(new_ref) * len(new_alt) > 1:
                # If more than one ref and alt base is changed it is a delins
                variant = CdnaDelins.copy_from(
                    cdna, start=start, end=end, refseq=new_ref, altseq=new_alt
                )
            else:
                raise NotImplementedError(f"Unknown ref/alt '{new_ref}/{new_alt}' for {cdna}")

            variant_list.append(variant)

        return variant_list

    def to_cdna(self) -> List[CdnaSmallVariant]:
        return self._data._cdna_to_cdna_variant(
            [self.transcript_id],
            self.start,
            self.start_offset,
            self.end,
            self.end_offset,
            [self.strand],
            self.refseq,
            self.altseq,
        )

    def to_dna(self) -> List[DnaSmallVariant]:
        return self._data._cdna_to_dna_variant(
            [self.transcript_id],
            self.start,
            self.start_offset,
            self.end,
            self.end_offset,
            [self.strand],
            self.refseq,
            self.altseq,
        )

    def to_exon(self) -> List[ExonSmallVariant]:
        raise NotImplementedError()

    def to_protein(self) -> List[ProteinSmallVariant]:
        return self._data._cdna_to_protein_variant(
            [self.transcript_id],
            self.start,
            self.start_offset,
            self.end,
            self.end_offset,
            [self.strand],
            self.refseq,
            self.altseq,
        )

    def to_rna(self) -> List[RnaSmallVariant]:
        return self._data._cdna_to_rna_variant(
            [self.transcript_id],
            self.start,
            self.start_offset,
            self.end,
            self.end_offset,
            [self.strand],
            self.refseq,
            self.altseq,
        )


@dataclass(eq=True, frozen=True)
class DnaSmallVariant(DnaPosition, SmallVariant):
    """Base class for DNA variant data objects."""

    @classmethod
    def from_dna(cls, dna: DnaMappablePosition, refseq: str, altseq: str) -> List[DnaSmallVariant]:
        """Convert a DNA position plus ref/alt sequences into a DNA variant data object."""
        variant_list = []

        ref_annotated = dna.sequence()
        for ref, alt in product(expand_nt(refseq), expand_nt(altseq)):
            # Assert that the given ref matches the annotated one
            if ref != ref_annotated:
                raise ValueError(
                    f"Given ref allele '{ref}' does not match annotated ref allele '{ref_annotated}'"
                )

            # Remove bases that are unchanged between the ref and alt alleles
            new_ref, new_alt, new_start, new_end = collapse_mutation(ref, alt)
            start = dna.start + new_start
            end = dna.end - new_end

            # Determine the type of variant
            if len(new_ref) == 1 and len(new_alt) == 1:
                # If exactly one base changes it is a substitution
                variant = DnaSubstitution.copy_from(
                    dna, start=start, end=end, refseq=new_ref, altseq=new_alt
                )
            elif len(new_ref) > 0 and len(new_alt) == 0:
                # If only the reference sequence is removed it is a deletion
                variant = DnaDeletion.copy_from(
                    dna, start=start, end=end, refseq=new_ref, altseq=new_alt
                )
            elif len(new_ref) == 0 and len(new_alt) > 0:
                # If only new bases are added it is an insertion
                variant = DnaInsertion.copy_from(
                    dna, start=start, end=end, refseq=new_ref, altseq=new_alt
                )
            elif new_alt == new_ref * 2:
                # If the bases are a copy of the bases immediately 5' it is an duplication
                variant = DnaDuplication.copy_from(
                    dna, start=start, end=end, refseq=new_ref, altseq=new_alt
                )
            elif len(new_ref) * len(new_alt) > 1:
                # If more than one ref and alt base is changed it is a delins
                variant = DnaDelins.copy_from(
                    dna, start=start, end=end, refseq=new_ref, altseq=new_alt
                )
            else:
                raise NotImplementedError(f"Unknown ref/alt '{new_ref}/{new_alt}' for {dna}")

            variant_list.append(variant)

        return variant_list

    def to_cdna(self) -> List[CdnaSmallVariant]:
        return self._data._dna_to_cdna_variant(
            [self.contig_id],
            self.start,
            self.start_offset,
            self.end,
            self.end_offset,
            [self.strand],
            self.refseq,
            self.altseq,
        )

    def to_dna(self) -> List[DnaSmallVariant]:
        return self._data._dna_to_dna_variant(
            [self.contig_id],
            self.start,
            self.start_offset,
            self.end,
            self.end_offset,
            [self.strand],
            self.refseq,
            self.altseq,
        )

    def to_exon(self) -> List[ExonSmallVariant]:
        raise NotImplementedError()

    def to_protein(self) -> List[ProteinSmallVariant]:
        return self._data._dna_to_protein_variant(
            [self.contig_id],
            self.start,
            self.start_offset,
            self.end,
            self.end_offset,
            [self.strand],
            self.refseq,
            self.altseq,
        )

    def to_rna(self) -> List[RnaSmallVariant]:
        return self._data._dna_to_rna_variant(
            [self.contig_id],
            self.start,
            self.start_offset,
            self.end,
            self.end_offset,
            [self.strand],
            self.refseq,
            self.altseq,
        )


@dataclass(eq=True, frozen=True)
class ExonSmallVariant(ExonPosition, SmallVariant):
    """Base class for exon variant data objects."""

    @classmethod
    def from_exon(
        cls, exon: ExonMappablePosition, refseq: str, altseq: str
    ) -> List[ExonSmallVariant]:
        """Convert an exon position plus ref/alt sequences into an exon variant data object."""
        raise NotImplementedError()


@dataclass(eq=True, frozen=True)
class ProteinSmallVariant(ProteinPosition, SmallVariant):
    """Base class for protein variant data objects."""

    @classmethod
    def from_cdna(
        cls, cdna: CdnaMappablePosition, protein: ProteinMappablePosition, refseq: str, altseq: str
    ) -> List[ProteinSmallVariant]:
        """Convert a cDNA position plus ref/alt sequences into a protein variant data object."""
        variant_list = []

        protein_refseq = protein.sequence()
        for cdna_alt in expand_nt(altseq):
            if is_frameshift(cdna.start, cdna.end, cdna_alt):
                # If the number of nucleotides added/deleted deletion is not divisible by 3, it results in a frame shift
                raise NotImplementedError()  # TODO
            else:
                for protein_alt in protein._data.mutate_cds_to_protein(
                    cdna.transcript_id, cdna.start, cdna.end, cdna_alt
                ):
                    # Remove bases that are unchanged between the ref and alt alleles
                    new_ref, new_alt, new_start, new_end = collapse_mutation(
                        protein_refseq, protein_alt
                    )
                    start = protein.start + new_start
                    end = protein.end - new_end

                    # Determine the type of variant
                    if len(new_ref) == 1 and len(new_alt) == 1:
                        # If exactly one base changes it is a substitution
                        variant = ProteinSubstitution.copy_from(
                            protein, start=start, end=end, refseq=new_ref, altseq=new_alt
                        )
                    elif len(new_ref) > 0 and len(new_alt) == 0:
                        # If only the reference sequence is removed it is a deletion
                        variant = ProteinDeletion.copy_from(
                            protein, start=start, end=end, refseq=new_ref, altseq=new_alt
                        )
                    elif len(new_ref) == 0 and len(new_alt) > 0:
                        # If only new bases are added it is an insertion
                        variant = ProteinInsertion.copy_from(
                            protein, start=start, end=end, refseq=new_ref, altseq=new_alt
                        )
                    elif new_alt == new_ref * 2:
                        # If the bases are a copy of the bases immediately 5' it is an duplication
                        variant = ProteinDuplication.copy_from(
                            protein, start=start, end=end, refseq=new_ref, altseq=new_alt
                        )
                    elif len(new_ref) * len(new_alt) > 1:
                        # If more than one ref and alt base is changed it is a delins
                        variant = ProteinDelins.copy_from(
                            protein, start=start, end=end, refseq=new_ref, altseq=new_alt
                        )
                    else:
                        raise NotImplementedError(
                            f"Unknown ref/alt '{new_ref}/{new_alt}' for {protein}"
                        )

                    variant_list.append(variant)

        return variant_list

    def to_cdna(self) -> List[CdnaSmallVariant]:
        return self._data._protein_to_cdna_variant(
            [self.transcript_id],
            self.start,
            self.start_offset,
            self.end,
            self.end_offset,
            [self.strand],
            self.refseq,
            self.altseq,
        )

    def to_dna(self) -> List[DnaSmallVariant]:
        return self._data._protein_to_dna_variant(
            [self.transcript_id],
            self.start,
            self.start_offset,
            self.end,
            self.end_offset,
            [self.strand],
            self.refseq,
            self.altseq,
        )

    def to_exon(self) -> List[ExonSmallVariant]:
        raise NotImplementedError()

    def to_protein(self) -> List[ProteinSmallVariant]:
        return self._data._protein_to_protein_variant(
            [self.transcript_id],
            self.start,
            self.start_offset,
            self.end,
            self.end_offset,
            [self.strand],
            self.refseq,
            self.altseq,
        )

    def to_rna(self) -> List[RnaSmallVariant]:
        return self._data._protein_to_rna_variant(
            [self.transcript_id],
            self.start,
            self.start_offset,
            self.end,
            self.end_offset,
            [self.strand],
            self.refseq,
            self.altseq,
        )


@dataclass(eq=True, frozen=True)
class RnaSmallVariant(RnaPosition, SmallVariant):
    """Base class for RNA variant data objects."""

    @classmethod
    def from_rna(cls, rna: RnaMappablePosition, refseq: str, altseq: str) -> List[RnaSmallVariant]:
        """Convert an RNA position plus ref/alt sequences into an RNA variant data object."""
        variant_list = []

        ref_annotated = rna.sequence()
        for ref, alt in product(expand_nt(refseq), expand_nt(altseq)):
            # Assert that the given ref matches the annotated one
            if ref != ref_annotated:
                raise ValueError(
                    f"Given ref allele '{ref}' does not match annotated ref allele '{ref_annotated}'"
                )

            # Remove bases that are unchanged between the ref and alt alleles
            new_ref, new_alt, new_start, new_end = collapse_mutation(ref, alt)
            start = rna.start + new_start
            end = rna.end - new_end

            # Determine the type of variant
            if len(new_ref) == 1 and len(new_alt) == 1:
                # If exactly one base changes it is a substitution
                variant = RnaSubstitution.copy_from(
                    rna, start=start, end=end, refseq=new_ref, altseq=new_alt
                )
            elif len(new_ref) > 0 and len(new_alt) == 0:
                # If only the reference sequence is removed it is a deletion
                variant = RnaDeletion.copy_from(
                    rna, start=start, end=end, refseq=new_ref, altseq=new_alt
                )
            elif len(new_ref) == 0 and len(new_alt) > 0:
                # If only new bases are added it is an insertion
                variant = RnaInsertion.copy_from(
                    rna, start=start, end=end, refseq=new_ref, altseq=new_alt
                )
            elif new_alt == new_ref * 2:
                # If the bases are a copy of the bases immediately 5' it is an duplication
                variant = RnaDuplication.copy_from(
                    rna, start=start, end=end, refseq=new_ref, altseq=new_alt
                )
            elif len(new_ref) * len(new_alt) > 1:
                # If more than one ref and alt base is changed it is a delins
                variant = RnaDelins.copy_from(
                    rna, start=start, end=end, refseq=new_ref, altseq=new_alt
                )
            else:
                raise NotImplementedError(f"Unknown ref/alt '{new_ref}/{new_alt}' for {rna}")

            variant_list.append(variant)

        return variant_list

    def to_cdna(self) -> List[CdnaSmallVariant]:
        return self._data._rna_to_cdna_variant(
            [self.transcript_id],
            self.start,
            self.start_offset,
            self.end,
            self.end_offset,
            [self.strand],
            self.refseq,
            self.altseq,
        )

    def to_dna(self) -> List[DnaSmallVariant]:
        return self._data._rna_to_dna_variant(
            [self.transcript_id],
            self.start,
            self.start_offset,
            self.end,
            self.end_offset,
            [self.strand],
            self.refseq,
            self.altseq,
        )

    def to_exon(self) -> List[ExonSmallVariant]:
        raise NotImplementedError()

    def to_protein(self) -> List[ProteinSmallVariant]:
        return self._data._rna_to_protein_variant(
            [self.transcript_id],
            self.start,
            self.start_offset,
            self.end,
            self.end_offset,
            [self.strand],
            self.refseq,
            self.altseq,
        )

    def to_rna(self) -> List[RnaSmallVariant]:
        return self._data._rna_to_rna_variant(
            [self.transcript_id],
            self.start,
            self.start_offset,
            self.end,
            self.end_offset,
            [self.strand],
            self.refseq,
            self.altseq,
        )


# -------------------------------------------------------------------------------------------------
# Deletion classes
# -------------------------------------------------------------------------------------------------
@dataclass(eq=True, frozen=True)
class Deletion(SmallVariant):
    """Base class for deletion variant data objects."""

    @property
    def type(self) -> str:
        return DELETION


@dataclass(eq=True, frozen=True)
class CdnaDeletion(Deletion, CdnaSmallVariant):
    """Stores information on a cDNA deletion variant and converts to other position types."""

    def __str__(self) -> str:
        start = format_hgvs_position(self.start, self.start_offset)
        end = format_hgvs_position(self.end, self.end_offset)
        if start == end:
            return f"{self.transcript_id}:c.{start}del"
        else:
            return f"{self.transcript_id}:c.{start}_{end}del"


@dataclass(eq=True, frozen=True)
class DnaDeletion(Deletion, DnaSmallVariant):
    """Stores information on a DNA deletion variant and converts to other position types."""

    def __str__(self) -> str:
        start = format_hgvs_position(self.start, self.start_offset)
        end = format_hgvs_position(self.end, self.end_offset)
        if start == end:
            return f"{self.contig_id}:g.{start}del"
        else:
            return f"{self.contig_id}:g.{start}_{end}del"


@dataclass(eq=True, frozen=True)
class ProteinDeletion(Deletion, ProteinSmallVariant):
    """Stores information on a protein deletion variant and converts to other position types."""

    def __str__(self) -> str:
        start = format_hgvs_position(self.start, self.start_offset)
        end = format_hgvs_position(self.end, self.end_offset)
        start_seq = self.refseq[0]
        end_seq = self.refseq[-1]
        if start == end:
            return f"{self.protein_id}:p.{start_seq}{start}del"
        else:
            return f"{self.protein_id}:p.{start_seq}{start}_{end_seq}{end}del"


@dataclass(eq=True, frozen=True)
class RnaDeletion(Deletion, RnaSmallVariant):
    """Stores information on an RNA deletion variant and converts to other position types."""

    def __str__(self) -> str:
        start = format_hgvs_position(self.start, self.start_offset)
        end = format_hgvs_position(self.end, self.end_offset)
        if start == end:
            return f"{self.transcript_id}:r.{start}del"
        else:
            return f"{self.transcript_id}:r.{start}_{end}del"


# -------------------------------------------------------------------------------------------------
# Delins classes
# -------------------------------------------------------------------------------------------------
@dataclass(eq=True, frozen=True)
class Delins(SmallVariant):
    """Base class for delins variant data objects."""

    @property
    def type(self) -> str:
        return DELINS


@dataclass(eq=True, frozen=True)
class CdnaDelins(Delins, CdnaSmallVariant):
    """Stores information on a cDNA delins variant and converts to other position types."""

    def __str__(self) -> str:
        start = format_hgvs_position(self.start, self.start_offset)
        end = format_hgvs_position(self.end, self.end_offset)
        if start == end:
            return f"{self.transcript_id}:c.{start}delins{self.altseq}"
        else:
            return f"{self.transcript_id}:c.{start}_{end}delins{self.altseq}"


@dataclass(eq=True, frozen=True)
class DnaDelins(Delins, DnaSmallVariant):
    """Stores information on a DNA delins variant and converts to other position types."""

    def __str__(self) -> str:
        start = format_hgvs_position(self.start, self.start_offset)
        end = format_hgvs_position(self.end, self.end_offset)
        if start == end:
            return f"{self.contig_id}:g.{start}delins{self.altseq}"
        else:
            return f"{self.contig_id}:g.{start}_{end}delins{self.altseq}"


@dataclass(eq=True, frozen=True)
class ProteinDelins(Delins, ProteinSmallVariant):
    """Stores information on a protein delins variant and converts to other position types."""

    def __str__(self) -> str:
        start = format_hgvs_position(self.start, self.start_offset)
        end = format_hgvs_position(self.end, self.end_offset)
        start_seq = self.refseq[0]
        end_seq = self.refseq[-1]
        if start == end:
            return f"{self.protein_id}:p.{start_seq}{start}delins{self.altseq}"
        else:
            return f"{self.protein_id}:p.{start_seq}{start}_{end_seq}{end}delins{self.altseq}"


@dataclass(eq=True, frozen=True)
class RnaDelins(Delins, RnaSmallVariant):
    """Stores information on an RNA delins variant and converts to other position types."""

    def __str__(self) -> str:
        start = format_hgvs_position(self.start, self.start_offset)
        end = format_hgvs_position(self.end, self.end_offset)
        if start == end:
            return f"{self.transcript_id}:r.{start}delins{self.altseq}"
        else:
            return f"{self.transcript_id}:r.{start}_{end}delins{self.altseq}"


# -------------------------------------------------------------------------------------------------
# Duplication classes
# -------------------------------------------------------------------------------------------------
@dataclass(eq=True, frozen=True)
class Duplication(SmallVariant):
    """Base class for duplication variant data objects."""

    @property
    def type(self) -> str:
        return DUPLICATION


@dataclass(eq=True, frozen=True)
class CdnaDuplication(Duplication, CdnaSmallVariant):
    """Stores information on a cDNA duplication variant and converts to other position types."""

    def __str__(self) -> str:
        start = format_hgvs_position(self.start, self.start_offset)
        end = format_hgvs_position(self.end, self.end_offset)
        if start == end:
            return f"{self.transcript_id}:c.{start}dup"
        else:
            return f"{self.transcript_id}:c.{start}_{end}dup"


@dataclass(eq=True, frozen=True)
class DnaDuplication(Duplication, DnaSmallVariant):
    """Stores information on a DNA duplication variant and converts to other position types."""

    def __str__(self) -> str:
        start = format_hgvs_position(self.start, self.start_offset)
        end = format_hgvs_position(self.end, self.end_offset)
        if start == end:
            return f"{self.contig_id}:g.{start}dup"
        else:
            return f"{self.contig_id}:g.{start}_{end}dup"


@dataclass(eq=True, frozen=True)
class ProteinDuplication(Duplication, ProteinSmallVariant):
    """Stores information on a protein duplication variant and converts to other position types."""

    def __str__(self) -> str:
        start = format_hgvs_position(self.start, self.start_offset)
        end = format_hgvs_position(self.end, self.end_offset)
        start_seq = self.refseq[0]
        end_seq = self.refseq[-1]
        if start == end:
            return f"{self.protein_id}:p.{start_seq}{start}dup"
        else:
            return f"{self.protein_id}:p.{start_seq}{start}_{end_seq}{end}dup"


@dataclass(eq=True, frozen=True)
class RnaDuplication(Duplication, RnaSmallVariant):
    """Stores information on an RNA duplication variant and converts to other position types."""

    def __str__(self) -> str:
        start = format_hgvs_position(self.start, self.start_offset)
        end = format_hgvs_position(self.end, self.end_offset)
        if start == end:
            return f"{self.transcript_id}:r.{start}dup"
        else:
            return f"{self.transcript_id}:r.{start}_{end}dup"


# -------------------------------------------------------------------------------------------------
# Frameshift classes
# -------------------------------------------------------------------------------------------------
@dataclass(eq=True, frozen=True)
class Frameshift(SmallVariant):
    """Base class for insertion variant data objects."""

    @property
    def type(self) -> str:
        return FRAMESHIFT


@dataclass(eq=True, frozen=True)
class ProteinFrameshift(Frameshift, ProteinSmallVariant):
    """Stores information on a protein frameshift variant and converts to other position types."""

    def to_cdna(self) -> List[CdnaSmallVariant]:
        raise NotImplementedError()  # TODO

    def to_dna(self) -> List[DnaSmallVariant]:
        raise NotImplementedError()  # TODO

    def to_exon(self) -> List[ExonSmallVariant]:
        raise NotImplementedError()  # TODO

    def to_protein(self) -> List[ProteinSmallVariant]:
        raise NotImplementedError()  # TODO

    def to_rna(self) -> List[RnaSmallVariant]:
        raise NotImplementedError()  # TODO


# -------------------------------------------------------------------------------------------------
# Insertion classes
# -------------------------------------------------------------------------------------------------
@dataclass(eq=True, frozen=True)
class Insertion(SmallVariant):
    """Base class for insertion variant data objects."""

    @property
    def type(self) -> str:
        return INSERTION


@dataclass(eq=True, frozen=True)
class CdnaInsertion(Insertion, CdnaSmallVariant):
    """Stores information on a cDNA insertion variant and converts to other position types."""


@dataclass(eq=True, frozen=True)
class DnaInsertion(Insertion, DnaSmallVariant):
    """Stores information on a DNA insertion variant and converts to other position types."""


@dataclass(eq=True, frozen=True)
class ProteinInsertion(Insertion, ProteinSmallVariant):
    """Stores information on a protein insertion variant and converts to other position types."""


@dataclass(eq=True, frozen=True)
class RnaInsertion(Insertion, RnaSmallVariant):
    """Stores information on an RNA insertion variant and converts to other position types."""


# -------------------------------------------------------------------------------------------------
# Substitution classes
# -------------------------------------------------------------------------------------------------
@dataclass(eq=True, frozen=True)
class Substitution:
    """Base class for substitution variant data objects."""

    refseq: str
    altseq: str

    @property
    def type(self) -> str:
        return SUBTITUTION


@dataclass(eq=True, frozen=True)
class CdnaSubstitution(Substitution, CdnaSmallVariant):
    """Stores information on a cDNA substitution variant and converts to other position types."""

    def __str__(self) -> str:
        start = format_hgvs_position(self.start, self.start_offset)
        return f"{self.transcript_id}:c.{start}{self.refseq}>{self.altseq}"


@dataclass(eq=True, frozen=True)
class DnaSubstitution(Substitution, DnaSmallVariant):
    """Stores information on a DNA substitution variant and converts to other position types."""

    def __str__(self) -> str:
        start = format_hgvs_position(self.start, self.start_offset)
        return f"{self.contig_id}:g.{start}{self.refseq}>{self.altseq}"


@dataclass(eq=True, frozen=True)
class ProteinSubstitution(Substitution, ProteinSmallVariant):
    """Stores information on a protein substitution variant and converts to other position types."""

    def __str__(self) -> str:
        start = format_hgvs_position(self.start, self.start_offset)
        return f"{self.protein_id}:p.{self.refseq}{start}{self.altseq}"


@dataclass(eq=True, frozen=True)
class RnaSubstitution(Substitution, RnaSmallVariant):
    """Stores information on an RNA substitution variant and converts to other position types."""

    def __str__(self) -> str:
        start = format_hgvs_position(self.start, self.start_offset)
        return f"{self.transcript_id}:r.{start}{self.refseq}>{self.altseq}"


# -------------------------------------------------------------------------------------------------
# Type hints
# -------------------------------------------------------------------------------------------------
PositionType = TypeVar(
    "PositionType",
    bound=Union[
        CdnaMappablePosition,
        DnaMappablePosition,
        ExonMappablePosition,
        ProteinMappablePosition,
        RnaMappablePosition,
    ],
)


# -------------------------------------------------------------------------------------------------
# Core classes and methods
# -------------------------------------------------------------------------------------------------
class CachedCore(type):
    """Metaclass that allows 'Core' objects, and objects that inherit from 'Core', to be treated as
    singletons.
    """

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
    """Core class that handles converting between position types and retrieving information on
    biological features.
    """

    def __init__(
        self,
        gtf: str,
        cds: List[str],
        dna: List[str],
        peptide: List[str],
        rna: List[str],
        canonical_transcript: str = "",
        contig_alias: str = "",
        exon_alias: str = "",
        gene_alias: str = "",
        protein_alias: str = "",
        transcript_alias: str = "",
    ):
        """
        Args:
            gtf: path to a GTF files with feature annotations
            cds: list of paths to a FASTA files of CDS sequences
            dna: list of paths to a FASTA files of DNA sequences
            peptide: list of paths to a FASTA files of peptide sequences
            rna: list of paths to a FASTA files of RNA sequences
            canonical_transcript: path to a text file of canonical transcript IDs
            contig_alias: path to a TSV file mapping contig IDs to their alias(es)
            exon_alias: path to a TSV file mapping contig IDs to their alias(es)
            gene_alias: path to a TSV file mapping contig IDs to their alias(es)
            protein_alias: path to a TSV file mapping contig IDs to their alias(es)
            transcript_alias: path to a TSV file mapping contig IDs to their alias(es)
        """
        self.df = read_gtf(gtf)
        self.cds = [read_fasta(i) for i in cds]
        self.dna = [read_fasta(i) for i in dna]
        self.peptide = [read_fasta(i) for i in peptide]
        self.rna = [read_fasta(i) for i in rna]
        self.canonical_transcript = txt_to_list(canonical_transcript)
        self.contig_alias = tsv_to_dict(contig_alias)
        self.exon_alias = tsv_to_dict(exon_alias)
        self.gene_alias = tsv_to_dict(gene_alias)
        self.protein_alias = tsv_to_dict(protein_alias)
        self.transcript_alias = tsv_to_dict(transcript_alias)

    # ---------------------------------------------------------------------------------------------
    # all_<feature_symbol>s
    # ---------------------------------------------------------------------------------------------
    def all_contig_ids(self) -> List[str]:
        """Return a list of all contig (chromosome) IDs."""
        return self._uniquify_series(self.df[CONTIG_ID])

    def all_exon_ids(self) -> List[str]:
        """Return a list of all exon IDs."""
        return self._uniquify_series(self.df[EXON_ID])

    def all_gene_ids(self) -> List[str]:
        """Return a list of all gene IDs."""
        return self._uniquify_series(self.df[GENE_ID])

    def all_gene_names(self) -> List[str]:
        """Return a list of all gene names."""
        return self._uniquify_series(self.df[GENE_NAME])

    def all_protein_ids(self) -> List[str]:
        """Return a list of all protein IDs."""
        return self._uniquify_series(self.df[PROTEIN_ID])

    def all_transcript_ids(self) -> List[str]:
        """Return a list of all transcript IDs."""
        return self._uniquify_series(self.df[TRANSCRIPT_ID])

    def all_transcript_names(self) -> List[str]:
        """Return a list of all transcript names."""
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
        """Given a feature symbol and a desired ID type, return the corresponding ID(s)."""
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
        """Given a contig ID, return a dataframe where all values have the ID."""
        return self._query(feature, CONTIG_ID)

    def _query_by_exon_id(self, feature: Union[List[str], str]) -> pd.DataFrame:
        """Given an exon ID, return a dataframe where all values have the ID."""
        return self._query(feature, EXON_ID)

    def _query_by_gene_id(self, feature: Union[List[str], str]) -> pd.DataFrame:
        """Given a gene ID, return a dataframe where all values have the ID."""
        return self._query(feature, GENE_ID)

    def _query_by_gene_name(self, feature: Union[List[str], str]) -> pd.DataFrame:
        """Given a gene name, return a dataframe where all values have the ID."""
        return self._query(feature, GENE_NAME)

    def _query_by_protein_id(self, feature: Union[List[str], str]) -> pd.DataFrame:
        """Given a protein ID, return a dataframe where all values have the ID."""
        return self._query(feature, PROTEIN_ID)

    def _query_by_transcript_id(self, feature: Union[List[str], str]) -> pd.DataFrame:
        """Given a transcript ID, return a dataframe where all values have the ID."""
        return self._query(feature, TRANSCRIPT_ID)

    def _query_by_transcript_name(self, feature: Union[List[str], str]) -> pd.DataFrame:
        """Given a transcript name, return a dataframe where all values have the ID."""
        return self._query(feature, TRANSCRIPT_NAME)

    def _query(self, feature: Union[List[str], str], col: str) -> pd.DataFrame:
        """Generic function for filtering the dataframe."""
        feature = [feature] if isinstance(feature, str) else feature
        sudbf = self.df.loc[self.df[col].isin(feature)]

        return sudbf

    # ---------------------------------------------------------------------------------------------
    # normalize_feature
    # ---------------------------------------------------------------------------------------------
    def normalize_feature(self, feature: str, feature_type: str = "") -> List[Tuple[str, str]]:
        """Normalize a feature ID to one or more representations."""
        normalized = []

        feature_type, result = self._normalize_feature(feature, feature_type=feature_type)
        if feature_type:
            normalized = self._uniquify_series(result[feature_type])

        return [(i, feature_type) for i in normalized]

    def _normalize_feature(self, feature: str, feature_type: str = "") -> Tuple[str, pd.DataFrame]:
        """Normalize a feature ID, return its type and a dataframe of matching values."""
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
        """Normalize a contig ID, return a dataframe of matching values."""
        featurel = [feature] + self.get_contig_alias(feature)

        return self._query_by_contig_id(featurel)

    def _normalize_exon_id(self, feature: str) -> pd.DataFrame:
        """Normalize an exon ID, return a dataframe of matching values."""
        featurel = [feature] + self.get_exon_alias(feature)

        return self._query_by_exon_id(featurel)

    def _normalize_gene_id(self, feature: str) -> pd.DataFrame:
        """Normalize a gene ID, return a dataframe of matching values."""
        featurel = [feature] + self.get_gene_alias(feature)

        return self._query_by_gene_id(featurel)

    def _normalize_gene_name(self, feature: str) -> pd.DataFrame:
        """Normalize a gene name, return a dataframe of matching values."""
        featurel = [feature] + self.get_gene_alias(feature)

        return self._query_by_gene_name(featurel)

    def _normalize_protein_id(self, feature: str) -> pd.DataFrame:
        """Normalize a protein ID, return a dataframe of matching values."""
        featurel = [feature] + self.get_protein_alias(feature)

        return self._query_by_protein_id(featurel)

    def _normalize_transcript_id(self, feature: str) -> pd.DataFrame:
        """Normalize a transcript ID, return a dataframe of matching values."""
        featurel = [feature] + self.get_transcript_alias(feature)

        return self._query_by_transcript_id(featurel)

    def _normalize_transcript_name(self, feature: str) -> pd.DataFrame:
        """Normalize a transcript name, return a dataframe of matching values."""
        featurel = [feature] + self.get_transcript_alias(feature)

        return self._query_by_transcript_name(featurel)

    # ---------------------------------------------------------------------------------------------
    # is_<feature>(feature)
    # ---------------------------------------------------------------------------------------------
    def is_contig(self, feature: str) -> bool:
        """Return True if the given feature is a contig ID."""
        return any((i[1] == CONTIG_ID for i in self.normalize_feature(feature, CONTIG_ID)))

    def is_exon(self, feature: str) -> bool:
        """Return True if the given feature is an exon ID."""
        return any((i[1] == EXON_ID for i in self.normalize_feature(feature, EXON_ID)))

    def is_gene(self, feature: str) -> bool:
        """Return True if the given feature is a gene ID or name."""
        return any((i[1] == GENE_ID for i in self.normalize_feature(feature, GENE_ID))) or any(
            (i[1] == GENE_NAME for i in self.normalize_feature(feature, GENE_NAME))
        )

    def is_protein(self, feature: str) -> bool:
        """Return True if the given feature is a protein ID."""
        return any((i[1] == PROTEIN_ID for i in self.normalize_feature(feature, PROTEIN_ID)))

    def is_transcript(self, feature: str) -> bool:
        """Return True if the given feature is a transcript ID or name."""
        return any(
            (i[1] == TRANSCRIPT_ID for i in self.normalize_feature(feature, TRANSCRIPT_ID))
        ) or any(
            (i[1] == TRANSCRIPT_NAME for i in self.normalize_feature(feature, TRANSCRIPT_NAME))
        )

    # ---------------------------------------------------------------------------------------------
    # <feature>_sequence
    # ---------------------------------------------------------------------------------------------
    def cds_sequence(
        self, transcript_id: str, start: Optional[int] = None, end: Optional[int] = None
    ) -> str:
        """Return the nucleotide sequence at the given CDS coordinates."""
        return self._get_sequence(self.cds, transcript_id, start=start, end=end)

    def dna_sequence(
        self,
        contig_id: str,
        start: Optional[int] = None,
        end: Optional[int] = None,
        strand: str = "+",
    ) -> str:
        """Return the nucleotide sequence at the given contig coordinates."""
        return self._get_sequence(self.dna, contig_id, start=start, end=end, strand=strand)

    def peptide_sequence(
        self, protein_id: str, start: Optional[int] = None, end: Optional[int] = None
    ) -> str:
        """Return the amino acid sequence at the given peptide coordinates."""
        return self._get_sequence(self.peptide, protein_id, start=start, end=end)

    def rna_sequence(
        self, transcript_id: str, start: Optional[int] = None, end: Optional[int] = None
    ) -> str:
        """Return the nucleotide sequence at the given cDNA or ncRNA coordinates."""
        return self._get_sequence(self.rna, transcript_id, start=start, end=end)

    def _get_sequence(
        self,
        fasta: List[Fasta],
        ref: str,
        start: Optional[int] = None,
        end: Optional[int] = None,
        strand: Optional[str] = None,
    ):
        """Return the sequence between the given positions (inclusive)."""
        for f in fasta:
            try:
                seq = f[ref]
                break
            except KeyError:
                continue
        else:
            raise KeyError(f"Sequence '{ref}' not found")

        # if no positions are given, return the whole sequence
        seqlen = len(seq)
        start = start if start is not None else 1
        end = end if end is not None else seqlen

        # validate that the given positions fall within the sequence
        if not (0 <= start < seqlen):
            raise ValueError(f"Start must be from 1 to {seqlen} ({start})")
        if not (1 < end <= seqlen):
            raise ValueError(f"End must be from 2 to {seqlen + 1} ({end})")

        # sanity check that the end position is after the start
        if end < start:
            raise ValueError(f"End must be >= start ({end} < {start})")

        subseq = seq[start - 1 : end]
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
        """Return True if the given transcript ID is a canonical transcript."""
        return feature in self.canonical_transcript

    # ---------------------------------------------------------------------------------------------
    # get_<feature>
    # ---------------------------------------------------------------------------------------------
    def get_cdna(self, feature: str) -> List[CdnaMappablePosition]:
        """Return the cDNA position(s) of the given feature."""
        result = []

        transcript_ids = self.transcript_ids(feature)
        mask = (self.df[TRANSCRIPT_ID].isin(transcript_ids)) & (self.df["feature"] == "cdna")
        for _, cdna in self.df[mask].iterrows():
            result.append(
                CdnaMappablePosition(
                    _data=self,
                    contig_id=cdna.contig_id,
                    start=cdna.cdna_start,
                    start_offset=0,
                    end=cdna.cdna_end,
                    end_offset=0,
                    strand=cdna.strand,
                    gene_id=cdna.gene_id,
                    gene_name=cdna.gene_name,
                    transcript_id=cdna.transcript_id,
                    transcript_name=cdna.transcript_name,
                    protein_id=cdna.protein_id,
                )
            )

        return sorted(set(result))

    def get_dna(self, feature: str) -> List[DnaMappablePosition]:
        """Return the DNA position(s) of the given feature."""
        result = []

        # get the strand of the original feature
        _, df = self._normalize_feature(feature)
        strand_list = self._uniquify_series(df["strand"])

        for contig_id in self.contig_ids(feature):
            for fasta in self.dna:
                try:
                    contig_seq = fasta[contig_id]
                    break
                except KeyError:
                    continue
            else:
                raise KeyError(f"Sequence '{contig_id}' not found")

            start = 1
            end = len(contig_seq)
            for strand in strand_list:
                result.append(
                    DnaMappablePosition(
                        _data=self,
                        contig_id=contig_id,
                        start=start,
                        start_offset=0,
                        end=end,
                        end_offset=0,
                        strand=strand,
                    )
                )

        return sorted(set(result))

    def get_exons(self, feature: str) -> List[ExonMappablePosition]:
        """Return the exon position(s) of the given feature."""
        result = []

        exon_ids = self.exon_ids(feature)
        mask = (self.df[EXON_ID].isin(exon_ids)) & (self.df["feature"] == "exon")
        for _, exon in self.df[mask].iterrows():
            result.append(
                ExonMappablePosition(
                    _data=self,
                    contig_id=exon.contig_id,
                    start=int(exon.exon_number),
                    start_offset=0,
                    end=int(exon.exon_number),
                    end_offset=0,
                    strand=exon.strand,
                    gene_id=exon.gene_id,
                    gene_name=exon.gene_name,
                    transcript_id=exon.transcript_id,
                    transcript_name=exon.transcript_name,
                    exon_id=exon.exon_id,
                )
            )

        return sorted(set(result))

    def get_genes(self, feature: str) -> List[DnaMappablePosition]:
        """Return the gene position(s) of the given feature."""
        result = []

        gene_ids = self.gene_ids(feature)
        mask = (self.df[GENE_ID].isin(gene_ids)) & (self.df["feature"] == "gene")
        for _, gene in self.df[mask].iterrows():
            result.append(
                DnaMappablePosition(
                    _data=self,
                    contig_id=gene.contig_id,
                    start=gene.start,
                    start_offset=0,
                    end=gene.end,
                    end_offset=0,
                    strand=gene.strand,
                )
            )

        return sorted(set(result))

    def get_transcripts(self, feature: str) -> List[RnaMappablePosition]:
        """Return the transcript position(s) of the given feature."""
        result = []

        transcript_ids = self.transcript_ids(feature)
        mask = (self.df[TRANSCRIPT_ID].isin(transcript_ids)) & (self.df["feature"] == "transcript")
        for _, transcript in self.df[mask].iterrows():
            result.append(
                RnaMappablePosition(
                    _data=self,
                    contig_id=transcript.contig_id,
                    start=transcript.transcript_start,
                    start_offset=0,
                    end=transcript.transcript_end,
                    end_offset=0,
                    strand=transcript.strand,
                    gene_id=transcript.gene_id,
                    gene_name=transcript.gene_name,
                    transcript_id=transcript.transcript_id,
                    transcript_name=transcript.transcript_name,
                )
            )

        return sorted(set(result))

    # ---------------------------------------------------------------------------------------------
    # <feature>_to_<feature>
    # ---------------------------------------------------------------------------------------------
    def cdna_to_cdna(
        self,
        feature: str,
        start: int,
        end: Optional[int] = None,
        strand: Optional[str] = None,
        start_offset: Optional[int] = None,
        end_offset: Optional[int] = None,
    ) -> List[CdnaMappablePosition]:
        """Map a cDNA position to zero or more cDNA positions."""
        return self._map(
            self.transcript_ids,
            self._cdna_to_cdna,
            feature,
            start,
            start_offset=start_offset,
            end=end,
            end_offset=end_offset,
            strand=strand,
        )

    def cdna_to_dna(
        self,
        feature: str,
        start: int,
        end: Optional[int] = None,
        strand: Optional[str] = None,
        start_offset: Optional[int] = None,
        end_offset: Optional[int] = None,
    ) -> List[DnaMappablePosition]:
        """Map a cDNA position to zero or more DNA positions."""
        return self._map(
            self.transcript_ids,
            self._cdna_to_dna,
            feature,
            start,
            start_offset=start_offset,
            end=end,
            end_offset=end_offset,
            strand=strand,
        )

    def cdna_to_exon(
        self,
        feature: str,
        start: int,
        end: Optional[int] = None,
        strand: Optional[str] = None,
        start_offset: Optional[int] = None,
        end_offset: Optional[int] = None,
    ) -> List[ExonMappablePosition]:
        """Map a cDNA position to an exon positions."""
        return self._map(
            self.transcript_ids,
            self._cdna_to_exon,
            feature,
            start,
            start_offset=start_offset,
            end=end,
            end_offset=end_offset,
            strand=strand,
        )

    def cdna_to_protein(
        self,
        feature: str,
        start: int,
        end: Optional[int] = None,
        strand: Optional[str] = None,
        start_offset: Optional[int] = None,
        end_offset: Optional[int] = None,
    ) -> List[ProteinMappablePosition]:
        """Map a cDNA position to zero or more protein positions."""
        return self._map(
            self.transcript_ids,
            self._cdna_to_protein,
            feature,
            start,
            start_offset=start_offset,
            end=end,
            end_offset=end_offset,
            strand=strand,
        )

    def cdna_to_rna(
        self,
        feature: str,
        start: int,
        end: Optional[int] = None,
        strand: Optional[str] = None,
        start_offset: Optional[int] = None,
        end_offset: Optional[int] = None,
    ) -> List[RnaMappablePosition]:
        """Map a cDNA position to zero or more RNA positions."""
        return self._map(
            self.transcript_ids,
            self._cdna_to_rna,
            feature,
            start,
            start_offset=start_offset,
            end=end,
            end_offset=end_offset,
            strand=strand,
        )

    def dna_to_cdna(
        self,
        feature: str,
        start: int,
        end: Optional[int] = None,
        strand: Optional[str] = None,
        start_offset: Optional[int] = None,
        end_offset: Optional[int] = None,
    ) -> List[CdnaMappablePosition]:
        """Map a DNA position to zero or more cDNA positions."""
        return self._map(
            self.contig_ids,
            self._dna_to_cdna,
            feature,
            start,
            start_offset=start_offset,
            end=end,
            end_offset=end_offset,
            strand=strand,
        )

    def dna_to_dna(
        self,
        feature: str,
        start: int,
        end: Optional[int] = None,
        strand: Optional[str] = None,
        start_offset: Optional[int] = None,
        end_offset: Optional[int] = None,
    ) -> List[DnaMappablePosition]:
        """Map a DNA position to zero or more DNA positions."""
        return self._map(
            self.contig_ids,
            self._dna_to_dna,
            feature,
            start,
            start_offset=start_offset,
            end=end,
            end_offset=end_offset,
            strand=strand,
        )

    def dna_to_exon(
        self,
        feature: str,
        start: int,
        end: Optional[int] = None,
        strand: Optional[str] = None,
        start_offset: Optional[int] = None,
        end_offset: Optional[int] = None,
    ) -> List[ExonMappablePosition]:
        """Map a DNA position to zero or more exon positions."""
        return self._map(
            self.contig_ids,
            self._dna_to_exon,
            feature,
            start,
            start_offset=start_offset,
            end=end,
            end_offset=end_offset,
            strand=strand,
        )

    def dna_to_protein(
        self,
        feature: str,
        start: int,
        end: Optional[int] = None,
        strand: Optional[str] = None,
        start_offset: Optional[int] = None,
        end_offset: Optional[int] = None,
    ) -> List[ProteinMappablePosition]:
        """Map a DNA position to zero or more protein positions."""
        return self._map(
            self.contig_ids,
            self._dna_to_protein,
            feature,
            start,
            start_offset=start_offset,
            end=end,
            end_offset=end_offset,
            strand=strand,
        )

    def dna_to_rna(
        self,
        feature: str,
        start: int,
        end: Optional[int] = None,
        strand: Optional[str] = None,
        start_offset: Optional[int] = None,
        end_offset: Optional[int] = None,
    ) -> List[RnaMappablePosition]:
        """Map a DNA position to zero or more RNA positions."""
        return self._map(
            self.contig_ids,
            self._dna_to_rna,
            feature,
            start,
            start_offset=start_offset,
            end=end,
            end_offset=end_offset,
            strand=strand,
        )

    def exon_to_cdna(
        self,
        feature: str,
        start: Optional[int] = None,
        end: Optional[int] = None,
        strand: Optional[str] = None,
        start_offset: Optional[int] = None,
        end_offset: Optional[int] = None,
    ) -> List[CdnaMappablePosition]:
        """Map an exon position to zero or more cDNA positions."""
        return self._map_exon(
            self._exon_to_cdna,
            feature,
            start=start,
            start_offset=start_offset,
            end=end,
            end_offset=end_offset,
            strand=strand,
        )

    def exon_to_dna(
        self,
        feature: str,
        start: Optional[int] = None,
        end: Optional[int] = None,
        strand: Optional[str] = None,
        start_offset: Optional[int] = None,
        end_offset: Optional[int] = None,
    ) -> List[DnaMappablePosition]:
        """Map an exon position to zero or more DNA positions."""
        return self._map_exon(
            self._exon_to_dna,
            feature,
            start=start,
            start_offset=start_offset,
            end=end,
            end_offset=end_offset,
            strand=strand,
        )

    def exon_to_exon(
        self,
        feature: str,
        start: Optional[int] = None,
        end: Optional[int] = None,
        strand: Optional[str] = None,
        start_offset: Optional[int] = None,
        end_offset: Optional[int] = None,
    ) -> List[ExonMappablePosition]:
        """Map an exon position to an exon positions."""
        return self._map_exon(
            self._exon_to_exon,
            feature,
            start=start,
            start_offset=start_offset,
            end=end,
            end_offset=end_offset,
            strand=strand,
        )

    def exon_to_protein(
        self,
        feature: str,
        start: Optional[int] = None,
        end: Optional[int] = None,
        strand: Optional[str] = None,
        start_offset: Optional[int] = None,
        end_offset: Optional[int] = None,
    ) -> List[ProteinMappablePosition]:
        """Map an exon position to zero or more protein positions."""
        return self._map_exon(
            self._exon_to_protein,
            feature,
            start=start,
            start_offset=start_offset,
            end=end,
            end_offset=end_offset,
            strand=strand,
        )

    def exon_to_rna(
        self,
        feature: str,
        start: Optional[int] = None,
        end: Optional[int] = None,
        strand: Optional[str] = None,
        start_offset: Optional[int] = None,
        end_offset: Optional[int] = None,
    ) -> List[RnaMappablePosition]:
        """Map an exon position to zero or more RNA positions."""
        return self._map_exon(
            self._exon_to_rna,
            feature,
            start=start,
            start_offset=start_offset,
            end=end,
            end_offset=end_offset,
            strand=strand,
        )

    def protein_to_cdna(
        self,
        feature: str,
        start: int,
        end: Optional[int] = None,
        strand: Optional[str] = None,
        start_offset: Optional[int] = None,
        end_offset: Optional[int] = None,
    ) -> List[CdnaMappablePosition]:
        """Map a protein position to zero or more cDNA positions."""
        return self._map(
            self.transcript_ids,
            self._protein_to_cdna,
            feature,
            start,
            start_offset=start_offset,
            end=end,
            end_offset=end_offset,
            strand=strand,
        )

    def protein_to_dna(
        self,
        feature: str,
        start: int,
        end: Optional[int] = None,
        strand: Optional[str] = None,
        start_offset: Optional[int] = None,
        end_offset: Optional[int] = None,
    ) -> List[DnaMappablePosition]:
        """Map a protein position to zero or more DNA positions."""
        return self._map(
            self.transcript_ids,
            self._protein_to_dna,
            feature,
            start,
            start_offset=start_offset,
            end=end,
            end_offset=end_offset,
            strand=strand,
        )

    def protein_to_exon(
        self,
        feature: str,
        start: int,
        end: Optional[int] = None,
        strand: Optional[str] = None,
        start_offset: Optional[int] = None,
        end_offset: Optional[int] = None,
    ) -> List[ExonMappablePosition]:
        """Map a protein position to zero or more exon positions."""
        return self._map(
            self.transcript_ids,
            self._protein_to_exon,
            feature,
            start,
            start_offset=start_offset,
            end=end,
            end_offset=end_offset,
            strand=strand,
        )

    def protein_to_protein(
        self,
        feature: str,
        start: int,
        end: Optional[int] = None,
        strand: Optional[str] = None,
        start_offset: Optional[int] = None,
        end_offset: Optional[int] = None,
    ) -> List[ProteinMappablePosition]:
        """Map a protein position to zero or more protein positions."""
        return self._map(
            self.transcript_ids,
            self._protein_to_protein,
            feature,
            start,
            start_offset=start_offset,
            end=end,
            end_offset=end_offset,
            strand=strand,
        )

    def protein_to_rna(
        self,
        feature: str,
        start: int,
        end: Optional[int] = None,
        strand: Optional[str] = None,
        start_offset: Optional[int] = None,
        end_offset: Optional[int] = None,
    ) -> List[RnaMappablePosition]:
        """Map a protein position to zero or more RNA positions."""
        return self._map(
            self.transcript_ids,
            self._protein_to_rna,
            feature,
            start,
            start_offset=start_offset,
            end=end,
            end_offset=end_offset,
            strand=strand,
        )

    def rna_to_cdna(
        self,
        feature: str,
        start: int,
        end: Optional[int] = None,
        strand: Optional[str] = None,
        start_offset: Optional[int] = None,
        end_offset: Optional[int] = None,
    ) -> List[CdnaMappablePosition]:
        """Map a RNA position to zero or more cDNA positions."""
        return self._map(
            self.transcript_ids,
            self._rna_to_cdna,
            feature,
            start,
            start_offset=start_offset,
            end=end,
            end_offset=end_offset,
            strand=strand,
        )

    def rna_to_dna(
        self,
        feature: str,
        start: int,
        end: Optional[int] = None,
        strand: Optional[str] = None,
        start_offset: Optional[int] = None,
        end_offset: Optional[int] = None,
    ) -> List[DnaMappablePosition]:
        """Map a RNA position to zero or more DNA positions."""
        return self._map(
            self.transcript_ids,
            self._rna_to_dna,
            feature,
            start,
            start_offset=start_offset,
            end=end,
            end_offset=end_offset,
            strand=strand,
        )

    def rna_to_exon(
        self,
        feature: str,
        start: int,
        end: Optional[int] = None,
        strand: Optional[str] = None,
        start_offset: Optional[int] = None,
        end_offset: Optional[int] = None,
    ) -> List[ExonMappablePosition]:
        """Map a RNA position to an exon positions."""
        return self._map(
            self.transcript_ids,
            self._rna_to_exon,
            feature,
            start,
            start_offset=start_offset,
            end=end,
            end_offset=end_offset,
            strand=strand,
        )

    def rna_to_protein(
        self,
        feature: str,
        start: int,
        end: Optional[int] = None,
        strand: Optional[str] = None,
        start_offset: Optional[int] = None,
        end_offset: Optional[int] = None,
    ) -> List[ProteinMappablePosition]:
        """Map a RNA position to zero or more protein positions."""
        return self._map(
            self.transcript_ids,
            self._rna_to_protein,
            feature,
            start,
            start_offset=start_offset,
            end=end,
            end_offset=end_offset,
            strand=strand,
        )

    def rna_to_rna(
        self,
        feature: str,
        start: int,
        end: Optional[int] = None,
        strand: Optional[str] = None,
        start_offset: Optional[int] = None,
        end_offset: Optional[int] = None,
    ) -> List[RnaMappablePosition]:
        """Map a RNA position to zero or more RNA positions."""
        return self._map(
            self.transcript_ids,
            self._rna_to_rna,
            feature,
            start,
            start_offset=start_offset,
            end=end,
            end_offset=end_offset,
            strand=strand,
        )

    def _cdna_to_cdna(
        self,
        transcript_ids: List[str],
        start: int,
        start_offset: int,
        end: int,
        end_offset: int,
        strand: List[str],
        include_stop: bool = True,
    ) -> List[CdnaMappablePosition]:
        feature = ["CDS", "stop_codon"] if include_stop else ["CDS"]

        def convert(position: int, offset: int):
            result = []

            # If the offset position can be normalized to a non-offset position, do so. Otherwise
            # just return an offset position.
            if offset:
                for dna in self._cdna_to_dna(
                    transcript_ids,
                    position,
                    offset,
                    position,
                    offset,
                    strand,
                    include_stop=include_stop,
                ):
                    for cdna in dna.to_cdna():
                        if cdna.transcript_id in transcript_ids:
                            result.append(cdna)

                if result:
                    return result

            mask = (
                (self.df[TRANSCRIPT_ID].isin(transcript_ids))
                & (self.df["cdna_start"] <= position)
                & (self.df["cdna_end"] >= position)
                & (self.df["strand"].isin(strand))
                & (self.df["feature"].isin(feature))
            )
            for _, cds in self.df[mask].iterrows():
                result.append(
                    CdnaMappablePosition(
                        _data=self,
                        contig_id=cds.contig_id,
                        start=start,
                        start_offset=offset,
                        end=end,
                        end_offset=offset,
                        strand=cds.strand,
                        gene_id=cds.gene_id,
                        gene_name=cds.gene_name,
                        transcript_id=cds.transcript_id,
                        transcript_name=cds.transcript_name,
                        protein_id=cds.protein_id,
                    )
                )

            return result

        result_start = convert(start, start_offset)
        if start == end:
            return sorted(result_start)
        else:
            result_end = convert(end, end_offset)
            return merge_positions(result_start, result_end, TRANSCRIPT_ID)

    def _cdna_to_cdna_variant(
        self,
        transcript_ids: List[str],
        start: int,
        start_offset: int,
        end: int,
        end_offset: int,
        strand: List[str],
        refseq: str,
        altseq: str,
        include_stop: bool = True,
    ) -> List[CdnaSmallVariant]:
        result: List[CdnaSmallVariant] = []

        for cdna in self._cdna_to_cdna(
            transcript_ids, start, start_offset, end, end_offset, strand, include_stop=include_stop
        ):
            result.extend(CdnaSmallVariant.from_cdna(cdna, refseq, altseq))

        return result

    def _cdna_to_dna(
        self,
        transcript_ids: List[str],
        start: int,
        start_offset: int,
        end: int,
        end_offset: int,
        strand: List[str],
        include_stop: bool = True,
    ) -> List[DnaMappablePosition]:
        feature = ["CDS", "stop_codon"] if include_stop else ["CDS"]

        def convert(position: int, offset: int):
            result = []

            mask = (
                (self.df[TRANSCRIPT_ID].isin(transcript_ids))
                & (self.df["cdna_start"] <= position)
                & (self.df["cdna_end"] >= position)
                & (self.df["strand"].isin(strand))
                & (self.df["feature"].isin(feature))
            )
            for _, cds in self.df[mask].iterrows():
                if cds.strand == "-":
                    new_start = new_end = cds.end - (position - cds.cdna_start) - offset
                else:
                    new_start = new_end = cds.start + (position - cds.cdna_start) + offset

                # TODO: Check that new new_start is actually on the contig
                result.append(
                    DnaMappablePosition(
                        _data=self,
                        contig_id=cds.contig_id,
                        start=new_start,
                        start_offset=0,
                        end=new_end,
                        end_offset=0,
                        strand=cds.strand,
                    )
                )

            return result

        result_start = convert(start, start_offset)
        if start == end:
            return sorted(result_start)
        else:
            result_end = convert(end, end_offset)
            return merge_positions(result_start, result_end, CONTIG_ID)

    def _cdna_to_dna_variant(
        self,
        transcript_ids: List[str],
        start: int,
        start_offset: int,
        end: int,
        end_offset: int,
        strand: List[str],
        refseq: str,
        altseq: str,
        include_stop: bool = True,
    ) -> List[DnaSmallVariant]:
        result: List[DnaSmallVariant] = []

        for dna in self._cdna_to_dna(
            transcript_ids, start, start_offset, end, end_offset, strand, include_stop=include_stop
        ):
            result.extend(DnaSmallVariant.from_dna(dna, refseq, altseq))

        return result

    def _cdna_to_exon(
        self,
        transcript_ids: List[str],
        start: int,
        start_offset: int,
        end: int,
        end_offset: int,
        strand: List[str],
        include_stop: bool = True,
    ) -> List[ExonMappablePosition]:
        feature = ["CDS", "stop_codon"] if include_stop else ["CDS"]

        def convert(position: int, offset: int):
            result = []

            # For an offset position, we need to calculate then equivalent DNA position then map
            # that to one or more exons. This is slower, so if there's no offset we can just map
            # directly to an exon.
            if offset:
                for dna in self._cdna_to_dna(
                    transcript_ids,
                    position,
                    offset,
                    position,
                    offset,
                    strand,
                    include_stop=include_stop,
                ):
                    for exon in dna.to_exon():
                        if exon.transcript_id in transcript_ids:
                            result.append(exon)

                return result

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
                for _, exon_row in self.df[mask_exon].iterrows():
                    result.append(
                        ExonMappablePosition(
                            _data=self,
                            contig_id=exon_row.contig_id,
                            start=int(exon_row.exon_number),
                            start_offset=0,
                            end=int(exon_row.exon_number),
                            end_offset=0,
                            strand=exon_row.strand,
                            gene_id=exon_row.gene_id,
                            gene_name=exon_row.gene_name,
                            transcript_id=exon_row.transcript_id,
                            transcript_name=exon_row.transcript_name,
                            exon_id=exon_row.exon_id,
                        )
                    )

            return result

        result_start = convert(start, start_offset)
        result_end = convert(end, end_offset)
        assert result_start == result_end  # TODO: mapping across introns?

        return sorted(result_start)

    def _cdna_to_exon_variant(
        self,
        transcript_ids: List[str],
        start: int,
        start_offset: int,
        end: int,
        end_offset: int,
        strand: List[str],
        refseq: str,
        altseq: str,
        include_stop: bool = True,
    ) -> List[ExonSmallVariant]:
        result: List[ExonSmallVariant] = []

        for exon in self._cdna_to_exon(
            transcript_ids, start, start_offset, end, end_offset, strand, include_stop=include_stop
        ):
            result.extend(ExonSmallVariant.from_exon(exon, refseq, altseq))

        return result

    def _cdna_to_protein(
        self,
        transcript_ids: List[str],
        start: int,
        start_offset: int,
        end: int,
        end_offset: int,
        strand: List[str],
    ) -> List[ProteinMappablePosition]:
        def convert(n: int) -> int:
            return floor((n - 1) / 3 + 1)

        result = []
        for cdna in self._cdna_to_cdna(
            transcript_ids, start, start_offset, end, end_offset, strand, include_stop=False
        ):
            # If the postion wasn't mapped to a non-offset position by `_cdna_to_cdna`, it means
            # that the position does not map to a protein.
            if cdna.start_offset or cdna.end_offset:
                continue

            # Convert the cDNA position to a protein position
            protein_start = convert(cdna.start)
            protein_end = convert(cdna.end)
            result.append(
                ProteinMappablePosition.copy_from(
                    cdna, start=protein_start, start_offset=0, end=protein_end, end_offset=0
                )
            )

        return result

    def _cdna_to_protein_variant(
        self,
        transcript_ids: List[str],
        start: int,
        start_offset: int,
        end: int,
        end_offset: int,
        strand: List[str],
        refseq: str,
        altseq: str,
    ) -> List[ProteinSmallVariant]:
        # TODO: A lot of this function is duplicated from _cdna_to_protein()
        def convert(n: int) -> int:
            return floor((n - 1) / 3 + 1)

        result: List[ProteinSmallVariant] = []
        for cdna in self._cdna_to_cdna(
            transcript_ids, start, start_offset, end, end_offset, strand, include_stop=False
        ):
            # If the postion wasn't mapped to a non-offset position by `_cdna_to_cdna`, it means
            # that the position does not map to a protein.
            if cdna.start_offset or cdna.end_offset:
                continue

            # Convert the cDNA position to a protein position
            protein_start = convert(cdna.start)
            protein_end = convert(cdna.end)
            protein = ProteinMappablePosition.copy_from(
                cdna, start=protein_start, start_offset=0, end=protein_end, end_offset=0
            )
            result.extend(ProteinSmallVariant.from_cdna(cdna, protein, refseq, altseq))

        return result

    def _cdna_to_rna(
        self,
        transcript_ids: List[str],
        start: int,
        start_offset: int,
        end: int,
        end_offset: int,
        strand: List[str],
        include_stop: bool = True,
    ) -> List[RnaMappablePosition]:
        feature = ["CDS", "stop_codon"] if include_stop else ["CDS"]

        def convert(position: int, offset: int):
            result = []

            # If the offset position can be normalized to a non-offset position, do so. Otherwise
            # just return an offset position.
            if offset:
                for dna in self._cdna_to_dna(
                    transcript_ids,
                    position,
                    offset,
                    position,
                    offset,
                    strand,
                    include_stop=include_stop,
                ):
                    for rna in dna.to_rna():
                        if rna.transcript_id in transcript_ids:
                            result.append(rna)

                if result:
                    return result

            mask = (
                (self.df[TRANSCRIPT_ID].isin(transcript_ids))
                & (self.df["cdna_start"] <= position)
                & (self.df["cdna_end"] >= position)
                & (self.df["strand"].isin(strand))
                & (self.df["feature"].isin(feature))
            )
            for _, cds in self.df[mask].iterrows():
                new_start = new_end = cds.transcript_start + (position - cds.cdna_start)
                result.append(
                    RnaMappablePosition(
                        _data=self,
                        contig_id=cds.contig_id,
                        start=new_start,
                        start_offset=offset,
                        end=new_end,
                        end_offset=offset,
                        strand=cds.strand,
                        gene_id=cds.gene_id,
                        gene_name=cds.gene_name,
                        transcript_id=cds.transcript_id,
                        transcript_name=cds.transcript_name,
                    )
                )

            return result

        result_start = convert(start, start_offset)
        if start == end:
            return sorted(result_start)
        else:
            result_end = convert(end, end_offset)
            return merge_positions(result_start, result_end, TRANSCRIPT_ID)

    def _cdna_to_rna_variant(
        self,
        transcript_ids: List[str],
        start: int,
        start_offset: int,
        end: int,
        end_offset: int,
        strand: List[str],
        refseq: str,
        altseq: str,
        include_stop: bool = True,
    ) -> List[RnaSmallVariant]:
        result: List[RnaSmallVariant] = []

        for rna in self._cdna_to_rna(
            transcript_ids, start, start_offset, end, end_offset, strand, include_stop=include_stop
        ):
            result.extend(RnaSmallVariant.from_rna(rna, refseq, altseq))

        return result

    def _dna_to_cdna(
        self,
        contig_ids: List[str],
        start: int,
        start_offset: int,
        end: int,
        end_offset: int,
        strand: List[str],
        include_stop: bool = True,
    ) -> List[CdnaMappablePosition]:
        feature = ["CDS", "stop_codon"] if include_stop else ["CDS"]

        def convert(position: int, offset: int):
            result = []

            for strand_ in strand:
                if strand_ == "-":
                    position_ = position - offset
                else:
                    position_ = position + offset

                mask = (
                    (self.df[CONTIG_ID].isin(contig_ids))
                    & (self.df["start"] <= position_)
                    & (self.df["end"] >= position_)
                    & (self.df["strand"] == strand_)
                    & (self.df["feature"].isin(feature))
                )
                for _, cds in self.df[mask].iterrows():
                    if cds.strand == "-":
                        new_start = new_end = cds.end - position_ + cds.cdna_start
                    else:
                        new_start = new_end = position_ - cds.start + cds.cdna_start

                    result.append(
                        CdnaMappablePosition(
                            _data=self,
                            contig_id=cds.contig_id,
                            start=new_start,
                            start_offset=0,
                            end=new_end,
                            end_offset=0,
                            strand=cds.strand,
                            gene_id=cds.gene_id,
                            gene_name=cds.gene_name,
                            transcript_id=cds.transcript_id,
                            transcript_name=cds.transcript_name,
                            protein_id=cds.protein_id,
                        )
                    )

            return result

        result_start = convert(start, start_offset)
        if start == end:
            return sorted(result_start)
        else:
            result_end = convert(end, end_offset)
            return merge_positions(result_start, result_end, TRANSCRIPT_ID)

    def _dna_to_cdna_variant(
        self,
        contig_ids: List[str],
        start: int,
        start_offset: int,
        end: int,
        end_offset: int,
        strand: List[str],
        refseq: str,
        altseq: str,
        include_stop: bool = True,
    ) -> List[CdnaSmallVariant]:
        result: List[CdnaSmallVariant] = []

        for cdna in self._dna_to_cdna(
            contig_ids, start, start_offset, end, end_offset, strand, include_stop=include_stop
        ):
            result.extend(CdnaSmallVariant.from_cdna(cdna, refseq, altseq))

        return result

    def _dna_to_dna(
        self,
        contig_ids: List[str],
        start: int,
        start_offset: int,
        end: int,
        end_offset: int,
        strand: List[str],
    ) -> List[DnaMappablePosition]:
        result = []

        # Sort the start and end positions after adjusting by offsets
        start, end = sorted([start, end])

        for contig_id_, strand_ in product(contig_ids, strand):
            if strand_ == "-":
                new_start = start - start_offset
                new_end = end - end_offset
            else:
                new_start = start + start_offset
                new_end = end + end_offset

            # TODO: Check that new new_start is actually on the contig
            result.append(
                DnaMappablePosition(
                    _data=self,
                    contig_id=contig_id_,
                    start=new_start,
                    start_offset=0,
                    end=new_end,
                    end_offset=0,
                    strand=strand_,
                )
            )

        return result

    def _dna_to_dna_variant(
        self,
        contig_ids: List[str],
        start: int,
        start_offset: int,
        end: int,
        end_offset: int,
        strand: List[str],
        refseq: str,
        altseq: str,
    ) -> List[DnaSmallVariant]:
        result: List[DnaSmallVariant] = []

        for dna in self._dna_to_dna(contig_ids, start, start_offset, end, end_offset, strand):
            result.extend(DnaSmallVariant.from_dna(dna, refseq, altseq))

        return result

    def _dna_to_exon(
        self,
        contig_ids: List[str],
        start: int,
        start_offset: int,
        end: int,
        end_offset: int,
        strand: List[str],
    ) -> List[ExonMappablePosition]:
        def convert(position: int, offset: int):
            result = []

            for strand_ in strand:
                if strand_ == "-":
                    position_ = position - offset
                else:
                    position_ = position + offset

                mask = (
                    (self.df[CONTIG_ID].isin(contig_ids))
                    & (self.df["start"] <= position_)
                    & (self.df["end"] >= position_)
                    & (self.df["strand"] == strand_)
                    & (self.df["feature"].isin(["exon"]))
                )
                for _, exon in self.df[mask].iterrows():
                    result.append(
                        ExonMappablePosition(
                            _data=self,
                            contig_id=exon.contig_id,
                            start=int(exon.exon_number),
                            start_offset=0,
                            end=int(exon.exon_number),
                            end_offset=0,
                            strand=exon.strand,
                            gene_id=exon.gene_id,
                            gene_name=exon.gene_name,
                            transcript_id=exon.transcript_id,
                            transcript_name=exon.transcript_name,
                            exon_id=exon.exon_id,
                        )
                    )

            return result

        result_start = convert(start, start_offset)
        result_end = convert(end, end_offset)
        assert result_start == result_end  # TODO: mapping across introns?

        return sorted(result_start)

    def _dna_to_exon_variant(
        self,
        contig_ids: List[str],
        start: int,
        start_offset: int,
        end: int,
        end_offset: int,
        strand: List[str],
        refseq: str,
        altseq: str,
    ) -> List[ExonSmallVariant]:
        result: List[ExonSmallVariant] = []

        for exon in self._dna_to_exon(contig_ids, start, start_offset, end, end_offset, strand):
            result.extend(ExonSmallVariant.from_exon(exon, refseq, altseq))

        return result

    def _dna_to_protein(
        self,
        contig_ids: List[str],
        start: int,
        start_offset: int,
        end: int,
        end_offset: int,
        strand: List[str],
    ) -> List[ProteinMappablePosition]:
        def convert(x):
            return floor((x - 1) / 3 + 1)

        protein = []
        for cdna in self._dna_to_cdna(
            contig_ids, start, start_offset, end, end_offset, strand, include_stop=False
        ):
            # Offset cDNA position are assumed to not map to a protein
            if cdna.start_offset or cdna.end_offset:
                continue

            pstart = convert(cdna.start)
            pend = convert(cdna.end)
            protein.append(
                ProteinMappablePosition.copy_from(
                    cdna, start=pstart, start_offset=0, end=pend, end_offset=0, _data=self
                )
            )

        return protein

    def _dna_to_protein_variant(
        self,
        contig_ids: List[str],
        start: int,
        start_offset: int,
        end: int,
        end_offset: int,
        strand: List[str],
        refseq: str,
        altseq: str,
    ) -> List[ProteinSmallVariant]:
        # TODO: A lot of this function is duplicated from _dna_to_protein()
        def convert(x):
            return floor((x - 1) / 3 + 1)

        result: List[ProteinSmallVariant] = []
        for cdna in self._dna_to_cdna(
            contig_ids, start, start_offset, end, end_offset, strand, include_stop=False
        ):
            # Offset cDNA position are assumed to not map to a protein
            if cdna.start_offset or cdna.end_offset:
                continue

            pstart = convert(cdna.start)
            pend = convert(cdna.end)
            protein = ProteinMappablePosition.copy_from(
                cdna, start=pstart, start_offset=0, end=pend, end_offset=0, _data=self
            )
            result.extend(ProteinSmallVariant.from_cdna(cdna, protein, refseq, altseq))

        return result

    def _dna_to_rna(
        self,
        contig_ids: List[str],
        start: int,
        start_offset: int,
        end: int,
        end_offset: int,
        strand: List[str],
    ) -> List[RnaMappablePosition]:
        def convert(position: int, offset: int):
            result = []

            for strand_ in strand:
                if strand_ == "-":
                    position_ = position - offset
                else:
                    position_ = position + offset

                mask = (
                    (self.df[CONTIG_ID].isin(contig_ids))
                    & (self.df["start"] <= position_)
                    & (self.df["end"] >= position_)
                    & (self.df["strand"] == strand_)
                    & (self.df["feature"] == "exon")
                )
                for _, exon in self.df[mask].iterrows():
                    if exon.strand == "-":
                        new_start = new_end = exon.end - position_ + exon.transcript_start
                    else:
                        new_start = new_end = position_ - exon.start + exon.transcript_start

                    result.append(
                        RnaMappablePosition(
                            _data=self,
                            contig_id=exon.contig_id,
                            start=new_start,
                            start_offset=0,
                            end=new_end,
                            end_offset=0,
                            strand=exon.strand,
                            gene_id=exon.gene_id,
                            gene_name=exon.gene_name,
                            transcript_id=exon.transcript_id,
                            transcript_name=exon.transcript_name,
                        )
                    )

            return result

        result_start = convert(start, start_offset)
        if start == end:
            return sorted(result_start)
        else:
            result_end = convert(end, end_offset)
            return merge_positions(result_start, result_end, TRANSCRIPT_ID)

    def _dna_to_rna_variant(
        self,
        contig_ids: List[str],
        start: int,
        start_offset: int,
        end: int,
        end_offset: int,
        strand: List[str],
        refseq: str,
        altseq: str,
    ) -> List[RnaSmallVariant]:
        result: List[RnaSmallVariant] = []

        for rna in self._dna_to_rna(contig_ids, start, start_offset, end, end_offset, strand):
            result.extend(RnaSmallVariant.from_rna(rna, refseq, altseq))

        return result

    def _exon_to_cdna(
        self,
        transcript_ids: List[str],
        start: int,
        start_offset: int,
        end: int,
        end_offset: int,
        strand: List[str],
        include_stop: bool = True,
    ) -> List[CdnaMappablePosition]:
        feature = ["CDS", "stop_codon"] if include_stop else ["CDS"]

        def convert(position: int, offset: int):
            result = []

            # TODO: Is there a reasonable case where an exon position would have an offset?
            assert not offset, offset

            mask = (
                (self.df[TRANSCRIPT_ID].isin(transcript_ids))
                & (self.df["exon_number"] == str(position))  # TODO: exon number should be int
                & (self.df["strand"].isin(strand))
                & (self.df["feature"].isin(feature))
            )
            for _, cds in self.df[mask].iterrows():
                result.append(
                    CdnaMappablePosition(
                        _data=self,
                        contig_id=cds.contig_id,
                        start=cds.cdna_start,
                        start_offset=0,
                        end=cds.cdna_end,
                        end_offset=0,
                        strand=cds.strand,
                        gene_id=cds.gene_id,
                        gene_name=cds.gene_name,
                        transcript_id=cds.transcript_id,
                        transcript_name=cds.transcript_name,
                        protein_id=cds.protein_id,
                    )
                )

            return result

        result_start = convert(start, start_offset)
        if start == end:
            return sorted(result_start)
        else:
            result_end = convert(end, end_offset)
            return merge_positions(result_start, result_end, TRANSCRIPT_ID)

    def _exon_to_dna(
        self,
        transcript_ids: List[str],
        start: int,
        start_offset: int,
        end: int,
        end_offset: int,
        strand: List[str],
    ) -> List[DnaMappablePosition]:
        def convert(position: int, offset: int):
            result = []

            # TODO: Is there a reasonable case where an exon position would have an offset?
            assert not offset, offset

            mask = (
                (self.df[TRANSCRIPT_ID].isin(transcript_ids))
                # TODO: exon number should be int
                & (self.df["exon_number"] == str(position))
                & (self.df["strand"].isin(strand))
                & (self.df["feature"] == "exon")
            )
            for _, exon in self.df[mask].iterrows():
                result.append(
                    DnaMappablePosition(
                        _data=self,
                        contig_id=exon.contig_id,
                        start=exon.start,
                        start_offset=0,
                        end=exon.end,
                        end_offset=0,
                        strand=exon.strand,
                    )
                )

            return result

        result_start = convert(start, start_offset)
        if start == end:
            return sorted(result_start)
        else:
            result_end = convert(end, end_offset)
            return merge_positions(result_start, result_end, CONTIG_ID)

    def _exon_to_exon(
        self,
        transcript_ids: List[str],
        start: int,
        start_offset: int,
        end: int,
        end_offset: int,
        strand: List[str],
    ) -> List[ExonMappablePosition]:
        def convert(position: int, offset):
            result = []

            # TODO: Is there a reasonable case where an exon position would have an offset?
            assert not offset, offset

            mask = (
                (self.df[TRANSCRIPT_ID].isin(transcript_ids))
                # TODO: exon number should be int
                & (self.df["exon_number"] == str(position))
                & (self.df["strand"].isin(strand))
                & (self.df["feature"] == "exon")
            )
            for _, exon in self.df[mask].iterrows():
                result.append(
                    ExonMappablePosition(
                        _data=self,
                        contig_id=exon.contig_id,
                        start=int(exon.exon_number),
                        start_offset=0,
                        end=int(exon.exon_number),
                        end_offset=0,
                        strand=exon.strand,
                        gene_id=exon.gene_id,
                        gene_name=exon.gene_name,
                        transcript_id=exon.transcript_id,
                        transcript_name=exon.transcript_name,
                        exon_id=exon.exon_id,
                    )
                )

            return result

        result_start = convert(start, start_offset)
        if start == end:
            return sorted(result_start)
        else:
            result_end = convert(end, end_offset)
            return merge_positions(result_start, result_end, TRANSCRIPT_ID)

    def _exon_to_protein(
        self,
        transcript_ids: List[str],
        start: int,
        start_offset: int,
        end: int,
        end_offset: int,
        strand: List[str],
    ) -> List[ProteinMappablePosition]:
        result = []
        for cdna in self._exon_to_cdna(
            transcript_ids, start, start_offset, end, end_offset, strand, include_stop=False
        ):
            result.extend(cdna.to_protein())

        return result

    def _exon_to_rna(
        self,
        transcript_ids: List[str],
        start: int,
        start_offset: int,
        end: int,
        end_offset: int,
        strand: List[str],
    ) -> List[RnaMappablePosition]:
        def convert(position: int, offset: int):
            result = []

            # TODO: Is there a reasonable case where an exon position would have an offset?
            assert not offset, offset

            mask = (
                (self.df[TRANSCRIPT_ID].isin(transcript_ids))
                & (self.df["exon_number"] == str(position))  # TODO: exon number should be int
                & (self.df["strand"].isin(strand))
                & (self.df["feature"] == "exon")
            )
            for _, exon in self.df[mask].iterrows():
                result.append(
                    RnaMappablePosition(
                        _data=self,
                        contig_id=exon.contig_id,
                        start=exon.transcript_start,
                        start_offset=0,
                        end=exon.transcript_end,
                        end_offset=0,
                        strand=exon.strand,
                        gene_id=exon.gene_id,
                        gene_name=exon.gene_name,
                        transcript_id=exon.transcript_id,
                        transcript_name=exon.transcript_name,
                    )
                )

            return result

        result_start = convert(start, start_offset)
        if start == end:
            return sorted(result_start)
        else:
            result_end = convert(end, end_offset)
            return merge_positions(result_start, result_end, TRANSCRIPT_ID)

    def _protein_to_cdna(
        self,
        transcript_ids: List[str],
        start: int,
        start_offset: int,
        end: int,
        end_offset: int,
        strand: List[str],
    ) -> List[CdnaMappablePosition]:
        def convert(position: int):
            return ((position - 1) * 3) + 1

        # TODO: Is there a reasonable case where an protein position would have an offset?
        assert not start_offset, start_offset
        assert not end_offset, end_offset

        cdna_start = convert(start)
        cdna_end = convert(end) + 2

        return self._cdna_to_cdna(
            transcript_ids, cdna_start, 0, cdna_end, 0, strand, include_stop=False
        )

    def _protein_to_cdna_variant(
        self,
        transcript_ids: List[str],
        start: int,
        start_offset: int,
        end: int,
        end_offset: int,
        strand: List[str],
        refseq: str,
        altseq: str,
    ) -> List[CdnaSmallVariant]:
        result: List[CdnaSmallVariant] = []

        for cdna in self._protein_to_cdna(
            transcript_ids, start, start_offset, end, end_offset, strand
        ):
            result.extend(CdnaSmallVariant.from_protein(cdna, refseq, altseq))

        return result

    def _protein_to_dna(
        self,
        transcript_ids: List[str],
        start: int,
        start_offset: int,
        end: int,
        end_offset: int,
        strand: List[str],
    ) -> List[DnaMappablePosition]:
        result = []
        for cdna in self._protein_to_cdna(
            transcript_ids, start, start_offset, end, end_offset, strand
        ):
            result.extend(cdna.to_dna())

        return sorted(set(result))

    def _protein_to_dna_variant(
        self,
        transcript_ids: List[str],
        start: int,
        start_offset: int,
        end: int,
        end_offset: int,
        strand: List[str],
        refseq: str,
        altseq: str,
    ) -> List[DnaSmallVariant]:
        result = []
        for cdna in self._protein_to_cdna_variant(
            transcript_ids, start, start_offset, end, end_offset, strand, refseq, altseq
        ):
            result.extend(cdna.to_dna())

        return sorted(set(result))

    def _protein_to_exon(
        self,
        transcript_ids: List[str],
        start: int,
        start_offset: int,
        end: int,
        end_offset: int,
        strand: List[str],
    ) -> List[ExonMappablePosition]:
        result = []
        for cdna in self._protein_to_cdna(
            transcript_ids, start, start_offset, end, end_offset, strand
        ):
            result.extend(cdna.to_exon())

        return sorted(set(result))

    def _protein_to_exon_variant(
        self,
        transcript_ids: List[str],
        start: int,
        start_offset: int,
        end: int,
        end_offset: int,
        strand: List[str],
        refseq: str,
        altseq: str,
    ) -> List[ExonSmallVariant]:
        result = []
        for cdna in self._protein_to_cdna_variant(
            transcript_ids, start, start_offset, end, end_offset, strand, refseq, altseq
        ):
            result.extend(cdna.to_exon())

        return sorted(set(result))

    def _protein_to_protein(
        self,
        transcript_ids: List[str],
        start: int,
        start_offset: int,
        end: int,
        end_offset: int,
        strand: List[str],
    ) -> List[ProteinMappablePosition]:
        result = []
        for cdna in self._protein_to_cdna(
            transcript_ids, start, start_offset, end, end_offset, strand
        ):
            result.extend(cdna.to_protein())

        return sorted(set(result))

    def _protein_to_protein_variant(
        self,
        transcript_ids: List[str],
        start: int,
        start_offset: int,
        end: int,
        end_offset: int,
        strand: List[str],
        refseq: str,
        altseq: str,
    ) -> List[ProteinSmallVariant]:
        result = []
        for cdna in self._protein_to_cdna_variant(
            transcript_ids, start, start_offset, end, end_offset, strand, refseq, altseq
        ):
            result.extend(cdna.to_protein())

        return sorted(set(result))

    def _protein_to_rna(
        self,
        transcript_ids: List[str],
        start: int,
        start_offset: int,
        end: int,
        end_offset: int,
        strand: List[str],
    ) -> List[RnaMappablePosition]:
        result = []
        for cdna in self._protein_to_cdna(
            transcript_ids, start, start_offset, end, end_offset, strand
        ):
            result.extend(cdna.to_rna())

        return sorted(set(result))

    def _protein_to_rna_variant(
        self,
        transcript_ids: List[str],
        start: int,
        start_offset: int,
        end: int,
        end_offset: int,
        strand: List[str],
        refseq: str,
        altseq: str,
    ) -> List[RnaSmallVariant]:
        result = []
        for cdna in self._protein_to_cdna_variant(
            transcript_ids, start, start_offset, end, end_offset, strand, refseq, altseq
        ):
            result.extend(cdna.to_rna())

        return sorted(set(result))

    def _rna_to_cdna(
        self,
        transcript_ids: List[str],
        start: int,
        start_offset: int,
        end: int,
        end_offset: int,
        strand: List[str],
        include_stop: bool = True,
    ) -> List[CdnaMappablePosition]:
        feature = ["CDS", "stop_codon"] if include_stop else ["CDS"]

        def convert(position: int, offset: int):
            result = []

            # If the offset position can be normalized to a non-offset position, do so. Otherwise
            # just return an offset position.
            if offset:
                for dna in self._rna_to_dna(
                    transcript_ids, position, offset, position, offset, strand
                ):
                    for cdna in dna.to_cdna():
                        if cdna.transcript_id in transcript_ids:
                            result.append(cdna)

                if result:
                    return result

            mask = (
                (self.df[TRANSCRIPT_ID].isin(transcript_ids))
                & (self.df["transcript_start"] <= position)
                & (self.df["transcript_end"] >= position)
                & (self.df["strand"].isin(strand))
                & (self.df["feature"].isin(feature))
            )
            for _, cds in self.df[mask].iterrows():
                new_start = new_end = cds.cdna_start + (position - cds.transcript_start)
                result.append(
                    CdnaMappablePosition(
                        _data=self,
                        contig_id=cds.contig_id,
                        start=new_start,
                        start_offset=offset,
                        end=new_end,
                        end_offset=offset,
                        strand=cds.strand,
                        gene_id=cds.gene_id,
                        gene_name=cds.gene_name,
                        transcript_id=cds.transcript_id,
                        transcript_name=cds.transcript_name,
                        protein_id=cds.protein_id,
                    )
                )

            return result

        result_start = convert(start, start_offset)
        if start == end:
            return sorted(result_start)
        else:
            result_end = convert(end, end_offset)
            return merge_positions(result_start, result_end, TRANSCRIPT_ID)

    def _rna_to_cdna_variant(
        self,
        transcript_ids: List[str],
        start: int,
        start_offset: int,
        end: int,
        end_offset: int,
        strand: List[str],
        refseq: str,
        altseq: str,
        include_stop: bool = True,
    ) -> List[CdnaSmallVariant]:
        result: List[CdnaSmallVariant] = []

        for cdna in self._rna_to_cdna(
            transcript_ids, start, start_offset, end, end_offset, strand, include_stop=include_stop
        ):
            result.extend(CdnaSmallVariant.from_cdna(cdna, refseq, altseq))

        return result

    def _rna_to_dna(
        self,
        transcript_ids: List[str],
        start: int,
        start_offset: int,
        end: int,
        end_offset: int,
        strand: List[str],
    ) -> List[DnaMappablePosition]:
        def convert(position: int, offset: int):
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
                if exon.strand == "-":
                    new_start = new_end = exon.end - (position - exon.transcript_start) - offset
                else:
                    new_start = new_end = exon.start + (position - exon.transcript_start) + offset

                # TODO: Check that new new_start is actually on the contig
                result.append(
                    DnaMappablePosition(
                        _data=self,
                        contig_id=exon.contig_id,
                        start=new_start,
                        start_offset=0,
                        end=new_end,
                        end_offset=0,
                        strand=exon.strand,
                    )
                )

            return result

        result_start = convert(start, start_offset)
        if start == end:
            return sorted(result_start)
        else:
            result_end = convert(end, end_offset)
            return merge_positions(result_start, result_end, CONTIG_ID)

    def _rna_to_dna_variant(
        self,
        transcript_ids: List[str],
        start: int,
        start_offset: int,
        end: int,
        end_offset: int,
        strand: List[str],
        refseq: str,
        altseq: str,
        include_stop: bool = True,
    ) -> List[DnaSmallVariant]:
        result: List[DnaSmallVariant] = []

        for dna in self._rna_to_dna(transcript_ids, start, start_offset, end, end_offset, strand):
            result.extend(DnaSmallVariant.from_dna(dna, refseq, altseq))

        return result

    def _rna_to_exon(
        self,
        transcript_ids: List[str],
        start: int,
        start_offset: int,
        end: int,
        end_offset: int,
        strand: List[str],
    ) -> List[ExonMappablePosition]:
        def convert(position: int, offset: int):
            result = []

            # For an offset position, we need to calculate then equivalent DNA position then map
            # that to one or more exons. This is slower, so if there's no offset we can just map
            # directly to an exon.
            if offset:
                for dna in self._rna_to_dna(
                    transcript_ids, position, offset, position, offset, strand
                ):
                    for exon in dna.to_exon():
                        if exon.transcript_id in transcript_ids:
                            result.append(exon)

                return result

            mask = (
                (self.df[TRANSCRIPT_ID].isin(transcript_ids))
                & (self.df["transcript_start"] <= position)
                & (self.df["transcript_end"] >= position)
                & (self.df["strand"].isin(strand))
                & (self.df["feature"].isin(["exon"]))
            )
            for _, exon_row in self.df[mask].iterrows():
                result.append(
                    ExonMappablePosition(
                        _data=self,
                        contig_id=exon_row.contig_id,
                        start=int(exon_row.exon_number),
                        start_offset=0,
                        end=int(exon_row.exon_number),
                        end_offset=0,
                        strand=exon_row.strand,
                        gene_id=exon_row.gene_id,
                        gene_name=exon_row.gene_name,
                        transcript_id=exon_row.transcript_id,
                        transcript_name=exon_row.transcript_name,
                        exon_id=exon_row.exon_id,
                    )
                )

            return result

        result_start = convert(start, start_offset)
        if start == end:
            return sorted(result_start)
        else:
            result_end = convert(end, end_offset)
            return merge_positions(result_start, result_end, TRANSCRIPT_ID)

    def _rna_to_exon_variant(
        self,
        transcript_ids: List[str],
        start: int,
        start_offset: int,
        end: int,
        end_offset: int,
        strand: List[str],
        refseq: str,
        altseq: str,
        include_stop: bool = True,
    ) -> List[ExonSmallVariant]:
        result: List[ExonSmallVariant] = []

        for exon in self._rna_to_exon(transcript_ids, start, start_offset, end, end_offset, strand):
            result.extend(ExonSmallVariant.from_exon(exon, refseq, altseq))

        return result

    def _rna_to_protein(
        self,
        transcript_ids: List[str],
        start: int,
        start_offset: int,
        end: int,
        end_offset: int,
        strand: List[str],
    ) -> List[ProteinMappablePosition]:
        result = []
        for cdna in self._rna_to_cdna(
            transcript_ids, start, start_offset, end, end_offset, strand, include_stop=False
        ):
            result.extend(cdna.to_protein())

        return sorted(set(result))

    def _rna_to_protein_variant(
        self,
        transcript_ids: List[str],
        start: int,
        start_offset: int,
        end: int,
        end_offset: int,
        strand: List[str],
        refseq: str,
        altseq: str,
    ) -> List[ProteinSmallVariant]:
        result = []
        for cdna in self._rna_to_cdna_variant(
            transcript_ids,
            start,
            start_offset,
            end,
            end_offset,
            strand,
            refseq,
            altseq,
            include_stop=False,
        ):
            result.extend(cdna.to_protein())

        return sorted(set(result))

    def _rna_to_rna(
        self,
        transcript_ids: List[str],
        start: int,
        start_offset: int,
        end: int,
        end_offset: int,
        strand: List[str],
    ) -> List[RnaMappablePosition]:
        def convert(position: int, offset: int):
            result = []

            # If the offset position can be normalized to a non-offset position, do so. Otherwise
            # just return an offset position.
            if offset:
                for dna in self._rna_to_dna(
                    transcript_ids, position, offset, position, offset, strand
                ):
                    for rna in dna.to_rna():
                        if rna.transcript_id in transcript_ids:
                            result.append(rna)

                if result:
                    return result

            mask = (
                (self.df[TRANSCRIPT_ID].isin(transcript_ids))
                & (self.df["transcript_start"] <= position)
                & (self.df["transcript_end"] >= position)
                & (self.df["strand"].isin(strand))
                & (self.df["feature"] == "exon")
            )
            for _, exon in self.df[mask].iterrows():
                result.append(
                    RnaMappablePosition(
                        _data=self,
                        contig_id=exon.contig_id,
                        start=start,
                        start_offset=offset,
                        end=end,
                        end_offset=offset,
                        strand=exon.strand,
                        gene_id=exon.gene_id,
                        gene_name=exon.gene_name,
                        transcript_id=exon.transcript_id,
                        transcript_name=exon.transcript_name,
                    )
                )

            return result

        result_start = convert(start, start_offset)
        if start == end:
            return sorted(result_start)
        else:
            result_end = convert(end, end_offset)
            return merge_positions(result_start, result_end, TRANSCRIPT_ID)

    def _rna_to_rna_variant(
        self,
        transcript_ids: List[str],
        start: int,
        start_offset: int,
        end: int,
        end_offset: int,
        strand: List[str],
        refseq: str,
        altseq: str,
        include_stop: bool = True,
    ) -> List[RnaSmallVariant]:
        result: List[RnaSmallVariant] = []

        for rna in self._rna_to_rna(transcript_ids, start, start_offset, end, end_offset, strand):
            result.extend(RnaSmallVariant.from_rna(rna, refseq, altseq))

        return result

    def _map(
        self,
        idfunc: Callable,
        mapfunc: Callable,
        feature: str,
        start: int,
        start_offset: Optional[int] = None,
        end: Optional[int] = None,
        end_offset: Optional[int] = None,
        strand: Optional[str] = None,
    ) -> List:
        end = end if end is not None else start
        start_offset = start_offset or 0
        end_offset = end_offset or 0
        strandl = [strand] if strand is not None else ["+", "-"]
        featurel = idfunc(feature)
        result = mapfunc(featurel, start, start_offset, end, end_offset, strandl)

        return sorted(set(result))

    def _map_exon(
        self,
        function: Callable,
        feature: str,
        start: Optional[int],
        start_offset: Optional[int],
        end: Optional[int],
        end_offset: Optional[int],
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
            result.extend(
                function([transcript_id], start2, start_offset, end2, end_offset, [strand2])
            )

        return sorted(set(result))

    def _get_exon_number(self, exon_id: str) -> List[Tuple[str, int, int, str]]:
        result = []

        mask = (self.df[EXON_ID] == exon_id) & (self.df["feature"] == "exon")
        for _, exon in self.df[mask].iterrows():
            result.append((exon.transcript_id, exon.exon_number, exon.exon_number, exon.strand))

        return sorted(set(result))

    # ---------------------------------------------------------------------------------------------
    # Utility functions
    # ---------------------------------------------------------------------------------------------
    def cds_offset(self, transcript_id: str) -> int:
        """Get the integer offset of the CDS from the start of the spliced RNA."""
        rna = self.cdna_to_rna(transcript_id, 1)
        assert len(rna) == 1, f"Ambiguous transcript start for '{transcript_id}'"
        offset = rna[0].start - 1
        assert offset >= 0

        return offset

    def mutate_cds_to_protein(
        self, transcript_id: str, cdna_start: int, cdna_end: int, cdna_alt: str
    ) -> List[str]:
        """Return the mutated protein sequence, given a cDNA position and alt allele."""
        pep_altseq_set = set()

        # If no alt, assume we're talking about a deletion variant
        if not cdna_alt:
            return [""]

        # Assert that the nucleotide change does not result in a frameshift
        assert not is_frameshift(cdna_start, cdna_end, cdna_alt)

        # Get the codon sequence
        codon_start_offset = (cdna_start - 1) % 3
        codon_start = cdna_start - codon_start_offset
        codon_end_offset = 2 - ((cdna_end - 1) % 3)
        codon_end = cdna_end + codon_end_offset
        codon_refseq = self.cds_sequence(transcript_id, codon_start, codon_end)

        # Mutate the codon sequence
        codon_refseq_left = codon_refseq[:codon_start_offset]
        codon_refseq_right = codon_refseq[-codon_end_offset:] if codon_end_offset else ""
        for i in expand_nt(cdna_alt):
            codon_altseq = codon_refseq_left + i + codon_refseq_right
            pep_altseq = "".join(AMINO_ACID_TABLE[codon] for codon in split_by_codon(codon_altseq))
            pep_altseq_set.add(pep_altseq)

        return sorted(pep_altseq_set)


class EnsemblRelease(Core):
    """Handles converting between position types and retrieving information on biological features
    based on an Ensembl release.
    """

    def __init__(
        self,
        species: str = DEFAULT_SPECIES,
        release: int = DEFAULT_ENSEMBL_RELEASE,
        cache_dir: str = "",
        canonical_transcript: str = "",
        contig_alias: str = "",
        exon_alias: str = "",
        gene_alias: str = "",
        protein_alias: str = "",
        transcript_alias: str = "",
    ):
        self.ensembl_cache = EnsemblCache(species, release, cache_dir=cache_dir)
        self.cache_dir = self.ensembl_cache.release_cache_dir
        self.reference = self.ensembl_cache.reference
        self.release = self.ensembl_cache.release
        self.species = self.ensembl_cache.species
        self.df = self.ensembl_cache.load_df()
        self.cds = []
        self.dna = [self.ensembl_cache.load_dna_fasta()]
        self.peptide = [self.ensembl_cache.load_pep_fasta()]
        self.rna = [self.ensembl_cache.load_cdna_fasta(), self.ensembl_cache.load_ncrna_fasta()]
        self.canonical_transcript = txt_to_list(canonical_transcript)
        self.contig_alias = tsv_to_dict(contig_alias)
        self.exon_alias = tsv_to_dict(exon_alias)
        self.gene_alias = tsv_to_dict(gene_alias)
        self.protein_alias = tsv_to_dict(protein_alias)
        self.transcript_alias = tsv_to_dict(transcript_alias)

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}(species={self.species}, release={self.release})"

    def install(
        self,
        clean: bool = True,
        recache: bool = False,
        redownload: bool = False,
        restrict_genes: List[str] = [],
    ):
        """Download missing data, process, and cache."""
        self.ensembl_cache.install(
            clean=clean, recache=recache, redownload=redownload, restrict_genes=restrict_genes
        )

    def cds_sequence(
        self, transcript_id: str, start: Optional[int] = None, end: Optional[int] = None
    ) -> str:
        """Return the nucleotide sequence at the given CDS coordinates."""
        offset = self.cds_offset(transcript_id)
        cds_start = start + offset if start is not None else offset
        cds_end = end + offset if end is not None else offset

        return self._get_sequence(self.rna, transcript_id, start=cds_start, end=cds_end)


def merge_positions(
    start_positions: List[PositionType], end_positions: List[PositionType], key: str
) -> List:
    """Merge two list of position objects into one list by the given key (e.g. 'transcript_id')"""
    result = set()

    for start_pos, end_pos in product(start_positions, end_positions):
        if start_pos.start > end_pos.start:
            start_pos, end_pos = end_pos, start_pos

        start_key = getattr(start_pos, key)
        end_key = getattr(end_pos, key)
        if start_key == end_key:
            assert start_pos.__class__ == end_pos.__class__
            kwargs = start_pos.asdict()
            kwargs.update(
                {
                    "start": start_pos.start,
                    "start_offset": start_pos.start_offset,
                    "end": end_pos.end,
                    "end_offset": end_pos.end_offset,
                }
            )
            new_position = start_pos.__class__(**kwargs)
            result.add(new_position)

    return sorted(result)
