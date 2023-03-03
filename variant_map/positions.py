from __future__ import annotations

from dataclasses import dataclass, fields
from itertools import product
from typing import TYPE_CHECKING, Any, Dict, List, TypeVar, Union

if TYPE_CHECKING:
    from .core import Core

from .constants import DELETION, DELINS, DUPLICATION, FRAMESHIFT, FUSION, INSERTION, SUBTITUTION
from .utils import (
    collapse_seq_change,
    expand_nt,
    format_hgvs_position,
    is_deletion,
    is_delins,
    is_duplication,
    is_frameshift,
    is_insertion,
    is_substitution,
    reverse_translate,
)


# -------------------------------------------------------------------------------------------------
# Position classes
# -------------------------------------------------------------------------------------------------
@dataclass(eq=True, frozen=True)
class _Position:
    """Base class for positions."""

    _data: Core
    contig_id: str
    start: int
    start_offset: int
    end: int
    end_offset: int
    strand: str

    @classmethod
    def copy_from(cls, position: _Position, **kwargs):
        """Initialize a new position object by copying values from another position object.

        Args:
            position (_Position): Position object to copy attributes from
            **kwargs: Keyword arguments with keys that match the attribute names of this position
                class will override attributes from `position`

        Returns:
            a new position object of the same class as the class that calls this method
        """
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

    def __lt__(self, other: _Position) -> bool:
        return str(self) < str(other)

    def __str__(self) -> str:
        raise NotImplementedError()  # Defined by inheriting class

    @property
    def is_cdna(self) -> bool:
        """Check if this fusion is between cDNA.

        Returns:
            bool: True if the fusion is between cDNA else False
        """
        return False

    @property
    def is_dna(self) -> bool:
        """Check if this position is on DNA.

        Returns:
            bool: True if the position is on DNA else False
        """
        return False

    @property
    def is_exon(self) -> bool:
        """Check if this fusion is between exon.

        Returns:
            bool: True if the fusion is between exon else False
        """
        return False

    @property
    def is_protein(self) -> bool:
        """Check if this fusion is between protein.

        Returns:
            bool: True if the fusion is between protein else False
        """
        return False

    @property
    def is_rna(self) -> bool:
        """Check if this fusion is between RNA.

        Returns:
            bool: True if the fusion is between RNA else False
        """
        return False

    @property
    def on_negative_strand(self) -> bool:
        """Check if this position originates from the negative strand of the genome.

        Returns:
            bool: True if the position originates from the negative strand of the genome else False
        """
        return self.strand == "-"

    @property
    def on_positive_strand(self) -> bool:
        """Check if this position originates from the positive strand of the genome.

        Returns:
            bool: True if the position originates from the positive strand of the genome else False
        """
        return self.strand == "+"

    def asdict(self) -> Dict[str, Any]:
        """Convert this position's attributes to a dictionary.

        Returns:
            Dict[str, Any]: Dictionary of attribute names and corresponding values
        """
        return {f.name: self[f.name] for f in fields(self)}

    def sequence(self) -> str:
        """Return the reference sequence from this position's start to it's end, inclusive.

        Returns:
            str: The reference sequence
        """
        raise NotImplementedError()  # Defined by inheriting class


@dataclass(eq=True, frozen=True)
class _CdnaPosition(_Position):
    """Base class for cDNA positions."""

    gene_id: str
    gene_name: str
    transcript_id: str
    transcript_name: str
    protein_id: str

    def __str__(self) -> str:
        start = format_hgvs_position(self.start, self.start_offset)
        end = format_hgvs_position(self.end, self.end_offset)
        if start == end:
            return f"{self.transcript_id}:c.{start}"
        else:
            return f"{self.transcript_id}:c.{start}_{end}"

    @property
    def is_cdna(self) -> bool:
        """Check if this position is on a cDNA.

        Returns:
            bool: True if the position is on a cDNA else False
        """
        return True

    def sequence(self) -> str:
        """Return the reference sequence from this position's start to it's end, inclusive.

        Returns:
            str: The reference sequence
        """
        if self.start_offset or self.end_offset:
            raise ValueError(
                f"Unable to get cDNA sequence of offset position ({self.start_offset}, {self.end_offset})"
            )

        return self._data.cds_sequence(self.transcript_id, self.start, self.end)


@dataclass(eq=True, frozen=True)
class _DnaPosition(_Position):
    """Base class for DNA positions."""

    def __str__(self) -> str:
        start = format_hgvs_position(self.start, self.start_offset)
        end = format_hgvs_position(self.end, self.end_offset)
        if start == end:
            return f"{self.contig_id}:g.{start}"
        else:
            return f"{self.contig_id}:g.{start}_{end}"

    @property
    def is_dna(self) -> bool:
        """Check if this position is on DNA.

        Returns:
            bool: True if the position is on DNA else False
        """
        return True

    def sequence(self) -> str:
        """Return the reference sequence from this position's start to it's end, inclusive.

        Returns:
            str: The reference sequence
        """
        if self.start_offset or self.end_offset:
            raise ValueError(
                f"Unable to get DNA sequence of offset position ({self.start_offset}, {self.end_offset})"
            )

        return self._data.dna_sequence(self.contig_id, self.start, self.end, self.strand)


@dataclass(eq=True, frozen=True)
class _ExonPosition(_Position):
    """Base class for exon positions."""

    gene_id: str
    gene_name: str
    transcript_id: str
    transcript_name: str
    # NOTE: If the start and end are different, the exon ID of the start is listed as the `exon_id`
    exon_id: str

    def __str__(self) -> str:
        start = format_hgvs_position(self.start, self.start_offset)
        end = format_hgvs_position(self.end, self.end_offset)
        if start == end:
            return f"{self.transcript_id}:e.{start}"
        else:
            return f"{self.transcript_id}:e.{start}_{end}"

    @property
    def is_exon(self) -> bool:
        """Check if this position is on an exon.

        Returns:
            bool: True if the position is on an exon else False
        """
        return True

    def sequence(self) -> str:
        """Return the reference sequence from this position's start to it's end, inclusive.

        Returns:
            str: The reference sequence
        """
        raise NotImplementedError()  # TODO


@dataclass(eq=True, frozen=True)
class _ProteinPosition(_Position):
    """Base class for protein positions."""

    gene_id: str
    gene_name: str
    transcript_id: str
    transcript_name: str
    protein_id: str

    def __str__(self) -> str:
        start = format_hgvs_position(self.start, self.start_offset)
        end = format_hgvs_position(self.end, self.end_offset)
        if start == end:
            return f"{self.protein_id}:p.{start}"
        else:
            return f"{self.protein_id}:p.{start}_{end}"

    @property
    def is_protein(self) -> bool:
        """Check if this position is on a protein.

        Returns:
            bool: True if the position is on a protein else False
        """
        return True

    def sequence(self) -> str:
        """Return the reference sequence from this position's start to it's end, inclusive.

        Returns:
            str: The reference sequence
        """
        if self.start_offset or self.end_offset:
            raise ValueError(
                f"Unable to get peptide sequence of offset position ({self.start_offset}, {self.end_offset})"
            )

        return self._data.protein_sequence(self.protein_id, self.start, self.end)


@dataclass(eq=True, frozen=True)
class _RnaPosition(_Position):
    """Base class for RNA positions."""

    gene_id: str
    gene_name: str
    transcript_id: str
    transcript_name: str

    def __str__(self) -> str:
        start = format_hgvs_position(self.start, self.start_offset)
        end = format_hgvs_position(self.end, self.end_offset)
        if start == end:
            return f"{self.transcript_id}:r.{start}"
        else:
            return f"{self.transcript_id}:r.{start}_{end}"

    @property
    def is_rna(self) -> bool:
        """Check if this position is on an RNA.

        Returns:
            bool: True if the position is on an RNA else False
        """
        return True

    def sequence(self) -> str:
        """Return the reference sequence from this position's start to it's end, inclusive.

        Returns:
            str: The reference sequence
        """
        if self.start_offset or self.end_offset:
            raise ValueError(
                f"Unable to get RNA sequence of offset position ({self.start_offset}, {self.end_offset})"
            )

        return self._data.rna_sequence(self.transcript_id, self.start, self.end)


# -------------------------------------------------------------------------------------------------
# MappablePosition classes
# -------------------------------------------------------------------------------------------------
@dataclass(eq=True, frozen=True)
class _MappablePosition:
    """Base class for mappable positions."""

    def to_cdna(self, canonical: bool = False) -> List[CdnaPosition]:
        """Map this position to one more cDNA positions.

        Args:
            canonical (bool, optional): Only map to the canonical transcript. Defaults to False.

        Returns:
            List[CdnaPosition]: One or more cDNA positions.
        """
        raise NotImplementedError()  # Defined by inheriting class

    def to_dna(self, canonical: bool = False) -> List[DnaPosition]:
        """Map this position to one more DNA positions.

        Args:
            canonical (bool, optional): Only map to the canonical transcript. Defaults to False.

        Returns:
            List[DnaPosition]: One or more DNA positions.
        """
        raise NotImplementedError()  # Defined by inheriting class

    def to_exon(self, canonical: bool = False) -> List[ExonPosition]:
        """Map this position to one more exon positions.

        Args:
            canonical (bool, optional): Only map to the canonical transcript. Defaults to False.

        Returns:
            List[ExonPosition]: One or more exon positions.
        """
        raise NotImplementedError()  # Defined by inheriting class

    def to_protein(self, canonical: bool = False) -> List[ProteinPosition]:
        """Map this position to one more protein positions.

        Args:
            canonical (bool, optional): Only map to the canonical transcript. Defaults to False.

        Returns:
            List[ProteinPosition]: One or more protein positions.
        """
        raise NotImplementedError()  # Defined by inheriting class

    def to_rna(self, canonical: bool = False) -> List[RnaPosition]:
        """Map this position to one more RNA positions.

        Args:
            canonical (bool, optional): Only map to the canonical transcript. Defaults to False.

        Returns:
            List[RnaPosition]: One or more RNA positions.
        """
        raise NotImplementedError()  # Defined by inheriting class


@dataclass(eq=True, frozen=True)
class CdnaPosition(_CdnaPosition, _MappablePosition):
    """Stores information on a cDNA position and maps to other position types."""

    def to_cdna(self, canonical: bool = False) -> List[CdnaPosition]:
        """Map this position to one more cDNA positions.

        Args:
            canonical (bool, optional): Only map to the canonical transcript. Defaults to False.

        Returns:
            List[CdnaPosition]: One or more cDNA positions.
        """
        return self._data._cdna_to_cdna(
            [self.transcript_id],
            self.start,
            self.start_offset,
            self.end,
            self.end_offset,
            [self.strand],
            canonical,
        )

    def to_dna(self, canonical: bool = False) -> List[DnaPosition]:
        """Map this position to one more DNA positions.

        Args:
            canonical (bool, optional): Only map to the canonical transcript. Defaults to False.

        Returns:
            List[DnaPosition]: One or more DNA positions.
        """
        return self._data._cdna_to_dna(
            [self.transcript_id],
            self.start,
            self.start_offset,
            self.end,
            self.end_offset,
            [self.strand],
            canonical,
        )

    def to_exon(self, canonical: bool = False) -> List[ExonPosition]:
        """Map this position to one more exon positions.

        Args:
            canonical (bool, optional): Only map to the canonical transcript. Defaults to False.

        Returns:
            List[ExonPosition]: One or more exon positions.
        """
        return self._data._cdna_to_exon(
            [self.transcript_id],
            self.start,
            self.start_offset,
            self.end,
            self.end_offset,
            [self.strand],
            canonical,
        )

    def to_protein(self, canonical: bool = False) -> List[ProteinPosition]:
        """Map this position to one more protein positions.

        Args:
            canonical (bool, optional): Only map to the canonical transcript. Defaults to False.

        Returns:
            List[ProteinPosition]: One or more protein positions.
        """
        return self._data._cdna_to_protein(
            [self.transcript_id],
            self.start,
            self.start_offset,
            self.end,
            self.end_offset,
            [self.strand],
            canonical,
        )

    def to_rna(self, canonical: bool = False) -> List[RnaPosition]:
        """Map this position to one more RNA positions.

        Args:
            canonical (bool, optional): Only map to the canonical transcript. Defaults to False.

        Returns:
            List[RnaPosition]: One or more RNA positions.
        """
        return self._data._cdna_to_rna(
            [self.transcript_id],
            self.start,
            self.start_offset,
            self.end,
            self.end_offset,
            [self.strand],
            canonical,
        )


@dataclass(eq=True, frozen=True)
class DnaPosition(_DnaPosition, _MappablePosition):
    """Stores information on a DNA position and maps to other position types."""

    def to_cdna(self, canonical: bool = False) -> List[CdnaPosition]:
        """Map this position to one more cDNA positions.

        Args:
            canonical (bool, optional): Only map to the canonical transcript. Defaults to False.

        Returns:
            List[CdnaPosition]: One or more cDNA positions.
        """
        return self._data._dna_to_cdna(
            [self.contig_id],
            self.start,
            self.start_offset,
            self.end,
            self.end_offset,
            [self.strand],
            canonical,
        )

    def to_dna(self, canonical: bool = False) -> List[DnaPosition]:
        """Map this position to one more DNA positions.

        Args:
            canonical (bool, optional): Only map to the canonical transcript. Defaults to False.

        Returns:
            List[DnaPosition]: One or more DNA positions.
        """
        return self._data._dna_to_dna(
            [self.contig_id],
            self.start,
            self.start_offset,
            self.end,
            self.end_offset,
            [self.strand],
            canonical,
        )

    def to_exon(self, canonical: bool = False) -> List[ExonPosition]:
        """Map this position to one more exon positions.

        Args:
            canonical (bool, optional): Only map to the canonical transcript. Defaults to False.

        Returns:
            List[ExonPosition]: One or more exon positions.
        """
        return self._data._dna_to_exon(
            [self.contig_id],
            self.start,
            self.start_offset,
            self.end,
            self.end_offset,
            [self.strand],
            canonical,
        )

    def to_protein(self, canonical: bool = False) -> List[ProteinPosition]:
        """Map this position to one more protein positions.

        Args:
            canonical (bool, optional): Only map to the canonical transcript. Defaults to False.

        Returns:
            List[ProteinPosition]: One or more protein positions.
        """
        return self._data._dna_to_protein(
            [self.contig_id],
            self.start,
            self.start_offset,
            self.end,
            self.end_offset,
            [self.strand],
            canonical,
        )

    def to_rna(self, canonical: bool = False) -> List[RnaPosition]:
        """Map this position to one more RNA positions.

        Args:
            canonical (bool, optional): Only map to the canonical transcript. Defaults to False.

        Returns:
            List[RnaPosition]: One or more RNA positions.
        """
        return self._data._dna_to_rna(
            [self.contig_id],
            self.start,
            self.start_offset,
            self.end,
            self.end_offset,
            [self.strand],
            canonical,
        )


@dataclass(eq=True, frozen=True)
class ExonPosition(_ExonPosition, _MappablePosition):
    """Stores information on an exon position and maps to other position types."""

    def to_cdna(self, canonical: bool = False) -> List[CdnaPosition]:
        """Map this position to one more cDNA positions.

        Args:
            canonical (bool, optional): Only map to the canonical transcript. Defaults to False.

        Returns:
            List[CdnaPosition]: One or more cDNA positions.
        """
        return self._data._exon_to_cdna(
            [self.transcript_id],
            self.start,
            self.start_offset,
            self.end,
            self.end_offset,
            [self.strand],
            canonical,
        )

    def to_dna(self, canonical: bool = False) -> List[DnaPosition]:
        """Map this position to one more DNA positions.

        Args:
            canonical (bool, optional): Only map to the canonical transcript. Defaults to False.

        Returns:
            List[DnaPosition]: One or more DNA positions.
        """
        return self._data._exon_to_dna(
            [self.transcript_id],
            self.start,
            self.start_offset,
            self.end,
            self.end_offset,
            [self.strand],
            canonical,
        )

    def to_exon(self, canonical: bool = False) -> List[ExonPosition]:
        """Map this position to one more exon positions.

        Args:
            canonical (bool, optional): Only map to the canonical transcript. Defaults to False.

        Returns:
            List[ExonPosition]: One or more exon positions.
        """
        return self._data._exon_to_exon(
            [self.transcript_id],
            self.start,
            self.start_offset,
            self.end,
            self.end_offset,
            [self.strand],
            canonical,
        )

    def to_protein(self, canonical: bool = False) -> List[ProteinPosition]:
        """Map this position to one more protein positions.

        Args:
            canonical (bool, optional): Only map to the canonical transcript. Defaults to False.

        Returns:
            List[ProteinPosition]: One or more protein positions.
        """
        return self._data._exon_to_protein(
            [self.transcript_id],
            self.start,
            self.start_offset,
            self.end,
            self.end_offset,
            [self.strand],
            canonical,
        )

    def to_rna(self, canonical: bool = False) -> List[RnaPosition]:
        """Map this position to one more RNA positions.

        Args:
            canonical (bool, optional): Only map to the canonical transcript. Defaults to False.

        Returns:
            List[RnaPosition]: One or more RNA positions.
        """
        return self._data._exon_to_rna(
            [self.transcript_id],
            self.start,
            self.start_offset,
            self.end,
            self.end_offset,
            [self.strand],
            canonical,
        )


@dataclass(eq=True, frozen=True)
class ProteinPosition(_ProteinPosition, _MappablePosition):
    """Stores information on a protein position and maps to other position types."""

    def to_cdna(self, canonical: bool = False) -> List[CdnaPosition]:
        """Map this position to one more cDNA positions.

        Args:
            canonical (bool, optional): Only map to the canonical transcript. Defaults to False.

        Returns:
            List[CdnaPosition]: One or more cDNA positions.
        """
        return self._data._protein_to_cdna(
            [self.transcript_id],
            self.start,
            self.start_offset,
            self.end,
            self.end_offset,
            [self.strand],
            canonical,
        )

    def to_dna(self, canonical: bool = False) -> List[DnaPosition]:
        """Map this position to one more DNA positions.

        Args:
            canonical (bool, optional): Only map to the canonical transcript. Defaults to False.

        Returns:
            List[DnaPosition]: One or more DNA positions.
        """
        return self._data._protein_to_dna(
            [self.transcript_id],
            self.start,
            self.start_offset,
            self.end,
            self.end_offset,
            [self.strand],
            canonical,
        )

    def to_exon(self, canonical: bool = False) -> List[ExonPosition]:
        """Map this position to one more exon positions.

        Args:
            canonical (bool, optional): Only map to the canonical transcript. Defaults to False.

        Returns:
            List[ExonPosition]: One or more exon positions.
        """
        return self._data._protein_to_exon(
            [self.transcript_id],
            self.start,
            self.start_offset,
            self.end,
            self.end_offset,
            [self.strand],
            canonical,
        )

    def to_protein(self, canonical: bool = False) -> List[ProteinPosition]:
        """Map this position to one more protein positions.

        Args:
            canonical (bool, optional): Only map to the canonical transcript. Defaults to False.

        Returns:
            List[ProteinPosition]: One or more protein positions.
        """
        return self._data._protein_to_protein(
            [self.transcript_id],
            self.start,
            self.start_offset,
            self.end,
            self.end_offset,
            [self.strand],
            canonical,
        )

    def to_rna(self, canonical: bool = False) -> List[RnaPosition]:
        """Map this position to one more RNA positions.

        Args:
            canonical (bool, optional): Only map to the canonical transcript. Defaults to False.

        Returns:
            List[RnaPosition]: One or more RNA positions.
        """
        return self._data._protein_to_rna(
            [self.transcript_id],
            self.start,
            self.start_offset,
            self.end,
            self.end_offset,
            [self.strand],
            canonical,
        )


@dataclass(eq=True, frozen=True)
class RnaPosition(_RnaPosition, _MappablePosition):
    """Stores information on an RNA position and maps to other position types."""

    def to_cdna(self, canonical: bool = False) -> List[CdnaPosition]:
        """Map this position to one more cDNA positions.

        Args:
            canonical (bool, optional): Only map to the canonical transcript. Defaults to False.

        Returns:
            List[CdnaPosition]: One or more cDNA positions.
        """
        return self._data._rna_to_cdna(
            [self.transcript_id],
            self.start,
            self.start_offset,
            self.end,
            self.end_offset,
            [self.strand],
            canonical,
        )

    def to_dna(self, canonical: bool = False) -> List[DnaPosition]:
        """Map this position to one more DNA positions.

        Args:
            canonical (bool, optional): Only map to the canonical transcript. Defaults to False.

        Returns:
            List[DnaPosition]: One or more DNA positions.
        """
        return self._data._rna_to_dna(
            [self.transcript_id],
            self.start,
            self.start_offset,
            self.end,
            self.end_offset,
            [self.strand],
            canonical,
        )

    def to_exon(self, canonical: bool = False) -> List[ExonPosition]:
        """Map this position to one more exon positions.

        Args:
            canonical (bool, optional): Only map to the canonical transcript. Defaults to False.

        Returns:
            List[ExonPosition]: One or more exon positions.
        """
        return self._data._rna_to_exon(
            [self.transcript_id],
            self.start,
            self.start_offset,
            self.end,
            self.end_offset,
            [self.strand],
            canonical,
        )

    def to_protein(self, canonical: bool = False) -> List[ProteinPosition]:
        """Map this position to one more protein positions.

        Args:
            canonical (bool, optional): Only map to the canonical transcript. Defaults to False.

        Returns:
            List[ProteinPosition]: One or more protein positions.
        """
        return self._data._rna_to_protein(
            [self.transcript_id],
            self.start,
            self.start_offset,
            self.end,
            self.end_offset,
            [self.strand],
            canonical,
        )

    def to_rna(self, canonical: bool = False) -> List[RnaPosition]:
        """Map this position to one more RNA positions.

        Args:
            canonical (bool, optional): Only map to the canonical transcript. Defaults to False.

        Returns:
            List[RnaPosition]: One or more RNA positions.
        """
        return self._data._rna_to_rna(
            [self.transcript_id],
            self.start,
            self.start_offset,
            self.end,
            self.end_offset,
            [self.strand],
            canonical,
        )


# -------------------------------------------------------------------------------------------------
# Variant classes
# -------------------------------------------------------------------------------------------------
@dataclass(eq=True, frozen=True)
class _Variant:
    """Base class for variants."""

    @property
    def is_deletion(self) -> bool:
        """Check if this variant is a deletion.

        Returns:
            bool: True if the position is a deletion else False
        """
        return self.type == DELETION

    @property
    def is_delins(self) -> bool:
        """Check if this variant is a delins (indel).

        Returns:
            bool: True if the position is a delins else False
        """
        return self.type == DELINS

    @property
    def is_duplication(self) -> bool:
        """Check if this variant is a duplication.

        Returns:
            bool: True if the position is a duplication else False
        """
        return self.type == DUPLICATION

    @property
    def is_frameshift(self) -> bool:
        """Check if this variant is a frameshift.

        Returns:
            bool: True if the position is a frameshift else False
        """
        return self.type == FRAMESHIFT

    @property
    def is_fusion(self) -> bool:
        """Check if this variant is a fusion.

        Returns:
            bool: True if the position is a fusion else False
        """
        return self.type == FUSION

    @property
    def is_insertion(self) -> bool:
        """Check if this variant is an insertion.

        Returns:
            bool: True if the position is an insertion else False
        """
        return self.type == INSERTION

    @property
    def is_substitution(self) -> bool:
        """Check if this variant is a substitution.

        Returns:
            bool: True if the position is a substitution else False
        """
        return self.type == SUBTITUTION

    @property
    def type(self) -> str:
        """Get the variant type.

        Returns:
            str: Type of variant
        """
        raise NotImplementedError()  # Defined by inheriting class

    def to_cdna(self) -> List:
        """Map this variant to one more cDNA variants.

        Args:
            canonical (bool, optional): Only map to the canonical transcript. Defaults to False.

        Returns:
            List: One or more cDNA variants.
        """
        raise NotImplementedError()  # Defined by inheriting class

    def to_dna(self) -> List:
        """Map this variant to one more DNA variants.

        Args:
            canonical (bool, optional): Only map to the canonical transcript. Defaults to False.

        Returns:
            List: One or more DNA variants.
        """
        raise NotImplementedError()  # Defined by inheriting class

    def to_exon(self) -> List:
        """Map this variant to one more exon variants.

        Args:
            canonical (bool, optional): Only map to the canonical transcript. Defaults to False.

        Returns:
            List: One or more exon variants.
        """
        raise NotImplementedError()  # Defined by inheriting class

    def to_protein(self) -> List:
        """Map this variant to one more protein variants.

        Args:
            canonical (bool, optional): Only map to the canonical transcript. Defaults to False.

        Returns:
            List: One or more protein variants.
        """
        raise NotImplementedError()  # Defined by inheriting class

    def to_rna(self) -> List:
        """Map this variant to one more RNA variants.

        Args:
            canonical (bool, optional): Only map to the canonical transcript. Defaults to False.

        Returns:
            List: One or more RNA variants.
        """
        raise NotImplementedError()  # Defined by inheriting class


# -------------------------------------------------------------------------------------------------
# SmallVariant classes
# -------------------------------------------------------------------------------------------------
@dataclass(eq=True, frozen=True)
class _SmallVariant(_Variant):
    """Base class for small variants."""

    refseq: str
    altseq: str


@dataclass(eq=True, frozen=True)
class _CdnaSmallVariant(_CdnaPosition, _SmallVariant):
    """Base class for cDNA small variants."""

    @classmethod
    def from_cdna(cls, cdna: CdnaPosition, refseq: str, altseq: str) -> List[_CdnaSmallVariant]:
        """Convert a cDNA position plus ref/alt nucleotides into a cDNA variant.

        Args:
            cdna (CdnaPosition): _Variant position
            refseq (str): reference allele
            altseq (str): alternate allele

        Raises:
            NotImplementedError: An unsupported combination of reference/alternate alleles was given
            ValueError: The given reference allele does not match the annotated reference allele

        Returns:
            List[_CdnaSmallVariant]: One or more cDNA variants
        """
        variant_list = []

        ref_annotated = cdna.sequence()
        for ref, alt in product(expand_nt(refseq), expand_nt(altseq)):
            # Assert that the given ref matches the annotated one
            if ref != ref_annotated:
                raise ValueError(
                    f"Given ref allele '{ref}' does not match annotated ref allele '{ref_annotated}'"
                )

            # Trim bases that are unchanged between the ref and alt alleles
            new_ref, new_alt, new_start, new_end = collapse_seq_change(ref, alt)
            start = cdna.start + new_start
            end = cdna.end - new_end

            # Determine the type of variant
            if is_deletion(new_ref, new_alt):
                variant = CdnaDeletion.copy_from(
                    cdna, start=start, end=end, refseq=new_ref, altseq=new_alt
                )
            elif is_delins(new_ref, new_alt):
                variant = CdnaDelins.copy_from(
                    cdna, start=start, end=end, refseq=new_ref, altseq=new_alt
                )
            elif is_duplication(new_ref, new_alt):
                variant = CdnaDuplication.copy_from(
                    cdna, start=start, end=end, refseq=new_ref, altseq=new_alt
                )
            elif is_insertion(new_ref, new_alt):
                variant = CdnaInsertion.copy_from(
                    cdna, start=start, end=end, refseq=new_ref, altseq=new_alt
                )
            elif is_substitution(new_ref, new_alt):
                variant = CdnaSubstitution.copy_from(
                    cdna, start=start, end=end, refseq=new_ref, altseq=new_alt
                )
            else:
                raise NotImplementedError(
                    f"Unrecognized sequence change '{new_ref}/{new_alt}' at {cdna}"
                )

            variant_list.append(variant)

        return variant_list

    @classmethod
    def from_protein(cls, cdna: CdnaPosition, refseq: str, altseq: str) -> List[_CdnaSmallVariant]:
        """Convert a cDNA position plus ref/alt amino acids into a cDNA variant.

        Args:
            cdna (CdnaPosition): cDNA position
            refseq (str): Reference allele
            altseq (str): Alternate allele

        Raises:
            NotImplementedError: An unsupported combination of reference/alternate alleles was given

        Returns:
            List[_CdnaSmallVariant]: One or more cDNA variants
        """
        variant_list = []

        ref_annotated = cdna.sequence()
        for ref, alt in product(reverse_translate(refseq), reverse_translate(altseq)):
            # Assert that the given ref matches the annotated one
            if ref != ref_annotated:
                continue

            # For insertions, check that the sequence flanking the inserted sequence matches the ref
            if is_insertion(refseq, altseq) and not is_insertion(ref, alt):
                continue

            # Trim bases that are unchanged between the ref and alt alleles
            new_ref, new_alt, new_start, new_end = collapse_seq_change(ref, alt)
            start = cdna.start + new_start
            end = cdna.end - new_end

            # Determine the type of variant
            if is_deletion(new_ref, new_alt):
                variant = CdnaDeletion.copy_from(
                    cdna, start=start, end=end, refseq=new_ref, altseq=new_alt
                )
            elif is_delins(new_ref, new_alt):
                variant = CdnaDelins.copy_from(
                    cdna, start=start, end=end, refseq=new_ref, altseq=new_alt
                )
            elif is_duplication(new_ref, new_alt):
                variant = CdnaDuplication.copy_from(
                    cdna, start=start, end=end, refseq=new_ref, altseq=new_alt
                )
            elif is_insertion(new_ref, new_alt):
                variant = CdnaInsertion.copy_from(
                    cdna, start=start, end=end, refseq=new_ref, altseq=new_alt
                )
            elif is_substitution(new_ref, new_alt):
                variant = CdnaSubstitution.copy_from(
                    cdna, start=start, end=end, refseq=new_ref, altseq=new_alt
                )
            else:
                raise NotImplementedError(
                    f"Unrecognized sequence change '{new_ref}/{new_alt}' at {cdna}"
                )

            variant_list.append(variant)

        return variant_list

    def to_cdna(self, canonical: bool = False) -> List[_CdnaSmallVariant]:
        """Map this variants to one more cDNA variants.

        Args:
            canonical (bool, optional): Only map to the canonical transcript. Defaults to False.

        Returns:
            List[_CdnaSmallVariant]: One or more cDNA variants.
        """
        return self._data._cdna_to_cdna_variant(
            [self.transcript_id],
            self.start,
            self.start_offset,
            self.end,
            self.end_offset,
            [self.strand],
            self.refseq,
            self.altseq,
            canonical,
        )

    def to_dna(self, canonical: bool = False) -> List[_DnaSmallVariant]:
        """Map this variants to one more DNA variants.

        Args:
            canonical (bool, optional): Only map to the canonical transcript. Defaults to False.

        Returns:
            List[_DnaSmallVariant]: One or more DNA variants.
        """
        return self._data._cdna_to_dna_variant(
            [self.transcript_id],
            self.start,
            self.start_offset,
            self.end,
            self.end_offset,
            [self.strand],
            self.refseq,
            self.altseq,
            canonical,
        )

    def to_exon(self, canonical: bool = False) -> List[_ExonSmallVariant]:
        """Map this variants to one more exon variants.

        Args:
            canonical (bool, optional): Only map to the canonical transcript. Defaults to False.

        Returns:
            List[_ExonSmallVariant]: One or more exon variants.
        """
        raise NotImplementedError()  # TODO

    def to_protein(self, canonical: bool = False) -> List[_ProteinSmallVariant]:
        """Map this variants to one more protein variants.

        Args:
            canonical (bool, optional): Only map to the canonical transcript. Defaults to False.

        Returns:
            List[_ProteinSmallVariant]: One or more protein variants.
        """
        return self._data._cdna_to_protein_variant(
            [self.transcript_id],
            self.start,
            self.start_offset,
            self.end,
            self.end_offset,
            [self.strand],
            self.refseq,
            self.altseq,
            canonical,
        )

    def to_rna(self, canonical: bool = False) -> List[_RnaSmallVariant]:
        """Map this variants to one more RNA variants.

        Args:
            canonical (bool, optional): Only map to the canonical transcript. Defaults to False.

        Returns:
            List[_RnaSmallVariant]: One or more RNA variants.
        """
        return self._data._cdna_to_rna_variant(
            [self.transcript_id],
            self.start,
            self.start_offset,
            self.end,
            self.end_offset,
            [self.strand],
            self.refseq,
            self.altseq,
            canonical,
        )


@dataclass(eq=True, frozen=True)
class _DnaSmallVariant(_DnaPosition, _SmallVariant):
    """Base class for DNA small variants."""

    @classmethod
    def from_dna(cls, dna: DnaPosition, refseq: str, altseq: str) -> List[_DnaSmallVariant]:
        """Convert a DNA position plus ref/alt nucleotides into a DNA variant.

        Args:
            dna (DnaPosition): DNA position
            refseq (str): Reference allele
            altseq (str): Alternate allele

        Raises:
            NotImplementedError: An unsupported combination of reference/alternate alleles was given

        Returns:
            List[_DnaSmallVariant]: One or more DNA variants
        """
        variant_list = []

        # TODO: DNA sequence get is slow
        # ref_annotated = dna.sequence()
        for ref, alt in product(expand_nt(refseq), expand_nt(altseq)):
            # # Assert that the given ref matches the annotated one
            # if ref != ref_annotated:
            #     raise ValueError(
            #         f"Given ref allele '{ref}' does not match annotated ref allele '{ref_annotated}' for {dna}"
            #     )

            # Trim bases that are unchanged between the ref and alt alleles
            new_ref, new_alt, new_start, new_end = collapse_seq_change(ref, alt)
            start = dna.start + new_start
            end = dna.end - new_end

            # Determine the type of variant
            if is_deletion(new_ref, new_alt):
                variant = DnaDeletion.copy_from(
                    dna, start=start, end=end, refseq=new_ref, altseq=new_alt
                )
            elif is_delins(new_ref, new_alt):
                variant = DnaDelins.copy_from(
                    dna, start=start, end=end, refseq=new_ref, altseq=new_alt
                )
            elif is_duplication(new_ref, new_alt):
                variant = DnaDuplication.copy_from(
                    dna, start=start, end=end, refseq=new_ref, altseq=new_alt
                )
            elif is_insertion(new_ref, new_alt):
                variant = DnaInsertion.copy_from(
                    dna, start=start, end=end, refseq=new_ref, altseq=new_alt
                )
            elif is_substitution(new_ref, new_alt):
                variant = DnaSubstitution.copy_from(
                    dna, start=start, end=end, refseq=new_ref, altseq=new_alt
                )
            else:
                raise NotImplementedError(
                    f"Unrecognized sequence change '{new_ref}/{new_alt}' at {dna}"
                )

            variant_list.append(variant)

        return variant_list

    def to_cdna(self, canonical: bool = False) -> List[_CdnaSmallVariant]:
        """Map this variants to one more cDNA variants.

        Args:
            canonical (bool, optional): Only map to the canonical transcript. Defaults to False.

        Returns:
            List[_CdnaSmallVariant]: One or more cDNA variants.
        """
        return self._data._dna_to_cdna_variant(
            [self.contig_id],
            self.start,
            self.start_offset,
            self.end,
            self.end_offset,
            [self.strand],
            self.refseq,
            self.altseq,
            canonical,
        )

    def to_dna(self, canonical: bool = False) -> List[_DnaSmallVariant]:
        """Map this variants to one more DNA variants.

        Args:
            canonical (bool, optional): Only map to the canonical transcript. Defaults to False.

        Returns:
            List[_DnaSmallVariant]: One or more DNA variants.
        """
        return self._data._dna_to_dna_variant(
            [self.contig_id],
            self.start,
            self.start_offset,
            self.end,
            self.end_offset,
            [self.strand],
            self.refseq,
            self.altseq,
            canonical,
        )

    def to_exon(self, canonical: bool = False) -> List[_ExonSmallVariant]:
        """Map this variants to one more exon variants.

        Args:
            canonical (bool, optional): Only map to the canonical transcript. Defaults to False.

        Returns:
            List[_ExonSmallVariant]: One or more exon variants.
        """
        raise NotImplementedError()  # TODO

    def to_protein(self, canonical: bool = False) -> List[_ProteinSmallVariant]:
        """Map this variants to one more protein variants.

        Args:
            canonical (bool, optional): Only map to the canonical transcript. Defaults to False.

        Returns:
            List[_ProteinSmallVariant]: One or more protein variants.
        """
        return self._data._dna_to_protein_variant(
            [self.contig_id],
            self.start,
            self.start_offset,
            self.end,
            self.end_offset,
            [self.strand],
            self.refseq,
            self.altseq,
            canonical,
        )

    def to_rna(self, canonical: bool = False) -> List[_RnaSmallVariant]:
        """Map this variants to one more RNA variants.

        Args:
            canonical (bool, optional): Only map to the canonical transcript. Defaults to False.

        Returns:
            List[_RnaSmallVariant]: One or more RNA variants.
        """
        return self._data._dna_to_rna_variant(
            [self.contig_id],
            self.start,
            self.start_offset,
            self.end,
            self.end_offset,
            [self.strand],
            self.refseq,
            self.altseq,
            canonical,
        )


@dataclass(eq=True, frozen=True)
class _ExonSmallVariant(_ExonPosition, _SmallVariant):
    """Base class for exon small variants."""

    @classmethod
    def from_exon(cls, exon: ExonPosition, refseq: str, altseq: str) -> List[_ExonSmallVariant]:
        """Convert an exon position plus ref/alt nucleotides into an exon variant.

        Args:
            exon (ExonPosition): exon position
            refseq (str): Reference allele
            altseq (str): Alternate allele

        Raises:
            NotImplementedError: An unsupported combination of reference/alternate alleles was given

        Returns:
            List[_ExonSmallVariant]: One or more exon variants
        """
        raise NotImplementedError()  # TODO


@dataclass(eq=True, frozen=True)
class _ProteinSmallVariant(_ProteinPosition, _SmallVariant):
    """Base class for protein small variants."""

    @classmethod
    def from_cdna(
        cls, cdna: CdnaPosition, protein: ProteinPosition, refseq: str, altseq: str
    ) -> List[_ProteinSmallVariant]:
        """Convert a cDNA position plus ref/alt amino acids into a protein variant.

        Args:
            cdna (CdnaPosition): cDNA position
            refseq (str): Reference allele
            altseq (str): Alternate allele

        Raises:
            NotImplementedError: An unsupported combination of reference/alternate alleles was given

        Returns:
            List[_ProteinSmallVariant]: One or more protein variants
        """
        variant_list = []

        protein_refseq = protein.sequence()
        for cdna_ref, cdna_alt in product(expand_nt(refseq), expand_nt(altseq)):
            if is_frameshift(cdna_ref, cdna_alt):
                raise NotImplementedError()  # TODO

            for protein_alt in protein._data.translate_cds_variant(
                cdna.transcript_id, cdna.start, cdna.end, cdna_alt
            ):
                # Trim bases that are unchanged between the ref and alt alleles
                new_ref, new_alt, new_start, new_end = collapse_seq_change(
                    protein_refseq, protein_alt
                )
                start = protein.start + new_start
                end = protein.end - new_end

                # Determine the type of variant
                if is_deletion(new_ref, new_alt):
                    variant = ProteinDeletion.copy_from(
                        protein, start=start, end=end, refseq=new_ref, altseq=new_alt
                    )
                elif is_delins(new_ref, new_alt):
                    variant = ProteinDelins.copy_from(
                        protein, start=start, end=end, refseq=new_ref, altseq=new_alt
                    )
                elif is_duplication(new_ref, new_alt):
                    variant = ProteinDuplication.copy_from(
                        protein, start=start, end=end, refseq=new_ref, altseq=new_alt
                    )
                elif is_insertion(new_ref, new_alt):
                    variant = ProteinInsertion.copy_from(
                        protein, start=start, end=end, refseq=new_ref, altseq=new_alt
                    )
                elif is_substitution(new_ref, new_alt):
                    variant = ProteinSubstitution.copy_from(
                        protein, start=start, end=end, refseq=new_ref, altseq=new_alt
                    )
                else:
                    raise NotImplementedError(
                        f"Unrecognized sequence change '{new_ref}/{new_alt}' at {protein}"
                    )

                variant_list.append(variant)

        return variant_list

    def to_cdna(self, canonical: bool = False) -> List[_CdnaSmallVariant]:
        """Map this variants to one more cDNA variants.

        Args:
            canonical (bool, optional): Only map to the canonical transcript. Defaults to False.

        Returns:
            List[_CdnaSmallVariant]: One or more cDNA variants.
        """
        return self._data._protein_to_cdna_variant(
            [self.transcript_id],
            self.start,
            self.start_offset,
            self.end,
            self.end_offset,
            [self.strand],
            self.refseq,
            self.altseq,
            canonical,
        )

    def to_dna(self, canonical: bool = False) -> List[_DnaSmallVariant]:
        """Map this variants to one more DNA variants.

        Args:
            canonical (bool, optional): Only map to the canonical transcript. Defaults to False.

        Returns:
            List[_DnaSmallVariant]: One or more DNA variants.
        """
        return self._data._protein_to_dna_variant(
            [self.transcript_id],
            self.start,
            self.start_offset,
            self.end,
            self.end_offset,
            [self.strand],
            self.refseq,
            self.altseq,
            canonical,
        )

    def to_exon(self, canonical: bool = False) -> List[_ExonSmallVariant]:
        """Map this variants to one more exon variants.

        Args:
            canonical (bool, optional): Only map to the canonical transcript. Defaults to False.

        Returns:
            List[_ExonSmallVariant]: One or more exon variants.
        """
        raise NotImplementedError()  # TODO

    def to_protein(self, canonical: bool = False) -> List[_ProteinSmallVariant]:
        """Map this variants to one more protein variants.

        Args:
            canonical (bool, optional): Only map to the canonical transcript. Defaults to False.

        Returns:
            List[_ProteinSmallVariant]: One or more protein variants.
        """
        return self._data._protein_to_protein_variant(
            [self.transcript_id],
            self.start,
            self.start_offset,
            self.end,
            self.end_offset,
            [self.strand],
            self.refseq,
            self.altseq,
            canonical,
        )

    def to_rna(self, canonical: bool = False) -> List[_RnaSmallVariant]:
        """Map this variants to one more RNA variants.

        Args:
            canonical (bool, optional): Only map to the canonical transcript. Defaults to False.

        Returns:
            List[_RnaSmallVariant]: One or more RNA variants.
        """
        return self._data._protein_to_rna_variant(
            [self.transcript_id],
            self.start,
            self.start_offset,
            self.end,
            self.end_offset,
            [self.strand],
            self.refseq,
            self.altseq,
            canonical,
        )


@dataclass(eq=True, frozen=True)
class _RnaSmallVariant(_RnaPosition, _SmallVariant):
    """Base class for RNA small variants."""

    @classmethod
    def from_rna(cls, rna: RnaPosition, refseq: str, altseq: str) -> List[_RnaSmallVariant]:
        """Convert an RNA position plus ref/alt nucleotides into an RNA variant.

        Args:
            rna (RnaPosition): RNA position
            refseq (str): Reference allele
            altseq (str): Alternate allele

        Raises:
            NotImplementedError: An unsupported combination of reference/alternate alleles was given
            ValueError: The given reference allele does not match the annotated reference allele

        Returns:
            List[_RnaSmallVariant]: One or more RNA variants
        """
        variant_list = []

        ref_annotated = rna.sequence()
        for ref, alt in product(expand_nt(refseq), expand_nt(altseq)):
            # Assert that the given ref matches the annotated one
            if ref != ref_annotated:
                raise ValueError(
                    f"Given ref allele '{ref}' does not match annotated ref allele '{ref_annotated}' for {rna}"
                )

            # Trim bases that are unchanged between the ref and alt alleles
            new_ref, new_alt, new_start, new_end = collapse_seq_change(ref, alt)
            start = rna.start + new_start
            end = rna.end - new_end

            # Determine the type of variant
            if is_deletion(new_ref, new_alt):
                variant = RnaDeletion.copy_from(
                    rna, start=start, end=end, refseq=new_ref, altseq=new_alt
                )
            elif is_delins(new_ref, new_alt):
                variant = RnaDelins.copy_from(
                    rna, start=start, end=end, refseq=new_ref, altseq=new_alt
                )
            elif is_duplication(new_ref, new_alt):
                variant = RnaDuplication.copy_from(
                    rna, start=start, end=end, refseq=new_ref, altseq=new_alt
                )
            elif is_insertion(new_ref, new_alt):
                variant = RnaInsertion.copy_from(
                    rna, start=start, end=end, refseq=new_ref, altseq=new_alt
                )
            elif is_substitution(new_ref, new_alt):
                variant = RnaSubstitution.copy_from(
                    rna, start=start, end=end, refseq=new_ref, altseq=new_alt
                )
            else:
                raise NotImplementedError(
                    f"Unrecognized sequence change '{new_ref}/{new_alt}' at {rna}"
                )

            variant_list.append(variant)

        return variant_list

    def to_cdna(self, canonical: bool = False) -> List[_CdnaSmallVariant]:
        """Map this variants to one more cDNA variants.

        Args:
            canonical (bool, optional): Only map to the canonical transcript. Defaults to False.

        Returns:
            List[_CdnaSmallVariant]: One or more cDNA variants.
        """
        return self._data._rna_to_cdna_variant(
            [self.transcript_id],
            self.start,
            self.start_offset,
            self.end,
            self.end_offset,
            [self.strand],
            self.refseq,
            self.altseq,
            canonical,
        )

    def to_dna(self, canonical: bool = False) -> List[_DnaSmallVariant]:
        """Map this variants to one more DNA variants.

        Args:
            canonical (bool, optional): Only map to the canonical transcript. Defaults to False.

        Returns:
            List[_DnaSmallVariant]: One or more DNA variants.
        """
        return self._data._rna_to_dna_variant(
            [self.transcript_id],
            self.start,
            self.start_offset,
            self.end,
            self.end_offset,
            [self.strand],
            self.refseq,
            self.altseq,
            canonical,
        )

    def to_exon(self, canonical: bool = False) -> List[_ExonSmallVariant]:
        """Map this variants to one more exon variants.

        Args:
            canonical (bool, optional): Only map to the canonical transcript. Defaults to False.

        Returns:
            List[_ExonSmallVariant]: One or more exon variants.
        """
        raise NotImplementedError()  # TODO

    def to_protein(self, canonical: bool = False) -> List[_ProteinSmallVariant]:
        """Map this variants to one more protein variants.

        Args:
            canonical (bool, optional): Only map to the canonical transcript. Defaults to False.

        Returns:
            List[_ProteinSmallVariant]: One or more protein variants.
        """
        return self._data._rna_to_protein_variant(
            [self.transcript_id],
            self.start,
            self.start_offset,
            self.end,
            self.end_offset,
            [self.strand],
            self.refseq,
            self.altseq,
            canonical,
        )

    def to_rna(self, canonical: bool = False) -> List[_RnaSmallVariant]:
        """Map this variants to one more RNA variants.

        Args:
            canonical (bool, optional): Only map to the canonical transcript. Defaults to False.

        Returns:
            List[_RnaSmallVariant]: One or more RNA variants.
        """
        return self._data._rna_to_rna_variant(
            [self.transcript_id],
            self.start,
            self.start_offset,
            self.end,
            self.end_offset,
            [self.strand],
            self.refseq,
            self.altseq,
            canonical,
        )


# -------------------------------------------------------------------------------------------------
# Deletion classes
# -------------------------------------------------------------------------------------------------
@dataclass(eq=True, frozen=True)
class _Deletion(_SmallVariant):
    """Base class for deletion variants."""

    @property
    def type(self) -> str:
        """Get the variant type.

        Returns:
            str: Type of variant
        """
        return DELETION


@dataclass(eq=True, frozen=True)
class CdnaDeletion(_Deletion, _CdnaSmallVariant):
    """Stores information on a cDNA deletion variant and maps to other position types."""

    def __str__(self) -> str:
        start = format_hgvs_position(self.start, self.start_offset)
        end = format_hgvs_position(self.end, self.end_offset)
        if start == end:
            return f"{self.transcript_id}:c.{start}del"
        else:
            return f"{self.transcript_id}:c.{start}_{end}del"


@dataclass(eq=True, frozen=True)
class DnaDeletion(_Deletion, _DnaSmallVariant):
    """Stores information on a DNA deletion variant and maps to other position types."""

    def __str__(self) -> str:
        start = format_hgvs_position(self.start, self.start_offset)
        end = format_hgvs_position(self.end, self.end_offset)
        if start == end:
            return f"{self.contig_id}:g.{start}del"
        else:
            return f"{self.contig_id}:g.{start}_{end}del"


@dataclass(eq=True, frozen=True)
class ProteinDeletion(_Deletion, _ProteinSmallVariant):
    """Stores information on a protein deletion variant and maps to other position types."""

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
class RnaDeletion(_Deletion, _RnaSmallVariant):
    """Stores information on an RNA deletion variant and maps to other position types."""

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
class _Delins(_SmallVariant):
    """Base class for delins variants."""

    @property
    def type(self) -> str:
        """Get the variant type.

        Returns:
            str: Type of variant
        """
        return DELINS


@dataclass(eq=True, frozen=True)
class CdnaDelins(_Delins, _CdnaSmallVariant):
    """Stores information on a cDNA delins variant and maps to other position types."""

    def __str__(self) -> str:
        start = format_hgvs_position(self.start, self.start_offset)
        end = format_hgvs_position(self.end, self.end_offset)
        if start == end:
            return f"{self.transcript_id}:c.{start}delins{self.altseq}"
        else:
            return f"{self.transcript_id}:c.{start}_{end}delins{self.altseq}"


@dataclass(eq=True, frozen=True)
class DnaDelins(_Delins, _DnaSmallVariant):
    """Stores information on a DNA delins variant and maps to other position types."""

    def __str__(self) -> str:
        start = format_hgvs_position(self.start, self.start_offset)
        end = format_hgvs_position(self.end, self.end_offset)
        if start == end:
            return f"{self.contig_id}:g.{start}delins{self.altseq}"
        else:
            return f"{self.contig_id}:g.{start}_{end}delins{self.altseq}"


@dataclass(eq=True, frozen=True)
class ProteinDelins(_Delins, _ProteinSmallVariant):
    """Stores information on a protein delins variant and maps to other position types."""

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
class RnaDelins(_Delins, _RnaSmallVariant):
    """Stores information on an RNA delins variant and maps to other position types."""

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
class _Duplication(_SmallVariant):
    """Base class for duplication variants."""

    @property
    def type(self) -> str:
        """Get the variant type.

        Returns:
            str: Type of variant
        """
        return DUPLICATION


@dataclass(eq=True, frozen=True)
class CdnaDuplication(_Duplication, _CdnaSmallVariant):
    """Stores information on a cDNA duplication variant and maps to other position types."""

    def __str__(self) -> str:
        start = format_hgvs_position(self.start, self.start_offset)
        end = format_hgvs_position(self.end, self.end_offset)
        if start == end:
            return f"{self.transcript_id}:c.{start}dup"
        else:
            return f"{self.transcript_id}:c.{start}_{end}dup"


@dataclass(eq=True, frozen=True)
class DnaDuplication(_Duplication, _DnaSmallVariant):
    """Stores information on a DNA duplication variant and maps to other position types."""

    def __str__(self) -> str:
        start = format_hgvs_position(self.start, self.start_offset)
        end = format_hgvs_position(self.end, self.end_offset)
        if start == end:
            return f"{self.contig_id}:g.{start}dup"
        else:
            return f"{self.contig_id}:g.{start}_{end}dup"


@dataclass(eq=True, frozen=True)
class ProteinDuplication(_Duplication, _ProteinSmallVariant):
    """Stores information on a protein duplication variant and maps to other position types."""

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
class RnaDuplication(_Duplication, _RnaSmallVariant):
    """Stores information on an RNA duplication variant and maps to other position types."""

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
class _Frameshift(_SmallVariant):
    """Base class for insertion variants."""

    @property
    def type(self) -> str:
        """Get the variant type.

        Returns:
            str: Type of variant
        """
        return FRAMESHIFT


@dataclass(eq=True, frozen=True)
class ProteinFrameshift(_Frameshift, _ProteinSmallVariant):
    """Stores information on a protein frameshift variant and maps to other position types."""

    def to_cdna(self, canonical: bool = False) -> List[_CdnaSmallVariant]:
        """Map this variants to one more cDNA variants.

        Args:
            canonical (bool, optional): Only map to the canonical transcript. Defaults to False.

        Returns:
            List[_CdnaSmallVariant]: One or more cDNA variants.
        """
        raise NotImplementedError()  # TODO

    def to_dna(self, canonical: bool = False) -> List[_DnaSmallVariant]:
        """Map this variants to one more DNA variants.

        Args:
            canonical (bool, optional): Only map to the canonical transcript. Defaults to False.

        Returns:
            List[_DnaSmallVariant]: One or more DNA variants.
        """
        raise NotImplementedError()  # TODO

    def to_exon(self, canonical: bool = False) -> List[_ExonSmallVariant]:
        """Map this variants to one more exon variants.

        Args:
            canonical (bool, optional): Only map to the canonical transcript. Defaults to False.

        Returns:
            List[_ExonSmallVariant]: One or more exon variants.
        """
        raise NotImplementedError()  # TODO

    def to_protein(self, canonical: bool = False) -> List[_ProteinSmallVariant]:
        """Map this variants to one more protein variants.

        Args:
            canonical (bool, optional): Only map to the canonical transcript. Defaults to False.

        Returns:
            List[_ProteinSmallVariant]: One or more protein variants.
        """
        raise NotImplementedError()  # TODO

    def to_rna(self, canonical: bool = False) -> List[_RnaSmallVariant]:
        """Map this variants to one more RNA variants.

        Args:
            canonical (bool, optional): Only map to the canonical transcript. Defaults to False.

        Returns:
            List[_RnaSmallVariant]: One or more RNA variants.
        """
        raise NotImplementedError()  # TODO


# -------------------------------------------------------------------------------------------------
# Insertion classes
# -------------------------------------------------------------------------------------------------
@dataclass(eq=True, frozen=True)
class _Insertion(_SmallVariant):
    """Base class for insertion variants."""

    @property
    def type(self) -> str:
        """Get the variant type.

        Returns:
            str: Type of variant
        """
        return INSERTION


@dataclass(eq=True, frozen=True)
class CdnaInsertion(_Insertion, _CdnaSmallVariant):
    """Stores information on a cDNA insertion variant and maps to other position types."""

    def __str__(self) -> str:
        start = format_hgvs_position(self.start, self.start_offset)
        end = format_hgvs_position(self.end, self.end_offset)
        altseq = self.altseq[1:-1]
        return f"{self.transcript_id}:c.{start}_{end}ins{altseq}"


@dataclass(eq=True, frozen=True)
class DnaInsertion(_Insertion, _DnaSmallVariant):
    """Stores information on a DNA insertion variant and maps to other position types."""

    def __str__(self) -> str:
        start = format_hgvs_position(self.start, self.start_offset)
        end = format_hgvs_position(self.end, self.end_offset)
        altseq = self.altseq[1:-1]
        return f"{self.contig_id}:g.{start}_{end}ins{altseq}"


@dataclass(eq=True, frozen=True)
class ProteinInsertion(_Insertion, _ProteinSmallVariant):
    """Stores information on a protein insertion variant and maps to other position types."""

    def __str__(self) -> str:
        start = format_hgvs_position(self.start, self.start_offset)
        end = format_hgvs_position(self.end, self.end_offset)
        start_seq = self.refseq[0]
        end_seq = self.refseq[-1]
        altseq = self.altseq[1:-1]
        return f"{self.protein_id}:p.{start_seq}{start}_{end_seq}{end}ins{altseq}"


@dataclass(eq=True, frozen=True)
class RnaInsertion(_Insertion, _RnaSmallVariant):
    """Stores information on an RNA insertion variant and maps to other position types."""

    def __str__(self) -> str:
        start = format_hgvs_position(self.start, self.start_offset)
        end = format_hgvs_position(self.end, self.end_offset)
        altseq = self.altseq[1:-1]
        return f"{self.transcript_id}:r.{start}_{end}ins{altseq}"


# -------------------------------------------------------------------------------------------------
# Substitution classes
# -------------------------------------------------------------------------------------------------
@dataclass(eq=True, frozen=True)
class _Substitution:
    """Base class for substitution variants."""

    refseq: str
    altseq: str

    @property
    def type(self) -> str:
        """Get the variant type.

        Returns:
            str: Type of variant
        """
        return SUBTITUTION


@dataclass(eq=True, frozen=True)
class CdnaSubstitution(_Substitution, _CdnaSmallVariant):
    """Stores information on a cDNA substitution variant and maps to other position types."""

    def __str__(self) -> str:
        start = format_hgvs_position(self.start, self.start_offset)
        return f"{self.transcript_id}:c.{start}{self.refseq}>{self.altseq}"


@dataclass(eq=True, frozen=True)
class DnaSubstitution(_Substitution, _DnaSmallVariant):
    """Stores information on a DNA substitution variant and maps to other position types."""

    def __str__(self) -> str:
        start = format_hgvs_position(self.start, self.start_offset)
        return f"{self.contig_id}:g.{start}{self.refseq}>{self.altseq}"


@dataclass(eq=True, frozen=True)
class ProteinSubstitution(_Substitution, _ProteinSmallVariant):
    """Stores information on a protein substitution variant and maps to other position types."""

    def __str__(self) -> str:
        start = format_hgvs_position(self.start, self.start_offset)
        return f"{self.protein_id}:p.{self.refseq}{start}{self.altseq}"


@dataclass(eq=True, frozen=True)
class RnaSubstitution(_Substitution, _RnaSmallVariant):
    """Stores information on an RNA substitution variant and maps to other position types."""

    def __str__(self) -> str:
        start = format_hgvs_position(self.start, self.start_offset)
        return f"{self.transcript_id}:r.{start}{self.refseq}>{self.altseq}"


# -------------------------------------------------------------------------------------------------
# Type hint variables
# -------------------------------------------------------------------------------------------------
MappablePositionOrSmallVariant = TypeVar(
    "MappablePositionOrSmallVariant",
    bound=Union[
        CdnaPosition,
        DnaPosition,
        ExonPosition,
        ProteinPosition,
        RnaPosition,
        _CdnaSmallVariant,
        _DnaSmallVariant,
        _ExonSmallVariant,
        _ProteinSmallVariant,
        _RnaSmallVariant,
    ],
)


# -------------------------------------------------------------------------------------------------
# Fusion classes
# -------------------------------------------------------------------------------------------------
@dataclass(eq=True, frozen=True)
class _Fusion(_Variant):
    """Base class for fusion variants."""

    breakpoint1: Union[CdnaPosition, DnaPosition, ExonPosition, ProteinPosition, RnaPosition]
    breakpoint2: Union[CdnaPosition, DnaPosition, ExonPosition, ProteinPosition, RnaPosition]

    @classmethod
    def copy_from(cls, fusion: _Fusion, **kwargs):
        """Initialize a new fusion object by copying values from another fusion object.

        Args:
            fusion (_Fusion): _Fusion object to copy attributes from
            **kwargs: Keyword arguments with keys that match the attribute names of this fusion
                class will override attributes from `fusion`

        Returns:
            a new fusion object of the same class as the class that calls this method
        """
        return cls(
            **{
                **{k: v for k, v in fusion.asdict().items() if k in [i.name for i in fields(cls)]},
                **kwargs,
            }
        )

    def __getitem__(self, item: Any) -> Any:
        return getattr(self, item)

    def __lt__(self, other: _Position) -> bool:
        return str(self) < str(other)

    def __str__(self) -> str:
        return f"{self.breakpoint1}::{self.breakpoint2}"

    @property
    def is_cdna(self) -> bool:
        """Check if this fusion is between cDNA.

        Returns:
            bool: True if the fusion is between cDNA else False
        """
        return False

    @property
    def is_dna(self) -> bool:
        """Check if this position is on DNA.

        Returns:
            bool: True if the position is on DNA else False
        """
        return False

    @property
    def is_exon(self) -> bool:
        """Check if this fusion is between exon.

        Returns:
            bool: True if the fusion is between exon else False
        """
        return False

    @property
    def is_protein(self) -> bool:
        """Check if this fusion is between protein.

        Returns:
            bool: True if the fusion is between protein else False
        """
        return False

    @property
    def is_rna(self) -> bool:
        """Check if this fusion is between RNA.

        Returns:
            bool: True if the fusion is between RNA else False
        """
        return False

    @property
    def on_negative_strand(self) -> bool:
        """Check if one or both breakpoints originate from the negative strand of the genome.

        Returns:
            bool: True if if one or both breakpoints are on the negative strand else False
        """
        return "-" in (self.breakpoint1.strand, self.breakpoint2.strand)

    @property
    def on_positive_strand(self) -> bool:
        """Check if one or both breakpoints originate from the positive strand of the genome.

        Returns:
            bool: True if if one or both breakpoints are on the positive strand else False
        """
        return "+" in (self.breakpoint1.strand, self.breakpoint2.strand)

    @property
    def type(self) -> str:
        """Get the variant type.

        Returns:
            str: Type of variant
        """
        return FUSION

    def asdict(self) -> Dict[str, Any]:
        """Convert this position's attributes to a dictionary.

        Returns:
            Dict[str, Any]: Dictionary of attribute names and corresponding values
        """
        return {f.name: self[f.name] for f in fields(self)}

    def sequence(self) -> str:
        """Return the reference sequence from this position's start to it's end, inclusive.

        Returns:
            str: The reference sequence
        """
        raise NotImplementedError()  # Defined by inheriting class

    def to_cdna(self, canonical: bool = False) -> List[CdnaFusion]:
        """Map this fusion to one more cDNA fusions.

        Args:
            canonical (bool, optional): Only map to the canonical transcript. Defaults to False.

        Returns:
            List[CdnaFusion]: One or more cDNA fusions.
        """
        result = []

        breakpoint1 = self.breakpoint1.to_cdna(canonical=canonical)
        breakpoint2 = self.breakpoint2.to_cdna(canonical=canonical)
        for b1, b2 in product(breakpoint1, breakpoint2):
            result.append(CdnaFusion(b1, b2))

        return result

    def to_dna(self, canonical: bool = False) -> List[DnaFusion]:
        """Map this fusion to one more DNA fusions.

        Args:
            canonical (bool, optional): Only map to the canonical transcript. Defaults to False.

        Returns:
            List[DnaFusion]: One or more DNA fusions.
        """
        result = []

        breakpoint1 = self.breakpoint1.to_dna(canonical=canonical)
        breakpoint2 = self.breakpoint2.to_dna(canonical=canonical)
        for b1, b2 in product(breakpoint1, breakpoint2):
            result.append(DnaFusion(b1, b2))

        return result

    def to_exon(self, canonical: bool = False) -> List[ExonFusion]:
        """Map this fusion to one more exon fusions.

        Args:
            canonical (bool, optional): Only map to the canonical transcript. Defaults to False.

        Returns:
            List[ExonFusion]: One or more exon fusions.
        """
        result = []

        breakpoint1 = self.breakpoint1.to_exon(canonical=canonical)
        breakpoint2 = self.breakpoint2.to_exon(canonical=canonical)
        for b1, b2 in product(breakpoint1, breakpoint2):
            result.append(ExonFusion(b1, b2))

        return result

    def to_protein(self, canonical: bool = False) -> List[ProteinFusion]:
        """Map this fusion to one more protein fusions.

        Args:
            canonical (bool, optional): Only map to the canonical transcript. Defaults to False.

        Returns:
            List[ProteinFusion]: One or more protein fusions.
        """
        result = []

        breakpoint1 = self.breakpoint1.to_protein(canonical=canonical)
        breakpoint2 = self.breakpoint2.to_protein(canonical=canonical)
        for b1, b2 in product(breakpoint1, breakpoint2):
            result.append(ProteinFusion(b1, b2))

        return result

    def to_rna(self, canonical: bool = False) -> List[RnaFusion]:
        """Map this fusion to one more RNA fusions.

        Args:
            canonical (bool, optional): Only map to the canonical transcript. Defaults to False.

        Returns:
            List[RnaFusion]: One or more RNA fusions.
        """
        result = []

        breakpoint1 = self.breakpoint1.to_rna(canonical=canonical)
        breakpoint2 = self.breakpoint2.to_rna(canonical=canonical)
        for b1, b2 in product(breakpoint1, breakpoint2):
            result.append(RnaFusion(b1, b2))

        return result


@dataclass(eq=True, frozen=True)
class CdnaFusion(_Fusion):
    """Stores information on a cDNA fusion variant and maps to other position types."""

    @property
    def is_cdna(self) -> bool:
        """Check if this fusion is between cDNA.

        Returns:
            bool: True if the fusion is between cDNA else False
        """
        return True


@dataclass(eq=True, frozen=True)
class DnaFusion(_Fusion):
    """Stores information on a DNA fusion variant and maps to other position types."""

    @property
    def is_dna(self) -> bool:
        """Check if this position is on DNA.

        Returns:
            bool: True if the position is on DNA else False
        """
        return True


@dataclass(eq=True, frozen=True)
class ExonFusion(_Fusion):
    """Stores information on an exon fusion variant and maps to other position types."""

    @property
    def is_exon(self) -> bool:
        """Check if this fusion is between exon.

        Returns:
            bool: True if the fusion is between exon else False
        """
        return True


@dataclass(eq=True, frozen=True)
class ProteinFusion(_Fusion):
    """Stores information on a protein fusion variant and maps to other position types."""

    @property
    def is_protein(self) -> bool:
        """Check if this fusion is between protein.

        Returns:
            bool: True if the fusion is between protein else False
        """
        return True


@dataclass(eq=True, frozen=True)
class RnaFusion(_Fusion):
    """Stores information on a RNA fusion variant and maps to other position types."""

    @property
    def is_rna(self) -> bool:
        """Check if this fusion is between RNA.

        Returns:
            bool: True if the fusion is between RNA else False
        """
        return True
