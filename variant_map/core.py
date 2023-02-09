from __future__ import annotations

from dataclasses import dataclass, fields
from itertools import product
from math import floor
from typing import Any, Callable, Dict, List, Optional, Tuple, TypeVar, Union

import pandas as pd
from gtfparse import read_gtf
from pyfaidx import Fasta

from .constants import (
    CONTIG_ID,
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
    calc_cdna_to_protein,
    collapse_seq_change,
    expand_nt,
    format_hgvs_position,
    is_deletion,
    is_delins,
    is_duplication,
    is_frameshift,
    is_insertion,
    is_substitution,
    reverse_complement,
    reverse_translate,
    split_by_codon,
    strip_version,
)


# -------------------------------------------------------------------------------------------------
# Position classes
# -------------------------------------------------------------------------------------------------
@dataclass(eq=True, frozen=True)
class _Position:
    """Base class for position objects."""

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
        raise NotImplementedError()  # Defined be inheriting classes

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
        raise NotImplementedError()  # Defined be inheriting classes


@dataclass(eq=True, frozen=True)
class CdnaPosition(_Position):
    """Stores information on cDNA position objects."""

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
class DnaPosition(_Position):
    """Stores information on DNA position objects."""

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
class ExonPosition(_Position):
    """Stores information on exon position objects."""

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
class ProteinPosition(_Position):
    """Stores information on protein position objects."""

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
class RnaPosition(_Position):
    """Stores information on RNA position objects."""

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
# _MappablePosition classes
# -------------------------------------------------------------------------------------------------
@dataclass(eq=True, frozen=True)
class _MappablePosition:
    """Base class for mappable position objects."""

    def to_cdna(self, canonical: bool = False) -> List[CdnaMappablePosition]:
        """Map this position to one more cDNA positions.

        Args:
            canonical (bool, optional): Only map to the canonical transcript. Defaults to False.

        Returns:
            List[CdnaMappablePosition]: One or more cDNA positions.
        """
        raise NotImplementedError()  # Defined be inheriting classes

    def to_dna(self, canonical: bool = False) -> List[DnaMappablePosition]:
        """Map this position to one more DNA positions.

        Args:
            canonical (bool, optional): Only map to the canonical transcript. Defaults to False.

        Returns:
            List[DnaMappablePosition]: One or more DNA positions.
        """
        raise NotImplementedError()  # Defined be inheriting classes

    def to_exon(self, canonical: bool = False) -> List[ExonMappablePosition]:
        """Map this position to one more exon positions.

        Args:
            canonical (bool, optional): Only map to the canonical transcript. Defaults to False.

        Returns:
            List[ExonMappablePosition]: One or more exon positions.
        """
        raise NotImplementedError()  # Defined be inheriting classes

    def to_protein(self, canonical: bool = False) -> List[ProteinMappablePosition]:
        """Map this position to one more protein positions.

        Args:
            canonical (bool, optional): Only map to the canonical transcript. Defaults to False.

        Returns:
            List[ProteinMappablePosition]: One or more protein positions.
        """
        raise NotImplementedError()  # Defined be inheriting classes

    def to_rna(self, canonical: bool = False) -> List[RnaMappablePosition]:
        """Map this position to one more RNA positions.

        Args:
            canonical (bool, optional): Only map to the canonical transcript. Defaults to False.

        Returns:
            List[RnaMappablePosition]: One or more RNA positions.
        """
        raise NotImplementedError()  # Defined be inheriting classes


@dataclass(eq=True, frozen=True)
class CdnaMappablePosition(CdnaPosition, _MappablePosition):
    """Stores information on a cDNA position and converts to other position types."""

    def to_cdna(self, canonical: bool = False) -> List[CdnaMappablePosition]:
        """Map this position to one more cDNA positions.

        Args:
            canonical (bool, optional): Only map to the canonical transcript. Defaults to False.

        Returns:
            List[CdnaMappablePosition]: One or more cDNA positions.
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

    def to_dna(self, canonical: bool = False) -> List[DnaMappablePosition]:
        """Map this position to one more DNA positions.

        Args:
            canonical (bool, optional): Only map to the canonical transcript. Defaults to False.

        Returns:
            List[DnaMappablePosition]: One or more DNA positions.
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

    def to_exon(self, canonical: bool = False) -> List[ExonMappablePosition]:
        """Map this position to one more exon positions.

        Args:
            canonical (bool, optional): Only map to the canonical transcript. Defaults to False.

        Returns:
            List[ExonMappablePosition]: One or more exon positions.
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

    def to_protein(self, canonical: bool = False) -> List[ProteinMappablePosition]:
        """Map this position to one more protein positions.

        Args:
            canonical (bool, optional): Only map to the canonical transcript. Defaults to False.

        Returns:
            List[ProteinMappablePosition]: One or more protein positions.
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

    def to_rna(self, canonical: bool = False) -> List[RnaMappablePosition]:
        """Map this position to one more RNA positions.

        Args:
            canonical (bool, optional): Only map to the canonical transcript. Defaults to False.

        Returns:
            List[RnaMappablePosition]: One or more RNA positions.
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
class DnaMappablePosition(DnaPosition, _MappablePosition):
    """Stores information on a DNA position and converts to other position types."""

    def to_cdna(self, canonical: bool = False) -> List[CdnaMappablePosition]:
        """Map this position to one more cDNA positions.

        Args:
            canonical (bool, optional): Only map to the canonical transcript. Defaults to False.

        Returns:
            List[CdnaMappablePosition]: One or more cDNA positions.
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

    def to_dna(self, canonical: bool = False) -> List[DnaMappablePosition]:
        """Map this position to one more DNA positions.

        Args:
            canonical (bool, optional): Only map to the canonical transcript. Defaults to False.

        Returns:
            List[DnaMappablePosition]: One or more DNA positions.
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

    def to_exon(self, canonical: bool = False) -> List[ExonMappablePosition]:
        """Map this position to one more exon positions.

        Args:
            canonical (bool, optional): Only map to the canonical transcript. Defaults to False.

        Returns:
            List[ExonMappablePosition]: One or more exon positions.
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

    def to_protein(self, canonical: bool = False) -> List[ProteinMappablePosition]:
        """Map this position to one more protein positions.

        Args:
            canonical (bool, optional): Only map to the canonical transcript. Defaults to False.

        Returns:
            List[ProteinMappablePosition]: One or more protein positions.
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

    def to_rna(self, canonical: bool = False) -> List[RnaMappablePosition]:
        """Map this position to one more RNA positions.

        Args:
            canonical (bool, optional): Only map to the canonical transcript. Defaults to False.

        Returns:
            List[RnaMappablePosition]: One or more RNA positions.
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
class ExonMappablePosition(ExonPosition, _MappablePosition):
    """Stores information on an exon position and converts to other position types."""

    def to_cdna(self, canonical: bool = False) -> List[CdnaMappablePosition]:
        """Map this position to one more cDNA positions.

        Args:
            canonical (bool, optional): Only map to the canonical transcript. Defaults to False.

        Returns:
            List[CdnaMappablePosition]: One or more cDNA positions.
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

    def to_dna(self, canonical: bool = False) -> List[DnaMappablePosition]:
        """Map this position to one more DNA positions.

        Args:
            canonical (bool, optional): Only map to the canonical transcript. Defaults to False.

        Returns:
            List[DnaMappablePosition]: One or more DNA positions.
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

    def to_exon(self, canonical: bool = False) -> List[ExonMappablePosition]:
        """Map this position to one more exon positions.

        Args:
            canonical (bool, optional): Only map to the canonical transcript. Defaults to False.

        Returns:
            List[ExonMappablePosition]: One or more exon positions.
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

    def to_protein(self, canonical: bool = False) -> List[ProteinMappablePosition]:
        """Map this position to one more protein positions.

        Args:
            canonical (bool, optional): Only map to the canonical transcript. Defaults to False.

        Returns:
            List[ProteinMappablePosition]: One or more protein positions.
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

    def to_rna(self, canonical: bool = False) -> List[RnaMappablePosition]:
        """Map this position to one more RNA positions.

        Args:
            canonical (bool, optional): Only map to the canonical transcript. Defaults to False.

        Returns:
            List[RnaMappablePosition]: One or more RNA positions.
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
class ProteinMappablePosition(ProteinPosition, _MappablePosition):
    """Stores information on a protein position and converts to other position types."""

    def to_cdna(self, canonical: bool = False) -> List[CdnaMappablePosition]:
        """Map this position to one more cDNA positions.

        Args:
            canonical (bool, optional): Only map to the canonical transcript. Defaults to False.

        Returns:
            List[CdnaMappablePosition]: One or more cDNA positions.
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

    def to_dna(self, canonical: bool = False) -> List[DnaMappablePosition]:
        """Map this position to one more DNA positions.

        Args:
            canonical (bool, optional): Only map to the canonical transcript. Defaults to False.

        Returns:
            List[DnaMappablePosition]: One or more DNA positions.
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

    def to_exon(self, canonical: bool = False) -> List[ExonMappablePosition]:
        """Map this position to one more exon positions.

        Args:
            canonical (bool, optional): Only map to the canonical transcript. Defaults to False.

        Returns:
            List[ExonMappablePosition]: One or more exon positions.
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

    def to_protein(self, canonical: bool = False) -> List[ProteinMappablePosition]:
        """Map this position to one more protein positions.

        Args:
            canonical (bool, optional): Only map to the canonical transcript. Defaults to False.

        Returns:
            List[ProteinMappablePosition]: One or more protein positions.
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

    def to_rna(self, canonical: bool = False) -> List[RnaMappablePosition]:
        """Map this position to one more RNA positions.

        Args:
            canonical (bool, optional): Only map to the canonical transcript. Defaults to False.

        Returns:
            List[RnaMappablePosition]: One or more RNA positions.
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
class RnaMappablePosition(RnaPosition, _MappablePosition):
    """Stores information on an RNA position and converts to other position types."""

    def to_cdna(self, canonical: bool = False) -> List[CdnaMappablePosition]:
        """Map this position to one more cDNA positions.

        Args:
            canonical (bool, optional): Only map to the canonical transcript. Defaults to False.

        Returns:
            List[CdnaMappablePosition]: One or more cDNA positions.
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

    def to_dna(self, canonical: bool = False) -> List[DnaMappablePosition]:
        """Map this position to one more DNA positions.

        Args:
            canonical (bool, optional): Only map to the canonical transcript. Defaults to False.

        Returns:
            List[DnaMappablePosition]: One or more DNA positions.
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

    def to_exon(self, canonical: bool = False) -> List[ExonMappablePosition]:
        """Map this position to one more exon positions.

        Args:
            canonical (bool, optional): Only map to the canonical transcript. Defaults to False.

        Returns:
            List[ExonMappablePosition]: One or more exon positions.
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

    def to_protein(self, canonical: bool = False) -> List[ProteinMappablePosition]:
        """Map this position to one more protein positions.

        Args:
            canonical (bool, optional): Only map to the canonical transcript. Defaults to False.

        Returns:
            List[ProteinMappablePosition]: One or more protein positions.
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

    def to_rna(self, canonical: bool = False) -> List[RnaMappablePosition]:
        """Map this position to one more RNA positions.

        Args:
            canonical (bool, optional): Only map to the canonical transcript. Defaults to False.

        Returns:
            List[RnaMappablePosition]: One or more RNA positions.
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
class Variant:
    """Base class for variant objects."""

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
        raise NotImplementedError()  # Defined be inheriting classes

    def to_cdna(self) -> List:
        """Map this variant to one more cDNA variants.

        Args:
            canonical (bool, optional): Only map to the canonical transcript. Defaults to False.

        Returns:
            List: One or more cDNA variants.
        """
        raise NotImplementedError()  # Defined be inheriting classes

    def to_dna(self) -> List:
        """Map this variant to one more DNA variants.

        Args:
            canonical (bool, optional): Only map to the canonical transcript. Defaults to False.

        Returns:
            List: One or more DNA variants.
        """
        raise NotImplementedError()  # Defined be inheriting classes

    def to_exon(self) -> List:
        """Map this variant to one more exon variants.

        Args:
            canonical (bool, optional): Only map to the canonical transcript. Defaults to False.

        Returns:
            List: One or more exon variants.
        """
        raise NotImplementedError()  # Defined be inheriting classes

    def to_protein(self) -> List:
        """Map this variant to one more protein variants.

        Args:
            canonical (bool, optional): Only map to the canonical transcript. Defaults to False.

        Returns:
            List: One or more protein variants.
        """
        raise NotImplementedError()  # Defined be inheriting classes

    def to_rna(self) -> List:
        """Map this variant to one more RNA variants.

        Args:
            canonical (bool, optional): Only map to the canonical transcript. Defaults to False.

        Returns:
            List: One or more RNA variants.
        """
        raise NotImplementedError()  # Defined be inheriting classes


# -------------------------------------------------------------------------------------------------
# SmallVariant classes
# -------------------------------------------------------------------------------------------------
@dataclass(eq=True, frozen=True)
class SmallVariant(Variant):
    """Base class for small variant objects."""

    refseq: str
    altseq: str


@dataclass(eq=True, frozen=True)
class CdnaSmallVariant(CdnaPosition, SmallVariant):
    """Base class for cDNA variant objects."""

    @classmethod
    def from_cdna(
        cls, cdna: CdnaMappablePosition, refseq: str, altseq: str
    ) -> List[CdnaSmallVariant]:
        """Convert a cDNA position plus ref/alt nucleotides into a cDNA variant object.

        Args:
            cdna (CdnaMappablePosition): Variant position
            refseq (str): reference allele
            altseq (str): alternate allele

        Raises:
            NotImplementedError: An unsupported combination of reference/alternate alleles was given
            ValueError: The given reference allele does not match the annotated reference allele

        Returns:
            List[CdnaSmallVariant]: One or more cDNA variants
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
    def from_protein(
        cls, cdna: CdnaMappablePosition, refseq: str, altseq: str
    ) -> List[CdnaSmallVariant]:
        """Convert a cDNA position plus ref/alt amino acids into a cDNA variant object.

        Args:
            cdna (CdnaMappablePosition): cDNA position
            refseq (str): Reference allele
            altseq (str): Alternate allele

        Raises:
            NotImplementedError: An unsupported combination of reference/alternate alleles was given

        Returns:
            List[CdnaSmallVariant]: One or more cDNA variants
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

    def to_cdna(self, canonical: bool = False) -> List[CdnaSmallVariant]:
        """Map this variants to one more cDNA variants.

        Args:
            canonical (bool, optional): Only map to the canonical transcript. Defaults to False.

        Returns:
            List[CdnaSmallVariant]: One or more cDNA variants.
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

    def to_dna(self, canonical: bool = False) -> List[DnaSmallVariant]:
        """Map this variants to one more DNA variants.

        Args:
            canonical (bool, optional): Only map to the canonical transcript. Defaults to False.

        Returns:
            List[DnaSmallVariant]: One or more DNA variants.
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

    def to_exon(self, canonical: bool = False) -> List[ExonSmallVariant]:
        """Map this variants to one more exon variants.

        Args:
            canonical (bool, optional): Only map to the canonical transcript. Defaults to False.

        Returns:
            List[ExonSmallVariant]: One or more exon variants.
        """
        raise NotImplementedError()  # TODO

    def to_protein(self, canonical: bool = False) -> List[ProteinSmallVariant]:
        """Map this variants to one more protein variants.

        Args:
            canonical (bool, optional): Only map to the canonical transcript. Defaults to False.

        Returns:
            List[ProteinSmallVariant]: One or more protein variants.
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

    def to_rna(self, canonical: bool = False) -> List[RnaSmallVariant]:
        """Map this variants to one more RNA variants.

        Args:
            canonical (bool, optional): Only map to the canonical transcript. Defaults to False.

        Returns:
            List[RnaSmallVariant]: One or more RNA variants.
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
class DnaSmallVariant(DnaPosition, SmallVariant):
    """Base class for DNA variant objects."""

    @classmethod
    def from_dna(cls, dna: DnaMappablePosition, refseq: str, altseq: str) -> List[DnaSmallVariant]:
        """Convert a DNA position plus ref/alt nucleotides into a DNA variant object.

        Args:
            dna (DnaMappablePosition): DNA position
            refseq (str): Reference allele
            altseq (str): Alternate allele

        Raises:
            NotImplementedError: An unsupported combination of reference/alternate alleles was given

        Returns:
            List[DnaSmallVariant]: One or more DNA variants
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

    def to_cdna(self, canonical: bool = False) -> List[CdnaSmallVariant]:
        """Map this variants to one more cDNA variants.

        Args:
            canonical (bool, optional): Only map to the canonical transcript. Defaults to False.

        Returns:
            List[CdnaSmallVariant]: One or more cDNA variants.
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

    def to_dna(self, canonical: bool = False) -> List[DnaSmallVariant]:
        """Map this variants to one more DNA variants.

        Args:
            canonical (bool, optional): Only map to the canonical transcript. Defaults to False.

        Returns:
            List[DnaSmallVariant]: One or more DNA variants.
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

    def to_exon(self, canonical: bool = False) -> List[ExonSmallVariant]:
        """Map this variants to one more exon variants.

        Args:
            canonical (bool, optional): Only map to the canonical transcript. Defaults to False.

        Returns:
            List[ExonSmallVariant]: One or more exon variants.
        """
        raise NotImplementedError()  # TODO

    def to_protein(self, canonical: bool = False) -> List[ProteinSmallVariant]:
        """Map this variants to one more protein variants.

        Args:
            canonical (bool, optional): Only map to the canonical transcript. Defaults to False.

        Returns:
            List[ProteinSmallVariant]: One or more protein variants.
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

    def to_rna(self, canonical: bool = False) -> List[RnaSmallVariant]:
        """Map this variants to one more RNA variants.

        Args:
            canonical (bool, optional): Only map to the canonical transcript. Defaults to False.

        Returns:
            List[RnaSmallVariant]: One or more RNA variants.
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
class ExonSmallVariant(ExonPosition, SmallVariant):
    """Base class for exon variant objects."""

    @classmethod
    def from_exon(
        cls, exon: ExonMappablePosition, refseq: str, altseq: str
    ) -> List[ExonSmallVariant]:
        """Convert an exon position plus ref/alt nucleotides into an exon variant object.

        Args:
            exon (ExonMappablePosition): exon position
            refseq (str): Reference allele
            altseq (str): Alternate allele

        Raises:
            NotImplementedError: An unsupported combination of reference/alternate alleles was given

        Returns:
            List[ExonSmallVariant]: One or more exon variants
        """
        raise NotImplementedError()  # TODO


@dataclass(eq=True, frozen=True)
class ProteinSmallVariant(ProteinPosition, SmallVariant):
    """Base class for protein variant objects."""

    @classmethod
    def from_cdna(
        cls, cdna: CdnaMappablePosition, protein: ProteinMappablePosition, refseq: str, altseq: str
    ) -> List[ProteinSmallVariant]:
        """Convert a cDNA position plus ref/alt amino acids into a protein variant object.

        Args:
            cdna (CdnaMappablePosition): cDNA position
            refseq (str): Reference allele
            altseq (str): Alternate allele

        Raises:
            NotImplementedError: An unsupported combination of reference/alternate alleles was given

        Returns:
            List[ProteinSmallVariant]: One or more protein variants
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

    def to_cdna(self, canonical: bool = False) -> List[CdnaSmallVariant]:
        """Map this variants to one more cDNA variants.

        Args:
            canonical (bool, optional): Only map to the canonical transcript. Defaults to False.

        Returns:
            List[CdnaSmallVariant]: One or more cDNA variants.
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

    def to_dna(self, canonical: bool = False) -> List[DnaSmallVariant]:
        """Map this variants to one more DNA variants.

        Args:
            canonical (bool, optional): Only map to the canonical transcript. Defaults to False.

        Returns:
            List[DnaSmallVariant]: One or more DNA variants.
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

    def to_exon(self, canonical: bool = False) -> List[ExonSmallVariant]:
        """Map this variants to one more exon variants.

        Args:
            canonical (bool, optional): Only map to the canonical transcript. Defaults to False.

        Returns:
            List[ExonSmallVariant]: One or more exon variants.
        """
        raise NotImplementedError()  # TODO

    def to_protein(self, canonical: bool = False) -> List[ProteinSmallVariant]:
        """Map this variants to one more protein variants.

        Args:
            canonical (bool, optional): Only map to the canonical transcript. Defaults to False.

        Returns:
            List[ProteinSmallVariant]: One or more protein variants.
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

    def to_rna(self, canonical: bool = False) -> List[RnaSmallVariant]:
        """Map this variants to one more RNA variants.

        Args:
            canonical (bool, optional): Only map to the canonical transcript. Defaults to False.

        Returns:
            List[RnaSmallVariant]: One or more RNA variants.
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
class RnaSmallVariant(RnaPosition, SmallVariant):
    """Base class for RNA variant objects."""

    @classmethod
    def from_rna(cls, rna: RnaMappablePosition, refseq: str, altseq: str) -> List[RnaSmallVariant]:
        """Convert an RNA position plus ref/alt nucleotides into an RNA variant object.

        Args:
            rna (RnaMappablePosition): RNA position
            refseq (str): Reference allele
            altseq (str): Alternate allele

        Raises:
            NotImplementedError: An unsupported combination of reference/alternate alleles was given
            ValueError: The given reference allele does not match the annotated reference allele

        Returns:
            List[RnaSmallVariant]: One or more RNA variants
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

    def to_cdna(self, canonical: bool = False) -> List[CdnaSmallVariant]:
        """Map this variants to one more cDNA variants.

        Args:
            canonical (bool, optional): Only map to the canonical transcript. Defaults to False.

        Returns:
            List[CdnaSmallVariant]: One or more cDNA variants.
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

    def to_dna(self, canonical: bool = False) -> List[DnaSmallVariant]:
        """Map this variants to one more DNA variants.

        Args:
            canonical (bool, optional): Only map to the canonical transcript. Defaults to False.

        Returns:
            List[DnaSmallVariant]: One or more DNA variants.
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

    def to_exon(self, canonical: bool = False) -> List[ExonSmallVariant]:
        """Map this variants to one more exon variants.

        Args:
            canonical (bool, optional): Only map to the canonical transcript. Defaults to False.

        Returns:
            List[ExonSmallVariant]: One or more exon variants.
        """
        raise NotImplementedError()  # TODO

    def to_protein(self, canonical: bool = False) -> List[ProteinSmallVariant]:
        """Map this variants to one more protein variants.

        Args:
            canonical (bool, optional): Only map to the canonical transcript. Defaults to False.

        Returns:
            List[ProteinSmallVariant]: One or more protein variants.
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

    def to_rna(self, canonical: bool = False) -> List[RnaSmallVariant]:
        """Map this variants to one more RNA variants.

        Args:
            canonical (bool, optional): Only map to the canonical transcript. Defaults to False.

        Returns:
            List[RnaSmallVariant]: One or more RNA variants.
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
class Deletion(SmallVariant):
    """Base class for deletion variant objects."""

    @property
    def type(self) -> str:
        """Get the variant type.

        Returns:
            str: Type of variant
        """
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
    """Base class for delins variant objects."""

    @property
    def type(self) -> str:
        """Get the variant type.

        Returns:
            str: Type of variant
        """
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
    """Base class for duplication variant objects."""

    @property
    def type(self) -> str:
        """Get the variant type.

        Returns:
            str: Type of variant
        """
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
    """Base class for insertion variant objects."""

    @property
    def type(self) -> str:
        """Get the variant type.

        Returns:
            str: Type of variant
        """
        return FRAMESHIFT


@dataclass(eq=True, frozen=True)
class ProteinFrameshift(Frameshift, ProteinSmallVariant):
    """Stores information on a protein frameshift variant and converts to other position types."""

    def to_cdna(self, canonical: bool = False) -> List[CdnaSmallVariant]:
        """Map this variants to one more cDNA variants.

        Args:
            canonical (bool, optional): Only map to the canonical transcript. Defaults to False.

        Returns:
            List[CdnaSmallVariant]: One or more cDNA variants.
        """
        raise NotImplementedError()  # TODO

    def to_dna(self, canonical: bool = False) -> List[DnaSmallVariant]:
        """Map this variants to one more DNA variants.

        Args:
            canonical (bool, optional): Only map to the canonical transcript. Defaults to False.

        Returns:
            List[DnaSmallVariant]: One or more DNA variants.
        """
        raise NotImplementedError()  # TODO

    def to_exon(self, canonical: bool = False) -> List[ExonSmallVariant]:
        """Map this variants to one more exon variants.

        Args:
            canonical (bool, optional): Only map to the canonical transcript. Defaults to False.

        Returns:
            List[ExonSmallVariant]: One or more exon variants.
        """
        raise NotImplementedError()  # TODO

    def to_protein(self, canonical: bool = False) -> List[ProteinSmallVariant]:
        """Map this variants to one more protein variants.

        Args:
            canonical (bool, optional): Only map to the canonical transcript. Defaults to False.

        Returns:
            List[ProteinSmallVariant]: One or more protein variants.
        """
        raise NotImplementedError()  # TODO

    def to_rna(self, canonical: bool = False) -> List[RnaSmallVariant]:
        """Map this variants to one more RNA variants.

        Args:
            canonical (bool, optional): Only map to the canonical transcript. Defaults to False.

        Returns:
            List[RnaSmallVariant]: One or more RNA variants.
        """
        raise NotImplementedError()  # TODO


# -------------------------------------------------------------------------------------------------
# Insertion classes
# -------------------------------------------------------------------------------------------------
@dataclass(eq=True, frozen=True)
class Insertion(SmallVariant):
    """Base class for insertion variant objects."""

    @property
    def type(self) -> str:
        """Get the variant type.

        Returns:
            str: Type of variant
        """
        return INSERTION


@dataclass(eq=True, frozen=True)
class CdnaInsertion(Insertion, CdnaSmallVariant):
    """Stores information on a cDNA insertion variant and converts to other position types."""

    def __str__(self) -> str:
        start = format_hgvs_position(self.start, self.start_offset)
        end = format_hgvs_position(self.end, self.end_offset)
        altseq = self.altseq[1:-1]
        return f"{self.transcript_id}:c.{start}_{end}ins{altseq}"


@dataclass(eq=True, frozen=True)
class DnaInsertion(Insertion, DnaSmallVariant):
    """Stores information on a DNA insertion variant and converts to other position types."""

    def __str__(self) -> str:
        start = format_hgvs_position(self.start, self.start_offset)
        end = format_hgvs_position(self.end, self.end_offset)
        altseq = self.altseq[1:-1]
        return f"{self.contig_id}:g.{start}_{end}ins{altseq}"


@dataclass(eq=True, frozen=True)
class ProteinInsertion(Insertion, ProteinSmallVariant):
    """Stores information on a protein insertion variant and converts to other position types."""

    def __str__(self) -> str:
        start = format_hgvs_position(self.start, self.start_offset)
        end = format_hgvs_position(self.end, self.end_offset)
        start_seq = self.refseq[0]
        end_seq = self.refseq[-1]
        altseq = self.altseq[1:-1]
        return f"{self.protein_id}:p.{start_seq}{start}_{end_seq}{end}ins{altseq}"


@dataclass(eq=True, frozen=True)
class RnaInsertion(Insertion, RnaSmallVariant):
    """Stores information on an RNA insertion variant and converts to other position types."""

    def __str__(self) -> str:
        start = format_hgvs_position(self.start, self.start_offset)
        end = format_hgvs_position(self.end, self.end_offset)
        altseq = self.altseq[1:-1]
        return f"{self.transcript_id}:r.{start}_{end}ins{altseq}"


# -------------------------------------------------------------------------------------------------
# Substitution classes
# -------------------------------------------------------------------------------------------------
@dataclass(eq=True, frozen=True)
class Substitution:
    """Base class for substitution variant objects."""

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
# Type hint variables
# -------------------------------------------------------------------------------------------------
MappablePositionOrSmallVariant = TypeVar(
    "MappablePositionOrSmallVariant",
    bound=Union[
        CdnaMappablePosition,
        DnaMappablePosition,
        ExonMappablePosition,
        ProteinMappablePosition,
        RnaMappablePosition,
        CdnaSmallVariant,
        DnaSmallVariant,
        ExonSmallVariant,
        ProteinSmallVariant,
        RnaSmallVariant,
    ],
)


# -------------------------------------------------------------------------------------------------
# Fusion classes
# -------------------------------------------------------------------------------------------------
@dataclass(eq=True, frozen=True)
class Fusion(Variant):
    """Base class for fusion variant objects."""

    breakpoint1: Union[
        CdnaMappablePosition,
        DnaMappablePosition,
        ExonMappablePosition,
        ProteinMappablePosition,
        RnaMappablePosition,
    ]
    breakpoint2: Union[
        CdnaMappablePosition,
        DnaMappablePosition,
        ExonMappablePosition,
        ProteinMappablePosition,
        RnaMappablePosition,
    ]

    @classmethod
    def copy_from(cls, fusion: Fusion, **kwargs):
        """Initialize a new fusion object by copying values from another fusion object.

        Args:
            fusion (Fusion): Fusion object to copy attributes from
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
        raise NotImplementedError()  # Defined be inheriting classes

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
class CdnaFusion(Fusion):
    """Stores information on a cDNA fusion variant and converts to other position types."""

    @property
    def is_cdna(self) -> bool:
        """Check if this fusion is between cDNA.

        Returns:
            bool: True if the fusion is between cDNA else False
        """
        return True


@dataclass(eq=True, frozen=True)
class DnaFusion(Fusion):
    """Stores information on a DNA fusion variant and converts to other position types."""

    @property
    def is_dna(self) -> bool:
        """Check if this position is on DNA.

        Returns:
            bool: True if the position is on DNA else False
        """
        return True


@dataclass(eq=True, frozen=True)
class ExonFusion(Fusion):
    """Stores information on an exon fusion variant and converts to other position types."""

    @property
    def is_exon(self) -> bool:
        """Check if this fusion is between exon.

        Returns:
            bool: True if the fusion is between exon else False
        """
        return True


@dataclass(eq=True, frozen=True)
class ProteinFusion(Fusion):
    """Stores information on a protein fusion variant and converts to other position types."""

    @property
    def is_protein(self) -> bool:
        """Check if this fusion is between protein.

        Returns:
            bool: True if the fusion is between protein else False
        """
        return True


@dataclass(eq=True, frozen=True)
class RnaFusion(Fusion):
    """Stores information on a RNA fusion variant and converts to other position types."""

    @property
    def is_rna(self) -> bool:
        """Check if this fusion is between RNA.

        Returns:
            bool: True if the fusion is between RNA else False
        """
        return True


# -------------------------------------------------------------------------------------------------
# Core classes and methods
# -------------------------------------------------------------------------------------------------
class Core:
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
            gtf (str): Path to a GTF files with feature annotations
            cds (List[str]): List of paths to a FASTA files of CDS sequences
            dna (List[str]): List of paths to a FASTA files of DNA sequences
            peptide (List[str]): List of paths to a FASTA files of peptide sequences
            rna (List[str]): List of paths to a FASTA files of RNA sequences
            canonical_transcript (str, optional): Path to a text file of canonical transcript IDs.
                Defaults to "".
            contig_alias (str, optional): Path to a TSV file mapping contig IDs to their alias(es).
                Defaults to "".
            exon_alias (str, optional): Path to a TSV file mapping contig IDs to their alias(es).
                Defaults to "".
            gene_alias (str, optional): Path to a TSV file mapping contig IDs to their alias(es).
                Defaults to "".
            protein_alias (str, optional): Path to a TSV file mapping contig IDs to their alias(es).
                Defaults to "".
            transcript_alias (str, optional): Path to a TSV file mapping contig IDs to their alias(es).
                Defaults to "".
        """
        self.df = read_gtf(gtf, result_type="pandas")  # TODO: switch to 'polars'?
        self.cds_fasta = [read_fasta(i) for i in cds]
        self.dna_fasta = [read_fasta(i) for i in dna]
        self.protein_fasta = [read_fasta(i) for i in peptide]
        self.rna_fasta = [read_fasta(i) for i in rna]
        self._canonical_transcript = txt_to_list(canonical_transcript)
        self._contig_alias = tsv_to_dict(contig_alias)
        self._exon_alias = tsv_to_dict(exon_alias)
        self._gene_alias = tsv_to_dict(gene_alias)
        self._protein_alias = tsv_to_dict(protein_alias)
        self._transcript_alias = tsv_to_dict(transcript_alias)

    # ---------------------------------------------------------------------------------------------
    # all_<feature_symbol>s
    # ---------------------------------------------------------------------------------------------
    def all_contig_ids(self) -> List[str]:
        """List all contig (chromosome) IDs.

        Returns:
            List[str]: Contig IDs
        """
        return self._uniquify_series(self.df[CONTIG_ID])

    def all_exon_ids(self) -> List[str]:
        """List all exon IDs.

        Returns:
            List[str]: Exon IDs
        """
        return self._uniquify_series(self.df[EXON_ID])

    def all_gene_ids(self) -> List[str]:
        """List all gene IDs.

        Returns:
            List[str]: Gene IDs
        """
        return self._uniquify_series(self.df[GENE_ID])

    def all_gene_names(self) -> List[str]:
        """List all gene names.

        Returns:
            List[str]: Gene names
        """
        return self._uniquify_series(self.df[GENE_NAME])

    def all_protein_ids(self) -> List[str]:
        """List all protein IDs.

        Returns:
            List[str]: Protein IDs
        """
        return self._uniquify_series(self.df[PROTEIN_ID])

    def all_transcript_ids(self) -> List[str]:
        """List all transcript IDs.

        Returns:
            List[str]: Transcript IDs
        """
        return self._uniquify_series(self.df[TRANSCRIPT_ID])

    def all_transcript_names(self) -> List[str]:
        """List all transcript names.

        Returns:
            List[str]: Transcript names
        """
        return self._uniquify_series(self.df[TRANSCRIPT_NAME])

    def _uniquify_series(self, series: pd.Series) -> List:
        return sorted(series.dropna().unique().tolist())

    # ---------------------------------------------------------------------------------------------
    # <feature_symbol>s
    # ---------------------------------------------------------------------------------------------
    def contig_ids(self, feature: str) -> List[str]:
        """Given an ID or name, return the corresponding contig IDs.

        Args:
            feature (str): Feature ID or name

        Returns:
            List[str]: Contig IDs
        """
        return self._query_feature(CONTIG_ID, feature)

    def exon_ids(self, feature: str) -> List[str]:
        """Given an ID or name, return the corresponding exon IDs.

        Args:
            feature (str): Feature ID or name

        Returns:
            List[str]: Exon IDs
        """
        return self._query_feature(EXON_ID, feature)

    def gene_ids(self, feature: str) -> List[str]:
        """Given an ID or name, return the corresponding gene IDs.

        Args:
            feature (str): Feature ID or name

        Returns:
            List[str]: Gene IDs
        """
        return self._query_feature(GENE_ID, feature)

    def gene_names(self, feature: str) -> List[str]:
        """Given an ID or name, return the corresponding gene names.

        Args:
            feature (str): Feature ID or name

        Returns:
            List[str]: Gene names
        """
        return self._query_feature(GENE_NAME, feature)

    def protein_ids(self, feature: str) -> List[str]:
        """Given an ID or name, return the corresponding protein IDs.

        Args:
            feature (str): Feature ID or name

        Returns:
            List[str]: Protein IDs
        """
        return self._query_feature(PROTEIN_ID, feature)

    def transcript_ids(self, feature: str) -> List[str]:
        """Given an ID or name, return the corresponding transcript IDs.

        Args:
            feature (str): Feature ID or name

        Returns:
            List[str]: Transcript IDs
        """
        return self._query_feature(TRANSCRIPT_ID, feature)

    def transcript_names(self, feature: str) -> List[str]:
        """Given an ID or name, return the corresponding transcript names.

        Args:
            feature (str): Feature ID or name

        Returns:
            List[str]: Transcript names
        """
        return self._query_feature(TRANSCRIPT_NAME, feature)

    def _query_feature(self, key: str, feature: str) -> List[str]:
        parts = []

        for feature, feature_type in self.normalize_id(feature):
            if feature_type == CONTIG_ID:
                func = self._query_contig_id
            elif feature_type == EXON_ID:
                func = self._query_exon_id
            elif feature_type == GENE_ID:
                func = self._query_gene_id
            elif feature_type == GENE_NAME:
                func = self._query_gene_name
            elif feature_type == PROTEIN_ID:
                func = self._query_protein_id
            elif feature_type == TRANSCRIPT_ID:
                func = self._query_transcript_id
            elif feature_type == TRANSCRIPT_NAME:
                func = self._query_transcript_name
            else:
                raise ValueError(f"Unable to get {key} for {feature} ({feature_type})")

            parts.append(func(feature)[key])

        if parts:
            result = pd.concat(parts)
        else:
            result = pd.Series(dtype="object")

        return self._uniquify_series(result)

    def _query_contig_id(self, feature: Union[List[str], str]) -> pd.DataFrame:
        return self._query(feature, CONTIG_ID)

    def _query_exon_id(self, feature: Union[List[str], str]) -> pd.DataFrame:
        return self._query(feature, EXON_ID)

    def _query_gene_id(self, feature: Union[List[str], str]) -> pd.DataFrame:
        return self._query(feature, GENE_ID)

    def _query_gene_name(self, feature: Union[List[str], str]) -> pd.DataFrame:
        return self._query(feature, GENE_NAME)

    def _query_protein_id(self, feature: Union[List[str], str]) -> pd.DataFrame:
        return self._query(feature, PROTEIN_ID)

    def _query_transcript_id(self, feature: Union[List[str], str]) -> pd.DataFrame:
        return self._query(feature, TRANSCRIPT_ID)

    def _query_transcript_name(self, feature: Union[List[str], str]) -> pd.DataFrame:
        return self._query(feature, TRANSCRIPT_NAME)

    def _query(self, feature: Union[List[str], str], col: str) -> pd.DataFrame:
        feature = [feature] if isinstance(feature, str) else feature
        sudbf = self.df.loc[self.df[col].isin(feature)]

        return sudbf

    # ---------------------------------------------------------------------------------------------
    # normalize_id
    # ---------------------------------------------------------------------------------------------
    def normalize_id(self, feature: str) -> List[Tuple[str, str]]:
        """Normalize an ID or name to the annotated equivalent(s).

        Args:
            feature (str): Feature ID or name

        Returns:
            List[Tuple[str, str]]: List of (normalized ID/name, normalized type (e.g. 'gene',
                'transcript', etc.))
        """
        normalized = []

        feature_type, result = self._normalize_id(feature)
        if feature_type:
            normalized = self._uniquify_series(result[feature_type])

        return [(i, feature_type) for i in normalized]

    def _normalize_id(self, feature: str) -> Tuple[str, pd.DataFrame]:
        feature_type = ""
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
            result = func(feature)
            if not result.empty:
                feature_type = key
                break

        return feature_type, result

    def _normalize_contig_id(self, feature: str) -> pd.DataFrame:
        featurel = [feature] + self.contig_alias(feature)

        return self._query_contig_id(featurel)

    def _normalize_exon_id(self, feature: str) -> pd.DataFrame:
        featurel = [feature] + self.exon_alias(feature)

        return self._query_exon_id(featurel)

    def _normalize_gene_id(self, feature: str) -> pd.DataFrame:
        featurel = [feature] + self.gene_alias(feature)

        return self._query_gene_id(featurel)

    def _normalize_gene_name(self, feature: str) -> pd.DataFrame:
        featurel = [feature] + self.gene_alias(feature)

        return self._query_gene_name(featurel)

    def _normalize_protein_id(self, feature: str) -> pd.DataFrame:
        featurel = [feature] + self.protein_alias(feature)

        return self._query_protein_id(featurel)

    def _normalize_transcript_id(self, feature: str) -> pd.DataFrame:
        featurel = [feature] + self.transcript_alias(feature)

        return self._query_transcript_id(featurel)

    def _normalize_transcript_name(self, feature: str) -> pd.DataFrame:
        featurel = [feature] + self.transcript_alias(feature)

        return self._query_transcript_name(featurel)

    # ---------------------------------------------------------------------------------------------
    # is_<feature>
    # ---------------------------------------------------------------------------------------------
    def is_contig(self, feature: str) -> bool:
        """Check if the given ID or name is a contig.

        Args:
            feature (str): Feature ID or name

        Returns:
            bool: True if is a contig else False
        """
        return any((i[1] == CONTIG_ID for i in self.normalize_id(feature)))

    def is_exon(self, feature: str) -> bool:
        """Check if the given ID or name is an exon.

        Args:
            feature (str): Feature ID or name

        Returns:
            bool: True if is an exon else False
        """
        return any((i[1] == EXON_ID for i in self.normalize_id(feature)))

    def is_gene(self, feature: str) -> bool:
        """Check if the given ID or name is a gene.

        Args:
            feature (str): Feature ID or name

        Returns:
            bool: True if is an gene else False
        """
        return any((i[1] in (GENE_ID, GENE_NAME) for i in self.normalize_id(feature)))

    def is_protein(self, feature: str) -> bool:
        """Check if the given ID or name is a protein.

        Args:
            feature (str): Feature ID or name

        Returns:
            bool: True if is an protein else False
        """
        return any((i[1] == PROTEIN_ID for i in self.normalize_id(feature)))

    def is_transcript(self, feature: str) -> bool:
        """Check if the given ID or name is a transcript.

        Args:
            feature (str): Feature ID or name

        Returns:
            bool: True if is an transcript else False
        """
        return any((i[1] in (TRANSCRIPT_ID, TRANSCRIPT_NAME) for i in self.normalize_id(feature)))

    # ---------------------------------------------------------------------------------------------
    # <feature>_sequence
    # ---------------------------------------------------------------------------------------------
    def cds_sequence(
        self, transcript_id: str, start: Optional[int] = None, end: Optional[int] = None
    ) -> str:
        """Return the CDS sequence between the given position(s), inclusive.

        Args:
            transcript_id (str): transcript ID
            start (Optional[int], optional): Start position of the sequence to return. Defaults to None.
            end (Optional[int], optional): End position of the sequence to return. Defaults to the same as `start`.

        Returns:
            str: The CDS sequence
        """
        return self._sequence(self.cds_fasta, transcript_id, start=start, end=end)

    def dna_sequence(
        self,
        contig_id: str,
        start: Optional[int] = None,
        end: Optional[int] = None,
        strand: str = "+",
    ) -> str:
        """Return the DNA sequence between the given position(s), inclusive.

        Args:
            contig_id (str): contig ID
            start (Optional[int], optional): Start position of the sequence to return. Defaults to None.
            end (Optional[int], optional): End position of the sequence to return. Defaults to the same as `start`.

        Returns:
            str: The DNA sequence
        """
        return self._sequence(self.dna_fasta, contig_id, start=start, end=end, strand=strand)

    def protein_sequence(
        self, protein_id: str, start: Optional[int] = None, end: Optional[int] = None
    ) -> str:
        """Return the protein sequence between the given position(s), inclusive.

        Args:
            protein_id (str): protein ID
            start (Optional[int], optional): Start position of the sequence to return. Defaults to None.
            end (Optional[int], optional): End position of the sequence to return. Defaults to the same as `start`.

        Returns:
            str: The protein sequence
        """
        return self._sequence(self.protein_fasta, protein_id, start=start, end=end)

    def rna_sequence(
        self, transcript_id: str, start: Optional[int] = None, end: Optional[int] = None
    ) -> str:
        """Return the RNA sequence between the given position(s), inclusive.

        Args:
            transcript_id (str): transcript ID
            start (Optional[int], optional): Start position of the sequence to return. Defaults to None.
            end (Optional[int], optional): End position of the sequence to return. Defaults to the same as `start`.

        Returns:
            str: The RNA sequence
        """
        return self._sequence(self.rna_fasta, transcript_id, start=start, end=end)

    def _sequence(
        self,
        fasta: List[Fasta],
        ref: str,
        start: Optional[int] = None,
        end: Optional[int] = None,
        strand: Optional[str] = None,
    ):
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
        if not (0 <= start <= seqlen):
            raise ValueError(f"Start must be from 1 to {seqlen} ({start})")
        if not (0 < end <= seqlen):
            raise ValueError(f"End must be from 1 to {seqlen} ({end})")

        # sanity check that the end position is after the start
        if end < start:
            raise ValueError(f"End must be >= start ({end} < {start})")

        subseq = seq[start - 1 : end]
        if strand == "-":
            subseq = reverse_complement(subseq)

        return subseq

    # ---------------------------------------------------------------------------------------------
    # <feature>_alias
    # ---------------------------------------------------------------------------------------------
    def contig_alias(self, contig_id: str) -> List[str]:
        """List all aliases of the given contig ID.

        Args:
            contig_id (str): contig ID

        Returns:
            List[str]: Any aliases of the given contig ID
        """
        return self._alias(contig_id, self._contig_alias)

    def exon_alias(self, exon_id: str) -> List[str]:
        """List all aliases of the given exon ID.

        Args:
            exon_id (str): exon ID

        Returns:
            List[str]: Any aliases of the given exon ID
        """
        return self._alias(exon_id, self._exon_alias)

    def gene_alias(self, gene_id: str) -> List[str]:
        """List all aliases of the given gene ID.

        Args:
            gene_id (str): gene ID

        Returns:
            List[str]: Any aliases of the given gene ID
        """
        return self._alias(gene_id, self._gene_alias)

    def protein_alias(self, protein_id: str) -> List[str]:
        """List all aliases of the given protein ID.

        Args:
            protein_id (str): protein ID

        Returns:
            List[str]: Any aliases of the given protein ID
        """
        return self._alias(protein_id, self._protein_alias)

    def transcript_alias(self, transcript_id: str) -> List[str]:
        """List all aliases of the given transcript ID.

        Args:
            transcript_id (str): transcript ID

        Returns:
            List[str]: Any aliases of the given transcript ID
        """
        return self._alias(transcript_id, self._transcript_alias)

    def _alias(self, feature: str, alias_dict: Dict[str, List[str]]) -> List[str]:
        for key in [feature, strip_version(feature)]:
            if alias := alias_dict.get(key, []):
                return alias
        else:
            return []

    # ---------------------------------------------------------------------------------------------
    # is_canonical_transcript
    # ---------------------------------------------------------------------------------------------
    def is_canonical_transcript(self, transcript_id: str) -> bool:
        """Check if the given transcript ID is the canonical transcript for its gene.

        Args:
            transcript_id (str): transcript ID

        Returns:
            bool: True if the transcript is the canonical transcript else False
        """
        return transcript_id in self._canonical_transcript

    # ---------------------------------------------------------------------------------------------
    # <feature>
    # ---------------------------------------------------------------------------------------------
    def cdna(self, feature: str, canonical: bool = False) -> List[CdnaMappablePosition]:
        """Return the cDNA position(s) matching the given feature ID or name.

        Args:
            feature (str): Feature ID or name
            canonical (bool, optional): Only map to the canonical transcript. Defaults to False.

        Returns:
            List[CdnaMappablePosition]: One or more cDNA positions.
        """
        result = []

        transcript_ids = self.transcript_ids(feature)
        mask = (self.df[TRANSCRIPT_ID].isin(transcript_ids)) & (self.df["feature"] == "cdna")
        for _, cdna in self.df[mask].iterrows():
            if canonical and not self.is_canonical_transcript(cdna.transcript_id):
                continue

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

    def dna(self, feature: str) -> List[DnaMappablePosition]:
        """Return the DNA position(s) matching the given feature ID or name.

        Args:
            feature (str): Feature ID or name
            canonical (bool, optional): Only map to the canonical transcript. Defaults to False.

        Returns:
            List[DnaMappablePosition]: One or more DNA positions.
        """
        result = []

        # get the strand of the original feature
        _, df = self._normalize_id(feature)
        strand_list = self._uniquify_series(df["strand"])

        for contig_id in self.contig_ids(feature):
            for fasta in self.dna_fasta:
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

    def exon(self, feature: str, canonical: bool = False) -> List[ExonMappablePosition]:
        """Return the exon position(s) matching the given feature ID or name.

        Args:
            feature (str): Feature ID or name
            canonical (bool, optional): Only map to the canonical transcript. Defaults to False.

        Returns:
            List[ExonMappablePosition]: One or more exon positions.
        """
        result = []

        exon_ids = self.exon_ids(feature)
        mask = (self.df[EXON_ID].isin(exon_ids)) & (self.df["feature"] == "exon")
        for _, exon in self.df[mask].iterrows():
            if canonical and not self.is_canonical_transcript(exon.transcript_id):
                continue

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

    def gene(self, feature: str) -> List[DnaMappablePosition]:
        """Return the gene position(s) matching the given feature ID or name.

        Args:
            feature (str): Feature ID or name
            canonical (bool, optional): Only map to the canonical transcript. Defaults to False.

        Returns:
            List[DnaMappablePosition]: One or more DNA positions.
        """
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

    def protein(self, feature: str, canonical: bool = False) -> List[ProteinMappablePosition]:
        """Return the protein position(s) matching the given feature ID or name.

        Args:
            feature (str): Feature ID or name
            canonical (bool, optional): Only map to the canonical transcript. Defaults to False.

        Returns:
            List[ProteinMappablePosition]: One or more protein positions.
        """
        result = []
        for cdna in self.cdna(feature, canonical=canonical):
            # Convert the cDNA position to a protein position
            protein_start = calc_cdna_to_protein(cdna.start)
            protein_end = calc_cdna_to_protein(cdna.end)
            result.append(
                ProteinMappablePosition.copy_from(
                    cdna, start=protein_start, start_offset=0, end=protein_end, end_offset=0
                )
            )

        return sorted(set(result))

    def rna(self, feature: str, canonical: bool = False) -> List[RnaMappablePosition]:
        """Return the RNA position(s) matching the given feature ID or name.

        Args:
            feature (str): Feature ID or name
            canonical (bool, optional): Only map to the canonical transcript. Defaults to False.

        Returns:
            List[RnaMappablePosition]: One or more RNA positions.
        """
        result = []

        transcript_ids = self.transcript_ids(feature)
        mask = (self.df[TRANSCRIPT_ID].isin(transcript_ids)) & (self.df["feature"] == "transcript")
        for _, transcript in self.df[mask].iterrows():
            if canonical and not self.is_canonical_transcript(transcript.transcript_id):
                continue

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
        refseq: str = "",
        altseq: str = "",
        canonical: bool = False,
    ) -> Union[List[CdnaMappablePosition], List[CdnaSmallVariant]]:
        """Map a cDNA position to zero or more cDNA positions.

        Args:
            feature (str): Feature ID or name
            start (int): Start position
            end (Optional[int], optional): End position. Defaults to the same as `start`.
            strand (Optional[str], optional): Strand ('+' or '-'). Defaults to checking either.
            start_offset (Optional[int], optional): Nucleotide offset from `start`. Defaults to 0.
            end_offset (Optional[int], optional): Nucleotide offset from `end`. Defaults to 0.
            refseq (str, optional): Reference allele. Defaults to "".
            altseq (str, optional): Alternate allele. Defaults to "".
            canonical (bool, optional): Only map to the canonical transcript. Defaults to False.

        Returns:
            Union[List[CdnaMappablePosition], List[CdnaSmallVariant]]:
                If `refseq` or `altseq` was given, returns a variant. Otherwise returns a position.
        """
        if refseq or altseq:
            return self._map_variant(
                self.transcript_ids,
                self._cdna_to_cdna_variant,
                feature,
                start,
                refseq,
                altseq,
                start_offset,
                end,
                end_offset,
                strand,
                canonical,
            )
        else:
            return self._map(
                self.transcript_ids,
                self._cdna_to_cdna,
                feature,
                start,
                start_offset,
                end,
                end_offset,
                strand,
                canonical,
            )

    def cdna_to_dna(
        self,
        feature: str,
        start: int,
        end: Optional[int] = None,
        strand: Optional[str] = None,
        start_offset: Optional[int] = None,
        end_offset: Optional[int] = None,
        refseq: str = "",
        altseq: str = "",
        canonical: bool = False,
    ) -> Union[List[DnaMappablePosition], List[DnaSmallVariant]]:
        """Map a cDNA position to zero or more DNA positions.

        Args:
            feature (str): Feature ID or name
            start (int): Start position
            end (Optional[int], optional): End position. Defaults to the same as `start`.
            strand (Optional[str], optional): Strand ('+' or '-'). Defaults to checking either.
            start_offset (Optional[int], optional): Nucleotide offset from `start`. Defaults to 0.
            end_offset (Optional[int], optional): Nucleotide offset from `end`. Defaults to 0.
            refseq (str, optional): Reference allele. Defaults to "".
            altseq (str, optional): Alternate allele. Defaults to "".
            canonical (bool, optional): Only map to the canonical transcript. Defaults to False.

        Returns:
            Union[List[DnaMappablePosition], List[DnaSmallVariant]]:
                If `refseq` or `altseq` was given, returns a variant. Otherwise returns a position.
        """
        if refseq or altseq:
            return self._map_variant(
                self.transcript_ids,
                self._cdna_to_dna_variant,
                feature,
                start,
                refseq,
                altseq,
                start_offset,
                end,
                end_offset,
                strand,
                canonical,
            )
        else:
            return self._map(
                self.transcript_ids,
                self._cdna_to_dna,
                feature,
                start,
                start_offset,
                end,
                end_offset,
                strand,
                canonical,
            )

    def cdna_to_exon(
        self,
        feature: str,
        start: int,
        end: Optional[int] = None,
        strand: Optional[str] = None,
        start_offset: Optional[int] = None,
        end_offset: Optional[int] = None,
        refseq: str = "",
        altseq: str = "",
        canonical: bool = False,
    ) -> Union[List[ExonMappablePosition], List[ExonSmallVariant]]:
        """Map a cDNA position to zero or more exon positions.

        Args:
            feature (str): Feature ID or name
            start (int): Start position
            end (Optional[int], optional): End position. Defaults to the same as `start`.
            strand (Optional[str], optional): Strand ('+' or '-'). Defaults to checking either.
            start_offset (Optional[int], optional): Nucleotide offset from `start`. Defaults to 0.
            end_offset (Optional[int], optional): Nucleotide offset from `end`. Defaults to 0.
            refseq (str, optional): Reference allele. Defaults to "".
            altseq (str, optional): Alternate allele. Defaults to "".
            canonical (bool, optional): Only map to the canonical transcript. Defaults to False.

        Returns:
            Union[List[ExonMappablePosition], List[ExonSmallVariant]]:
                If `refseq` or `altseq` was given, returns a variant. Otherwise returns a position.
        """
        if refseq or altseq:
            return self._map_variant(
                self.transcript_ids,
                self._cdna_to_exon_variant,
                feature,
                start,
                refseq,
                altseq,
                start_offset,
                end,
                end_offset,
                strand,
                canonical,
            )
        else:
            return self._map(
                self.transcript_ids,
                self._cdna_to_exon,
                feature,
                start,
                start_offset,
                end,
                end_offset,
                strand,
                canonical,
            )

    def cdna_to_protein(
        self,
        feature: str,
        start: int,
        end: Optional[int] = None,
        strand: Optional[str] = None,
        start_offset: Optional[int] = None,
        end_offset: Optional[int] = None,
        refseq: str = "",
        altseq: str = "",
        canonical: bool = False,
    ) -> Union[List[ProteinMappablePosition], List[ProteinSmallVariant]]:
        """Map a cDNA position to zero or more protein positions.

        Args:
            feature (str): Feature ID or name
            start (int): Start position
            end (Optional[int], optional): End position. Defaults to the same as `start`.
            strand (Optional[str], optional): Strand ('+' or '-'). Defaults to checking either.
            start_offset (Optional[int], optional): Nucleotide offset from `start`. Defaults to 0.
            end_offset (Optional[int], optional): Nucleotide offset from `end`. Defaults to 0.
            refseq (str, optional): Reference allele. Defaults to "".
            altseq (str, optional): Alternate allele. Defaults to "".
            canonical (bool, optional): Only map to the canonical transcript. Defaults to False.

        Returns:
            Union[List[ProteinMappablePosition], List[ProteinSmallVariant]]:
                If `refseq` or `altseq` was given, returns a variant. Otherwise returns a position.
        """
        if refseq or altseq:
            return self._map_variant(
                self.transcript_ids,
                self._cdna_to_protein_variant,
                feature,
                start,
                refseq,
                altseq,
                start_offset,
                end,
                end_offset,
                strand,
                canonical,
            )
        else:
            return self._map(
                self.transcript_ids,
                self._cdna_to_protein,
                feature,
                start,
                start_offset,
                end,
                end_offset,
                strand,
                canonical,
            )

    def cdna_to_rna(
        self,
        feature: str,
        start: int,
        end: Optional[int] = None,
        strand: Optional[str] = None,
        start_offset: Optional[int] = None,
        end_offset: Optional[int] = None,
        refseq: str = "",
        altseq: str = "",
        canonical: bool = False,
    ) -> Union[List[RnaMappablePosition], List[RnaSmallVariant]]:
        """Map a cDNA position to zero or more RNA positions.

        Args:
            feature (str): Feature ID or name
            start (int): Start position
            end (Optional[int], optional): End position. Defaults to the same as `start`.
            strand (Optional[str], optional): Strand ('+' or '-'). Defaults to checking either.
            start_offset (Optional[int], optional): Nucleotide offset from `start`. Defaults to 0.
            end_offset (Optional[int], optional): Nucleotide offset from `end`. Defaults to 0.
            refseq (str, optional): Reference allele. Defaults to "".
            altseq (str, optional): Alternate allele. Defaults to "".
            canonical (bool, optional): Only map to the canonical transcript. Defaults to False.

        Returns:
            Union[List[RnaMappablePosition], List[RnaSmallVariant]]:
                If `refseq` or `altseq` was given, returns a variant. Otherwise returns a position.
        """
        if refseq or altseq:
            return self._map_variant(
                self.transcript_ids,
                self._cdna_to_rna_variant,
                feature,
                start,
                refseq,
                altseq,
                start_offset,
                end,
                end_offset,
                strand,
                canonical,
            )
        else:
            return self._map(
                self.transcript_ids,
                self._cdna_to_rna,
                feature,
                start,
                start_offset,
                end,
                end_offset,
                strand,
                canonical,
            )

    def dna_to_cdna(
        self,
        feature: str,
        start: int,
        end: Optional[int] = None,
        strand: Optional[str] = None,
        start_offset: Optional[int] = None,
        end_offset: Optional[int] = None,
        refseq: str = "",
        altseq: str = "",
        canonical: bool = False,
    ) -> Union[List[CdnaMappablePosition], List[CdnaSmallVariant]]:
        """Map a DNA position to zero or more cDNA positions.

        Args:
            feature (str): Feature ID or name
            start (int): Start position
            end (Optional[int], optional): End position. Defaults to the same as `start`.
            strand (Optional[str], optional): Strand ('+' or '-'). Defaults to checking either.
            start_offset (Optional[int], optional): Nucleotide offset from `start`. Defaults to 0.
            end_offset (Optional[int], optional): Nucleotide offset from `end`. Defaults to 0.
            refseq (str, optional): Reference allele. Defaults to "".
            altseq (str, optional): Alternate allele. Defaults to "".
            canonical (bool, optional): Only map to the canonical transcript. Defaults to False.

        Returns:
            Union[List[CdnaMappablePosition], List[CdnaSmallVariant]]:
                If `refseq` or `altseq` was given, returns a variant. Otherwise returns a position.
        """
        if refseq or altseq:
            return self._map_variant(
                self.contig_ids,
                self._dna_to_cdna_variant,
                feature,
                start,
                refseq,
                altseq,
                start_offset,
                end,
                end_offset,
                strand,
                canonical,
            )
        else:
            return self._map(
                self.contig_ids,
                self._dna_to_cdna,
                feature,
                start,
                start_offset,
                end,
                end_offset,
                strand,
                canonical,
            )

    def dna_to_dna(
        self,
        feature: str,
        start: int,
        end: Optional[int] = None,
        strand: Optional[str] = None,
        start_offset: Optional[int] = None,
        end_offset: Optional[int] = None,
        refseq: str = "",
        altseq: str = "",
        canonical: bool = False,
    ) -> Union[List[DnaMappablePosition], List[DnaSmallVariant]]:
        """Map a DNA position to zero or more DNA positions.

        Args:
            feature (str): Feature ID or name
            start (int): Start position
            end (Optional[int], optional): End position. Defaults to the same as `start`.
            strand (Optional[str], optional): Strand ('+' or '-'). Defaults to checking either.
            start_offset (Optional[int], optional): Nucleotide offset from `start`. Defaults to 0.
            end_offset (Optional[int], optional): Nucleotide offset from `end`. Defaults to 0.
            refseq (str, optional): Reference allele. Defaults to "".
            altseq (str, optional): Alternate allele. Defaults to "".
            canonical (bool, optional): Only map to the canonical transcript. Defaults to False.

        Returns:
            Union[List[DnaMappablePosition], List[DnaSmallVariant]]:
                If `refseq` or `altseq` was given, returns a variant. Otherwise returns a position.
        """
        if refseq or altseq:
            return self._map_variant(
                self.contig_ids,
                self._dna_to_dna_variant,
                feature,
                start,
                refseq,
                altseq,
                start_offset,
                end,
                end_offset,
                strand,
                canonical,
            )
        else:
            return self._map(
                self.contig_ids,
                self._dna_to_dna,
                feature,
                start,
                start_offset,
                end,
                end_offset,
                strand,
                canonical,
            )

    def dna_to_exon(
        self,
        feature: str,
        start: int,
        end: Optional[int] = None,
        strand: Optional[str] = None,
        start_offset: Optional[int] = None,
        end_offset: Optional[int] = None,
        refseq: str = "",
        altseq: str = "",
        canonical: bool = False,
    ) -> Union[List[ExonMappablePosition], List[ExonSmallVariant]]:
        """Map a DNA position to zero or more exon positions.

        Args:
            feature (str): Feature ID or name
            start (int): Start position
            end (Optional[int], optional): End position. Defaults to the same as `start`.
            strand (Optional[str], optional): Strand ('+' or '-'). Defaults to checking either.
            start_offset (Optional[int], optional): Nucleotide offset from `start`. Defaults to 0.
            end_offset (Optional[int], optional): Nucleotide offset from `end`. Defaults to 0.
            refseq (str, optional): Reference allele. Defaults to "".
            altseq (str, optional): Alternate allele. Defaults to "".
            canonical (bool, optional): Only map to the canonical transcript. Defaults to False.

        Returns:
            Union[List[ExonMappablePosition], List[ExonSmallVariant]]:
                If `refseq` or `altseq` was given, returns a variant. Otherwise returns a position.
        """
        if refseq or altseq:
            return self._map_variant(
                self.contig_ids,
                self._dna_to_exon_variant,
                feature,
                start,
                refseq,
                altseq,
                start_offset,
                end,
                end_offset,
                strand,
                canonical,
            )
        else:
            return self._map(
                self.contig_ids,
                self._dna_to_exon,
                feature,
                start,
                start_offset,
                end,
                end_offset,
                strand,
                canonical,
            )

    def dna_to_protein(
        self,
        feature: str,
        start: int,
        end: Optional[int] = None,
        strand: Optional[str] = None,
        start_offset: Optional[int] = None,
        end_offset: Optional[int] = None,
        refseq: str = "",
        altseq: str = "",
        canonical: bool = False,
    ) -> Union[List[ProteinMappablePosition], List[ProteinSmallVariant]]:
        """Map a DNA position to zero or more protein positions.

        Args:
            feature (str): Feature ID or name
            start (int): Start position
            end (Optional[int], optional): End position. Defaults to the same as `start`.
            strand (Optional[str], optional): Strand ('+' or '-'). Defaults to checking either.
            start_offset (Optional[int], optional): Nucleotide offset from `start`. Defaults to 0.
            end_offset (Optional[int], optional): Nucleotide offset from `end`. Defaults to 0.
            refseq (str, optional): Reference allele. Defaults to "".
            altseq (str, optional): Alternate allele. Defaults to "".
            canonical (bool, optional): Only map to the canonical transcript. Defaults to False.

        Returns:
            Union[List[ProteinMappablePosition], List[ProteinSmallVariant]]:
                If `refseq` or `altseq` was given, returns a variant. Otherwise returns a position.
        """
        if refseq or altseq:
            return self._map_variant(
                self.contig_ids,
                self._dna_to_protein_variant,
                feature,
                start,
                refseq,
                altseq,
                start_offset,
                end,
                end_offset,
                strand,
                canonical,
            )
        else:
            return self._map(
                self.contig_ids,
                self._dna_to_protein,
                feature,
                start,
                start_offset,
                end,
                end_offset,
                strand,
                canonical,
            )

    def dna_to_rna(
        self,
        feature: str,
        start: int,
        end: Optional[int] = None,
        strand: Optional[str] = None,
        start_offset: Optional[int] = None,
        end_offset: Optional[int] = None,
        refseq: str = "",
        altseq: str = "",
        canonical: bool = False,
    ) -> Union[List[RnaMappablePosition], List[RnaSmallVariant]]:
        """Map a DNA position to zero or more RNA positions.

        Args:
            feature (str): Feature ID or name
            start (int): Start position
            end (Optional[int], optional): End position. Defaults to the same as `start`.
            strand (Optional[str], optional): Strand ('+' or '-'). Defaults to checking either.
            start_offset (Optional[int], optional): Nucleotide offset from `start`. Defaults to 0.
            end_offset (Optional[int], optional): Nucleotide offset from `end`. Defaults to 0.
            refseq (str, optional): Reference allele. Defaults to "".
            altseq (str, optional): Alternate allele. Defaults to "".
            canonical (bool, optional): Only map to the canonical transcript. Defaults to False.

        Returns:
            Union[List[RnaMappablePosition], List[RnaSmallVariant]]:
                If `refseq` or `altseq` was given, returns a variant. Otherwise returns a position.
        """
        if refseq or altseq:
            return self._map_variant(
                self.contig_ids,
                self._dna_to_rna_variant,
                feature,
                start,
                refseq,
                altseq,
                start_offset,
                end,
                end_offset,
                strand,
                canonical,
            )
        else:
            return self._map(
                self.contig_ids,
                self._dna_to_rna,
                feature,
                start,
                start_offset,
                end,
                end_offset,
                strand,
                canonical,
            )

    def exon_to_cdna(
        self,
        feature: str,
        start: Optional[int] = None,
        end: Optional[int] = None,
        strand: Optional[str] = None,
        start_offset: Optional[int] = None,
        end_offset: Optional[int] = None,
        canonical: bool = False,
    ) -> List[CdnaMappablePosition]:
        """Map an exon position to zero or more cDNA positions.

        Args:
            feature (str): Feature ID or name
            start (int): Start position
            end (Optional[int], optional): End position. Defaults to the same as `start`.
            strand (Optional[str], optional): Strand ('+' or '-'). Defaults to checking either.
            start_offset (Optional[int], optional): Nucleotide offset from `start`. Defaults to 0.
            end_offset (Optional[int], optional): Nucleotide offset from `end`. Defaults to 0.
            canonical (bool, optional): Only map to the canonical transcript. Defaults to False.

        Returns:
            List[CdnaMappablePosition]: Zero or more cDNA positions
        """
        return self._map_exon(
            self._exon_to_cdna, feature, start, start_offset, end, end_offset, strand, canonical
        )

    def exon_to_dna(
        self,
        feature: str,
        start: Optional[int] = None,
        end: Optional[int] = None,
        strand: Optional[str] = None,
        start_offset: Optional[int] = None,
        end_offset: Optional[int] = None,
        canonical: bool = False,
    ) -> List[DnaMappablePosition]:
        """Map an exon position to zero or more DNA positions.

        Args:
            feature (str): Feature ID or name
            start (int): Start position
            end (Optional[int], optional): End position. Defaults to the same as `start`.
            strand (Optional[str], optional): Strand ('+' or '-'). Defaults to checking either.
            start_offset (Optional[int], optional): Nucleotide offset from `start`. Defaults to 0.
            end_offset (Optional[int], optional): Nucleotide offset from `end`. Defaults to 0.
            canonical (bool, optional): Only map to the canonical transcript. Defaults to False.

        Returns:
            List[DnaMappablePosition]: Zero or more DNA positions
        """
        return self._map_exon(
            self._exon_to_dna, feature, start, start_offset, end, end_offset, strand, canonical
        )

    def exon_to_exon(
        self,
        feature: str,
        start: Optional[int] = None,
        end: Optional[int] = None,
        strand: Optional[str] = None,
        start_offset: Optional[int] = None,
        end_offset: Optional[int] = None,
        canonical: bool = False,
    ) -> List[ExonMappablePosition]:
        """Map an exon position to zero or more exon positions.

        Args:
            feature (str): Feature ID or name
            start (int): Start position
            end (Optional[int], optional): End position. Defaults to the same as `start`.
            strand (Optional[str], optional): Strand ('+' or '-'). Defaults to checking either.
            start_offset (Optional[int], optional): Nucleotide offset from `start`. Defaults to 0.
            end_offset (Optional[int], optional): Nucleotide offset from `end`. Defaults to 0.
            canonical (bool, optional): Only map to the canonical transcript. Defaults to False.

        Returns:
            List[ExonMappablePosition]: Zero or more exon positions
        """
        return self._map_exon(
            self._exon_to_exon, feature, start, start_offset, end, end_offset, strand, canonical
        )

    def exon_to_protein(
        self,
        feature: str,
        start: Optional[int] = None,
        end: Optional[int] = None,
        strand: Optional[str] = None,
        start_offset: Optional[int] = None,
        end_offset: Optional[int] = None,
        canonical: bool = False,
    ) -> List[ProteinMappablePosition]:
        """Map an exon position to zero or more protein positions.

        Args:
            feature (str): Feature ID or name
            start (int): Start position
            end (Optional[int], optional): End position. Defaults to the same as `start`.
            strand (Optional[str], optional): Strand ('+' or '-'). Defaults to checking either.
            start_offset (Optional[int], optional): Nucleotide offset from `start`. Defaults to 0.
            end_offset (Optional[int], optional): Nucleotide offset from `end`. Defaults to 0.
            canonical (bool, optional): Only map to the canonical transcript. Defaults to False.

        Returns:
            List[ProteinMappablePosition]: Zero or more protein positions
        """
        return self._map_exon(
            self._exon_to_protein, feature, start, start_offset, end, end_offset, strand, canonical
        )

    def exon_to_rna(
        self,
        feature: str,
        start: Optional[int] = None,
        end: Optional[int] = None,
        strand: Optional[str] = None,
        start_offset: Optional[int] = None,
        end_offset: Optional[int] = None,
        canonical: bool = False,
    ) -> List[RnaMappablePosition]:
        """Map an exon position to zero or more RNA positions.

        Args:
            feature (str): Feature ID or name
            start (int): Start position
            end (Optional[int], optional): End position. Defaults to the same as `start`.
            strand (Optional[str], optional): Strand ('+' or '-'). Defaults to checking either.
            start_offset (Optional[int], optional): Nucleotide offset from `start`. Defaults to 0.
            end_offset (Optional[int], optional): Nucleotide offset from `end`. Defaults to 0.
            canonical (bool, optional): Only map to the canonical transcript. Defaults to False.

        Returns:
            List[RnaMappablePosition]: Zero or more RNA positions
        """
        return self._map_exon(
            self._exon_to_rna, feature, start, start_offset, end, end_offset, strand, canonical
        )

    def protein_to_cdna(
        self,
        feature: str,
        start: int,
        end: Optional[int] = None,
        strand: Optional[str] = None,
        start_offset: Optional[int] = None,
        end_offset: Optional[int] = None,
        refseq: str = "",
        altseq: str = "",
        canonical: bool = False,
    ) -> Union[List[CdnaMappablePosition], List[CdnaSmallVariant]]:
        """Map a protein position to zero or more cDNA positions.

        Args:
            feature (str): Feature ID or name
            start (int): Start position
            end (Optional[int], optional): End position. Defaults to the same as `start`.
            strand (Optional[str], optional): Strand ('+' or '-'). Defaults to checking either.
            start_offset (Optional[int], optional): Nucleotide offset from `start`. Defaults to 0.
            end_offset (Optional[int], optional): Nucleotide offset from `end`. Defaults to 0.
            refseq (str, optional): Reference allele. Defaults to "".
            altseq (str, optional): Alternate allele. Defaults to "".
            canonical (bool, optional): Only map to the canonical transcript. Defaults to False.

        Returns:
            Union[List[CdnaMappablePosition], List[CdnaSmallVariant]]:
                If `refseq` or `altseq` was given, returns a variant. Otherwise returns a position.
        """
        if refseq or altseq:
            return self._map_variant(
                self.transcript_ids,
                self._protein_to_cdna_variant,
                feature,
                start,
                refseq,
                altseq,
                start_offset,
                end,
                end_offset,
                strand,
                canonical,
            )
        else:
            return self._map(
                self.transcript_ids,
                self._protein_to_cdna,
                feature,
                start,
                start_offset,
                end,
                end_offset,
                strand,
                canonical,
            )

    def protein_to_dna(
        self,
        feature: str,
        start: int,
        end: Optional[int] = None,
        strand: Optional[str] = None,
        start_offset: Optional[int] = None,
        end_offset: Optional[int] = None,
        refseq: str = "",
        altseq: str = "",
        canonical: bool = False,
    ) -> Union[List[DnaMappablePosition], List[DnaSmallVariant]]:
        """Map a protein position to zero or more DNA positions.

        Args:
            feature (str): Feature ID or name
            start (int): Start position
            end (Optional[int], optional): End position. Defaults to the same as `start`.
            strand (Optional[str], optional): Strand ('+' or '-'). Defaults to checking either.
            start_offset (Optional[int], optional): Nucleotide offset from `start`. Defaults to 0.
            end_offset (Optional[int], optional): Nucleotide offset from `end`. Defaults to 0.
            refseq (str, optional): Reference allele. Defaults to "".
            altseq (str, optional): Alternate allele. Defaults to "".
            canonical (bool, optional): Only map to the canonical transcript. Defaults to False.

        Returns:
            Union[List[DnaMappablePosition], List[DnaSmallVariant]]:
                If `refseq` or `altseq` was given, returns a variant. Otherwise returns a position.
        """
        if refseq or altseq:
            return self._map_variant(
                self.transcript_ids,
                self._protein_to_dna_variant,
                feature,
                start,
                refseq,
                altseq,
                start_offset,
                end,
                end_offset,
                strand,
                canonical,
            )
        else:
            return self._map(
                self.transcript_ids,
                self._protein_to_dna,
                feature,
                start,
                start_offset,
                end,
                end_offset,
                strand,
                canonical,
            )

    def protein_to_exon(
        self,
        feature: str,
        start: int,
        end: Optional[int] = None,
        strand: Optional[str] = None,
        start_offset: Optional[int] = None,
        end_offset: Optional[int] = None,
        refseq: str = "",
        altseq: str = "",
        canonical: bool = False,
    ) -> Union[List[ExonMappablePosition], List[ExonSmallVariant]]:
        """Map a protein position to zero or more exon positions.

        Args:
            feature (str): Feature ID or name
            start (int): Start position
            end (Optional[int], optional): End position. Defaults to the same as `start`.
            strand (Optional[str], optional): Strand ('+' or '-'). Defaults to checking either.
            start_offset (Optional[int], optional): Nucleotide offset from `start`. Defaults to 0.
            end_offset (Optional[int], optional): Nucleotide offset from `end`. Defaults to 0.
            refseq (str, optional): Reference allele. Defaults to "".
            altseq (str, optional): Alternate allele. Defaults to "".
            canonical (bool, optional): Only map to the canonical transcript. Defaults to False.

        Returns:
            Union[List[ExonMappablePosition], List[ExonSmallVariant]]:
                If `refseq` or `altseq` was given, returns a variant. Otherwise returns a position.
        """
        if refseq or altseq:
            return self._map_variant(
                self.transcript_ids,
                self._protein_to_exon_variant,
                feature,
                start,
                refseq,
                altseq,
                start_offset,
                end,
                end_offset,
                strand,
                canonical,
            )
        else:
            return self._map(
                self.transcript_ids,
                self._protein_to_exon,
                feature,
                start,
                start_offset,
                end,
                end_offset,
                strand,
                canonical,
            )

    def protein_to_protein(
        self,
        feature: str,
        start: int,
        end: Optional[int] = None,
        strand: Optional[str] = None,
        start_offset: Optional[int] = None,
        end_offset: Optional[int] = None,
        refseq: str = "",
        altseq: str = "",
        canonical: bool = False,
    ) -> Union[List[ProteinMappablePosition], List[ProteinSmallVariant]]:
        """Map a protein position to zero or more protein positions.

        Args:
            feature (str): Feature ID or name
            start (int): Start position
            end (Optional[int], optional): End position. Defaults to the same as `start`.
            strand (Optional[str], optional): Strand ('+' or '-'). Defaults to checking either.
            start_offset (Optional[int], optional): Nucleotide offset from `start`. Defaults to 0.
            end_offset (Optional[int], optional): Nucleotide offset from `end`. Defaults to 0.
            refseq (str, optional): Reference allele. Defaults to "".
            altseq (str, optional): Alternate allele. Defaults to "".
            canonical (bool, optional): Only map to the canonical transcript. Defaults to False.

        Returns:
            Union[List[ProteinMappablePosition], List[ProteinSmallVariant]]:
                If `refseq` or `altseq` was given, returns a variant. Otherwise returns a position.
        """
        if refseq or altseq:
            return self._map_variant(
                self.transcript_ids,
                self._protein_to_protein_variant,
                feature,
                start,
                refseq,
                altseq,
                start_offset,
                end,
                end_offset,
                strand,
                canonical,
            )
        else:
            return self._map(
                self.transcript_ids,
                self._protein_to_protein,
                feature,
                start,
                start_offset,
                end,
                end_offset,
                strand,
                canonical,
            )

    def protein_to_rna(
        self,
        feature: str,
        start: int,
        end: Optional[int] = None,
        strand: Optional[str] = None,
        start_offset: Optional[int] = None,
        end_offset: Optional[int] = None,
        refseq: str = "",
        altseq: str = "",
        canonical: bool = False,
    ) -> Union[List[RnaMappablePosition], List[RnaSmallVariant]]:
        """Map a protein position to zero or more RNA positions.

        Args:
            feature (str): Feature ID or name
            start (int): Start position
            end (Optional[int], optional): End position. Defaults to the same as `start`.
            strand (Optional[str], optional): Strand ('+' or '-'). Defaults to checking either.
            start_offset (Optional[int], optional): Nucleotide offset from `start`. Defaults to 0.
            end_offset (Optional[int], optional): Nucleotide offset from `end`. Defaults to 0.
            refseq (str, optional): Reference allele. Defaults to "".
            altseq (str, optional): Alternate allele. Defaults to "".
            canonical (bool, optional): Only map to the canonical transcript. Defaults to False.

        Returns:
            Union[List[RnaMappablePosition], List[RnaSmallVariant]]:
                If `refseq` or `altseq` was given, returns a variant. Otherwise returns a position.
        """
        if refseq or altseq:
            return self._map_variant(
                self.transcript_ids,
                self._protein_to_rna_variant,
                feature,
                start,
                refseq,
                altseq,
                start_offset,
                end,
                end_offset,
                strand,
                canonical,
            )
        else:
            return self._map(
                self.transcript_ids,
                self._protein_to_rna,
                feature,
                start,
                start_offset,
                end,
                end_offset,
                strand,
                canonical,
            )

    def rna_to_cdna(
        self,
        feature: str,
        start: int,
        end: Optional[int] = None,
        strand: Optional[str] = None,
        start_offset: Optional[int] = None,
        end_offset: Optional[int] = None,
        refseq: str = "",
        altseq: str = "",
        canonical: bool = False,
    ) -> Union[List[CdnaMappablePosition], List[CdnaSmallVariant]]:
        """Map an RNA position to zero or more cDNA positions.

        Args:
            feature (str): Feature ID or name
            start (int): Start position
            end (Optional[int], optional): End position. Defaults to the same as `start`.
            strand (Optional[str], optional): Strand ('+' or '-'). Defaults to checking either.
            start_offset (Optional[int], optional): Nucleotide offset from `start`. Defaults to 0.
            end_offset (Optional[int], optional): Nucleotide offset from `end`. Defaults to 0.
            refseq (str, optional): Reference allele. Defaults to "".
            altseq (str, optional): Alternate allele. Defaults to "".
            canonical (bool, optional): Only map to the canonical transcript. Defaults to False.

        Returns:
            Union[List[CdnaMappablePosition], List[CdnaSmallVariant]]:
                If `refseq` or `altseq` was given, returns a variant. Otherwise returns a position.
        """
        if refseq or altseq:
            return self._map_variant(
                self.transcript_ids,
                self._rna_to_cdna_variant,
                feature,
                start,
                refseq,
                altseq,
                start_offset,
                end,
                end_offset,
                strand,
                canonical,
            )
        else:
            return self._map(
                self.transcript_ids,
                self._rna_to_cdna,
                feature,
                start,
                start_offset,
                end,
                end_offset,
                strand,
                canonical,
            )

    def rna_to_dna(
        self,
        feature: str,
        start: int,
        end: Optional[int] = None,
        strand: Optional[str] = None,
        start_offset: Optional[int] = None,
        end_offset: Optional[int] = None,
        refseq: str = "",
        altseq: str = "",
        canonical: bool = False,
    ) -> Union[List[DnaMappablePosition], List[DnaSmallVariant]]:
        """Map an RNA position to zero or more DNA positions.

        Args:
            feature (str): Feature ID or name
            start (int): Start position
            end (Optional[int], optional): End position. Defaults to the same as `start`.
            strand (Optional[str], optional): Strand ('+' or '-'). Defaults to checking either.
            start_offset (Optional[int], optional): Nucleotide offset from `start`. Defaults to 0.
            end_offset (Optional[int], optional): Nucleotide offset from `end`. Defaults to 0.
            refseq (str, optional): Reference allele. Defaults to "".
            altseq (str, optional): Alternate allele. Defaults to "".
            canonical (bool, optional): Only map to the canonical transcript. Defaults to False.

        Returns:
            Union[List[DnaMappablePosition], List[DnaSmallVariant]]:
                If `refseq` or `altseq` was given, returns a variant. Otherwise returns a position.
        """
        if refseq or altseq:
            return self._map_variant(
                self.transcript_ids,
                self._rna_to_dna_variant,
                feature,
                start,
                refseq,
                altseq,
                start_offset,
                end,
                end_offset,
                strand,
                canonical,
            )
        else:
            return self._map(
                self.transcript_ids,
                self._rna_to_dna,
                feature,
                start,
                start_offset,
                end,
                end_offset,
                strand,
                canonical,
            )

    def rna_to_exon(
        self,
        feature: str,
        start: int,
        end: Optional[int] = None,
        strand: Optional[str] = None,
        start_offset: Optional[int] = None,
        end_offset: Optional[int] = None,
        refseq: str = "",
        altseq: str = "",
        canonical: bool = False,
    ) -> Union[List[ExonMappablePosition], List[ExonSmallVariant]]:
        """Map an RNA position to zero or more exon positions.

        Args:
            feature (str): Feature ID or name
            start (int): Start position
            end (Optional[int], optional): End position. Defaults to the same as `start`.
            strand (Optional[str], optional): Strand ('+' or '-'). Defaults to checking either.
            start_offset (Optional[int], optional): Nucleotide offset from `start`. Defaults to 0.
            end_offset (Optional[int], optional): Nucleotide offset from `end`. Defaults to 0.
            refseq (str, optional): Reference allele. Defaults to "".
            altseq (str, optional): Alternate allele. Defaults to "".
            canonical (bool, optional): Only map to the canonical transcript. Defaults to False.

        Returns:
            Union[List[ExonMappablePosition], List[ExonSmallVariant]]:
                If `refseq` or `altseq` was given, returns a variant. Otherwise returns a position.
        """
        if refseq or altseq:
            return self._map_variant(
                self.transcript_ids,
                self._rna_to_exon_variant,
                feature,
                start,
                refseq,
                altseq,
                start_offset,
                end,
                end_offset,
                strand,
                canonical,
            )
        else:
            return self._map(
                self.transcript_ids,
                self._rna_to_exon,
                feature,
                start,
                start_offset,
                end,
                end_offset,
                strand,
                canonical,
            )

    def rna_to_protein(
        self,
        feature: str,
        start: int,
        end: Optional[int] = None,
        strand: Optional[str] = None,
        start_offset: Optional[int] = None,
        end_offset: Optional[int] = None,
        refseq: str = "",
        altseq: str = "",
        canonical: bool = False,
    ) -> Union[List[ProteinMappablePosition], List[ProteinSmallVariant]]:
        """Map an RNA position to zero or more protein positions.

        Args:
            feature (str): Feature ID or name
            start (int): Start position
            end (Optional[int], optional): End position. Defaults to the same as `start`.
            strand (Optional[str], optional): Strand ('+' or '-'). Defaults to checking either.
            start_offset (Optional[int], optional): Nucleotide offset from `start`. Defaults to 0.
            end_offset (Optional[int], optional): Nucleotide offset from `end`. Defaults to 0.
            refseq (str, optional): Reference allele. Defaults to "".
            altseq (str, optional): Alternate allele. Defaults to "".
            canonical (bool, optional): Only map to the canonical transcript. Defaults to False.

        Returns:
            Union[List[ProteinMappablePosition], List[ProteinSmallVariant]]:
                If `refseq` or `altseq` was given, returns a variant. Otherwise returns a position.
        """
        if refseq or altseq:
            return self._map_variant(
                self.transcript_ids,
                self._rna_to_protein_variant,
                feature,
                start,
                refseq,
                altseq,
                start_offset,
                end,
                end_offset,
                strand,
                canonical,
            )
        else:
            return self._map(
                self.transcript_ids,
                self._rna_to_protein,
                feature,
                start,
                start_offset,
                end,
                end_offset,
                strand,
                canonical,
            )

    def rna_to_rna(
        self,
        feature: str,
        start: int,
        end: Optional[int] = None,
        strand: Optional[str] = None,
        start_offset: Optional[int] = None,
        end_offset: Optional[int] = None,
        refseq: str = "",
        altseq: str = "",
        canonical: bool = False,
    ) -> Union[List[RnaMappablePosition], List[RnaSmallVariant]]:
        """Map an RNA position to zero or more RNA positions.

        Args:
            feature (str): Feature ID or name
            start (int): Start position
            end (Optional[int], optional): End position. Defaults to the same as `start`.
            strand (Optional[str], optional): Strand ('+' or '-'). Defaults to checking either.
            start_offset (Optional[int], optional): Nucleotide offset from `start`. Defaults to 0.
            end_offset (Optional[int], optional): Nucleotide offset from `end`. Defaults to 0.
            refseq (str, optional): Reference allele. Defaults to "".
            altseq (str, optional): Alternate allele. Defaults to "".
            canonical (bool, optional): Only map to the canonical transcript. Defaults to False.

        Returns:
            Union[List[RnaMappablePosition], List[RnaSmallVariant]]:
                If `refseq` or `altseq` was given, returns a variant. Otherwise returns a position.
        """
        if refseq or altseq:
            return self._map_variant(
                self.transcript_ids,
                self._rna_to_rna_variant,
                feature,
                start,
                refseq,
                altseq,
                start_offset,
                end,
                end_offset,
                strand,
                canonical,
            )
        else:
            return self._map(
                self.transcript_ids,
                self._rna_to_rna,
                feature,
                start,
                start_offset,
                end,
                end_offset,
                strand,
                canonical,
            )

    def _cdna_to_cdna(
        self,
        transcript_ids: List[str],
        start: int,
        start_offset: int,
        end: int,
        end_offset: int,
        strand: List[str],
        canonical: bool,
        include_stop: bool = True,
    ) -> List[CdnaMappablePosition]:
        feature = ["CDS", "stop_codon"] if include_stop else ["CDS"]

        if canonical:
            transcript_ids = [i for i in transcript_ids if self.is_canonical_transcript(i)]

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
                    canonical,
                    include_stop=include_stop,
                ):
                    for cdna in dna.to_cdna(canonical=canonical):
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
            return join_positions(result_start, result_end, TRANSCRIPT_ID)

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
        canonical: bool,
        include_stop: bool = True,
    ) -> List[CdnaSmallVariant]:
        result: List[CdnaSmallVariant] = []

        for cdna in self._cdna_to_cdna(
            transcript_ids,
            start,
            start_offset,
            end,
            end_offset,
            strand,
            canonical,
            include_stop=include_stop,
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
        canonical: bool,
        include_stop: bool = True,
    ) -> List[DnaMappablePosition]:
        feature = ["CDS", "stop_codon"] if include_stop else ["CDS"]

        if canonical:
            transcript_ids = [i for i in transcript_ids if self.is_canonical_transcript(i)]

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
            return join_positions(result_start, result_end, CONTIG_ID)

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
        canonical: bool,
        include_stop: bool = True,
    ) -> List[DnaSmallVariant]:
        result: List[DnaSmallVariant] = []

        for dna in self._cdna_to_dna(
            transcript_ids,
            start,
            start_offset,
            end,
            end_offset,
            strand,
            canonical,
            include_stop=include_stop,
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
        canonical: bool,
        include_stop: bool = True,
    ) -> List[ExonMappablePosition]:
        feature = ["CDS", "stop_codon"] if include_stop else ["CDS"]

        if canonical:
            transcript_ids = [i for i in transcript_ids if self.is_canonical_transcript(i)]

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
                    canonical,
                    include_stop=include_stop,
                ):
                    for exon in dna.to_exon(canonical=canonical):
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
        if start == end:
            return sorted(result_start)
        else:
            result_end = convert(end, end_offset)
            return join_positions(result_start, result_end, TRANSCRIPT_ID)

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
        canonical: bool,
        include_stop: bool = True,
    ) -> List[ExonSmallVariant]:
        result: List[ExonSmallVariant] = []

        for exon in self._cdna_to_exon(
            transcript_ids,
            start,
            start_offset,
            end,
            end_offset,
            strand,
            canonical,
            include_stop=include_stop,
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
        canonical: bool,
    ) -> List[ProteinMappablePosition]:
        result = []
        for cdna in self._cdna_to_cdna(
            transcript_ids,
            start,
            start_offset,
            end,
            end_offset,
            strand,
            canonical,
            include_stop=False,
        ):
            # If the postion wasn't mapped to a non-offset position by `_cdna_to_cdna`, it means
            # that the position does not map to a protein.
            if cdna.start_offset or cdna.end_offset:
                continue

            # Convert the cDNA position to a protein position
            protein_start = calc_cdna_to_protein(cdna.start)
            protein_end = calc_cdna_to_protein(cdna.end)
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
        canonical: bool,
    ) -> List[ProteinSmallVariant]:
        # TODO: A lot of this function is duplicated from _cdna_to_protein()
        result: List[ProteinSmallVariant] = []
        for cdna in self._cdna_to_cdna(
            transcript_ids,
            start,
            start_offset,
            end,
            end_offset,
            strand,
            canonical,
            include_stop=False,
        ):
            # If the postion wasn't mapped to a non-offset position by `_cdna_to_cdna`, it means
            # that the position does not map to a protein.
            if cdna.start_offset or cdna.end_offset:
                continue

            # Convert the cDNA position to a protein position
            protein_start = calc_cdna_to_protein(cdna.start)
            protein_end = calc_cdna_to_protein(cdna.end)
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
        canonical: bool,
        include_stop: bool = True,
    ) -> List[RnaMappablePosition]:
        feature = ["CDS", "stop_codon"] if include_stop else ["CDS"]

        if canonical:
            transcript_ids = [i for i in transcript_ids if self.is_canonical_transcript(i)]

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
                    canonical,
                    include_stop=include_stop,
                ):
                    for rna in dna.to_rna(canonical=canonical):
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
            return join_positions(result_start, result_end, TRANSCRIPT_ID)

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
        canonical: bool,
        include_stop: bool = True,
    ) -> List[RnaSmallVariant]:
        result: List[RnaSmallVariant] = []

        for rna in self._cdna_to_rna(
            transcript_ids,
            start,
            start_offset,
            end,
            end_offset,
            strand,
            canonical,
            include_stop=include_stop,
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
        canonical: bool,
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
                    if canonical and not self.is_canonical_transcript(cds.transcript_id):
                        continue

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
            return join_positions(result_start, result_end, TRANSCRIPT_ID)

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
        canonical: bool,
        include_stop: bool = True,
    ) -> List[CdnaSmallVariant]:
        result: List[CdnaSmallVariant] = []

        for cdna in self._dna_to_cdna(
            contig_ids,
            start,
            start_offset,
            end,
            end_offset,
            strand,
            canonical,
            include_stop=include_stop,
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
        canonical: bool,
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
        canonical: bool,
    ) -> List[DnaSmallVariant]:
        result: List[DnaSmallVariant] = []

        for dna in self._dna_to_dna(
            contig_ids, start, start_offset, end, end_offset, strand, canonical
        ):
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
        canonical: bool,
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
                    if canonical and not self.is_canonical_transcript(exon.transcript_id):
                        continue

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
            return join_positions(result_start, result_end, TRANSCRIPT_ID)

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
        canonical: bool,
    ) -> List[ExonSmallVariant]:
        result: List[ExonSmallVariant] = []

        for exon in self._dna_to_exon(
            contig_ids, start, start_offset, end, end_offset, strand, canonical
        ):
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
        canonical: bool,
    ) -> List[ProteinMappablePosition]:
        def convert(x):
            return floor((x - 1) / 3 + 1)

        protein = []
        for cdna in self._dna_to_cdna(
            contig_ids, start, start_offset, end, end_offset, strand, canonical, include_stop=False
        ):
            # Offset cDNA position are assumed to not map to a protein
            if cdna.start_offset or cdna.end_offset:
                continue

            pstart = calc_cdna_to_protein(cdna.start)
            pend = calc_cdna_to_protein(cdna.end)
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
        canonical: bool,
    ) -> List[ProteinSmallVariant]:
        # TODO: A lot of this function is duplicated from _dna_to_protein()
        def convert(x):
            return floor((x - 1) / 3 + 1)

        result: List[ProteinSmallVariant] = []
        for cdna in self._dna_to_cdna(
            contig_ids, start, start_offset, end, end_offset, strand, canonical, include_stop=False
        ):
            # Offset cDNA position are assumed to not map to a protein
            if cdna.start_offset or cdna.end_offset:
                continue

            pstart = calc_cdna_to_protein(cdna.start)
            pend = calc_cdna_to_protein(cdna.end)
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
        canonical: bool,
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
                    if canonical and not self.is_canonical_transcript(exon.transcript_id):
                        continue

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
            return join_positions(result_start, result_end, TRANSCRIPT_ID)

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
        canonical: bool,
    ) -> List[RnaSmallVariant]:
        result: List[RnaSmallVariant] = []

        for rna in self._dna_to_rna(
            contig_ids, start, start_offset, end, end_offset, strand, canonical
        ):
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
        canonical: bool,
        include_stop: bool = True,
    ) -> List[CdnaMappablePosition]:
        feature = ["CDS", "stop_codon"] if include_stop else ["CDS"]

        if canonical:
            transcript_ids = [i for i in transcript_ids if self.is_canonical_transcript(i)]

        def convert(position: int, offset: int):
            result = []

            # TODO: Is there a reasonable case where an exon position would have an offset?
            assert not offset, offset

            mask = (
                (self.df[TRANSCRIPT_ID].isin(transcript_ids))
                & (self.df["exon_number"] == float(position))
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
            return join_positions(result_start, result_end, TRANSCRIPT_ID)

    def _exon_to_dna(
        self,
        transcript_ids: List[str],
        start: int,
        start_offset: int,
        end: int,
        end_offset: int,
        strand: List[str],
        canonical: bool,
    ) -> List[DnaMappablePosition]:
        if canonical:
            transcript_ids = [i for i in transcript_ids if self.is_canonical_transcript(i)]

        def convert(position: int, offset: int):
            result = []

            # TODO: Is there a reasonable case where an exon position would have an offset?
            assert not offset, offset

            mask = (
                (self.df[TRANSCRIPT_ID].isin(transcript_ids))
                & (self.df["exon_number"] == float(position))
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
            return join_positions(result_start, result_end, CONTIG_ID)

    def _exon_to_exon(
        self,
        transcript_ids: List[str],
        start: int,
        start_offset: int,
        end: int,
        end_offset: int,
        strand: List[str],
        canonical: bool,
    ) -> List[ExonMappablePosition]:
        if canonical:
            transcript_ids = [i for i in transcript_ids if self.is_canonical_transcript(i)]

        def convert(position: int, offset):
            result = []

            # TODO: Is there a reasonable case where an exon position would have an offset?
            assert not offset, offset

            mask = (
                (self.df[TRANSCRIPT_ID].isin(transcript_ids))
                & (self.df["exon_number"] == float(position))
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
            return join_positions(result_start, result_end, TRANSCRIPT_ID)

    def _exon_to_protein(
        self,
        transcript_ids: List[str],
        start: int,
        start_offset: int,
        end: int,
        end_offset: int,
        strand: List[str],
        canonical: bool,
    ) -> List[ProteinMappablePosition]:
        result = []
        for cdna in self._exon_to_cdna(
            transcript_ids,
            start,
            start_offset,
            end,
            end_offset,
            strand,
            canonical,
            include_stop=False,
        ):
            result.extend(cdna.to_protein(canonical=canonical))

        return result

    def _exon_to_rna(
        self,
        transcript_ids: List[str],
        start: int,
        start_offset: int,
        end: int,
        end_offset: int,
        strand: List[str],
        canonical: bool,
    ) -> List[RnaMappablePosition]:
        if canonical:
            transcript_ids = [i for i in transcript_ids if self.is_canonical_transcript(i)]

        def convert(position: int, offset: int):
            result = []

            # TODO: Is there a reasonable case where an exon position would have an offset?
            assert not offset, offset

            mask = (
                (self.df[TRANSCRIPT_ID].isin(transcript_ids))
                & (self.df["exon_number"] == float(position))
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
            return join_positions(result_start, result_end, TRANSCRIPT_ID)

    def _protein_to_cdna(
        self,
        transcript_ids: List[str],
        start: int,
        start_offset: int,
        end: int,
        end_offset: int,
        strand: List[str],
        canonical: bool,
    ) -> List[CdnaMappablePosition]:
        def convert(position: int):
            return ((position - 1) * 3) + 1

        # TODO: Is there a reasonable case where an protein position would have an offset?
        assert not start_offset, start_offset
        assert not end_offset, end_offset

        cdna_start = convert(start)
        cdna_end = convert(end) + 2

        return self._cdna_to_cdna(
            transcript_ids, cdna_start, 0, cdna_end, 0, strand, canonical, include_stop=False
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
        canonical: bool,
    ) -> List[CdnaSmallVariant]:
        result: List[CdnaSmallVariant] = []

        for cdna in self._protein_to_cdna(
            transcript_ids, start, start_offset, end, end_offset, strand, canonical
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
        canonical: bool,
    ) -> List[DnaMappablePosition]:
        result = []
        for cdna in self._protein_to_cdna(
            transcript_ids, start, start_offset, end, end_offset, strand, canonical
        ):
            result.extend(cdna.to_dna(canonical=canonical))

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
        canonical: bool,
    ) -> List[DnaSmallVariant]:
        result = []
        for cdna in self._protein_to_cdna_variant(
            transcript_ids, start, start_offset, end, end_offset, strand, refseq, altseq, canonical
        ):
            result.extend(cdna.to_dna(canonical=canonical))

        return sorted(set(result))

    def _protein_to_exon(
        self,
        transcript_ids: List[str],
        start: int,
        start_offset: int,
        end: int,
        end_offset: int,
        strand: List[str],
        canonical: bool,
    ) -> List[ExonMappablePosition]:
        result = []
        for cdna in self._protein_to_cdna(
            transcript_ids, start, start_offset, end, end_offset, strand, canonical
        ):
            result.extend(cdna.to_exon(canonical=canonical))

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
        canonical: bool,
    ) -> List[ExonSmallVariant]:
        result = []
        for cdna in self._protein_to_cdna_variant(
            transcript_ids, start, start_offset, end, end_offset, strand, refseq, altseq, canonical
        ):
            result.extend(cdna.to_exon(canonical=canonical))

        return sorted(set(result))

    def _protein_to_protein(
        self,
        transcript_ids: List[str],
        start: int,
        start_offset: int,
        end: int,
        end_offset: int,
        strand: List[str],
        canonical: bool,
    ) -> List[ProteinMappablePosition]:
        result = []
        for cdna in self._protein_to_cdna(
            transcript_ids, start, start_offset, end, end_offset, strand, canonical
        ):
            result.extend(cdna.to_protein(canonical=canonical))

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
        canonical: bool,
    ) -> List[ProteinSmallVariant]:
        result = []
        for cdna in self._protein_to_cdna_variant(
            transcript_ids, start, start_offset, end, end_offset, strand, refseq, altseq, canonical
        ):
            result.extend(cdna.to_protein(canonical=canonical))

        return sorted(set(result))

    def _protein_to_rna(
        self,
        transcript_ids: List[str],
        start: int,
        start_offset: int,
        end: int,
        end_offset: int,
        strand: List[str],
        canonical: bool,
    ) -> List[RnaMappablePosition]:
        result = []
        for cdna in self._protein_to_cdna(
            transcript_ids, start, start_offset, end, end_offset, strand, canonical
        ):
            result.extend(cdna.to_rna(canonical=canonical))

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
        canonical: bool,
    ) -> List[RnaSmallVariant]:
        result = []
        for cdna in self._protein_to_cdna_variant(
            transcript_ids, start, start_offset, end, end_offset, strand, refseq, altseq, canonical
        ):
            result.extend(cdna.to_rna(canonical=canonical))

        return sorted(set(result))

    def _rna_to_cdna(
        self,
        transcript_ids: List[str],
        start: int,
        start_offset: int,
        end: int,
        end_offset: int,
        strand: List[str],
        canonical: bool,
        include_stop: bool = True,
    ) -> List[CdnaMappablePosition]:
        feature = ["CDS", "stop_codon"] if include_stop else ["CDS"]

        if canonical:
            transcript_ids = [i for i in transcript_ids if self.is_canonical_transcript(i)]

        def convert(position: int, offset: int):
            result = []

            # If the offset position can be normalized to a non-offset position, do so. Otherwise
            # just return an offset position.
            if offset:
                for dna in self._rna_to_dna(
                    transcript_ids, position, offset, position, offset, strand, canonical
                ):
                    for cdna in dna.to_cdna(canonical=canonical):
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
            return join_positions(result_start, result_end, TRANSCRIPT_ID)

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
        canonical: bool,
        include_stop: bool = True,
    ) -> List[CdnaSmallVariant]:
        result: List[CdnaSmallVariant] = []

        for cdna in self._rna_to_cdna(
            transcript_ids,
            start,
            start_offset,
            end,
            end_offset,
            strand,
            canonical,
            include_stop=include_stop,
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
        canonical: bool,
    ) -> List[DnaMappablePosition]:
        if canonical:
            transcript_ids = [i for i in transcript_ids if self.is_canonical_transcript(i)]

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
            return join_positions(result_start, result_end, CONTIG_ID)

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
        canonical: bool,
    ) -> List[DnaSmallVariant]:
        result: List[DnaSmallVariant] = []

        for dna in self._rna_to_dna(
            transcript_ids, start, start_offset, end, end_offset, strand, canonical
        ):
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
        canonical: bool,
    ) -> List[ExonMappablePosition]:
        if canonical:
            transcript_ids = [i for i in transcript_ids if self.is_canonical_transcript(i)]

        def convert(position: int, offset: int):
            result = []

            # For an offset position, we need to calculate then equivalent DNA position then map
            # that to one or more exons. This is slower, so if there's no offset we can just map
            # directly to an exon.
            if offset:
                for dna in self._rna_to_dna(
                    transcript_ids, position, offset, position, offset, strand, canonical
                ):
                    for exon in dna.to_exon(canonical=canonical):
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
            return join_positions(result_start, result_end, TRANSCRIPT_ID)

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
        canonical: bool,
    ) -> List[ExonSmallVariant]:
        result: List[ExonSmallVariant] = []

        for exon in self._rna_to_exon(
            transcript_ids, start, start_offset, end, end_offset, strand, canonical
        ):
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
        canonical: bool,
    ) -> List[ProteinMappablePosition]:
        result = []
        for cdna in self._rna_to_cdna(
            transcript_ids,
            start,
            start_offset,
            end,
            end_offset,
            strand,
            canonical,
            include_stop=False,
        ):
            result.extend(cdna.to_protein(canonical=canonical))

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
        canonical: bool,
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
            canonical,
            include_stop=False,
        ):
            result.extend(cdna.to_protein(canonical=canonical))

        return sorted(set(result))

    def _rna_to_rna(
        self,
        transcript_ids: List[str],
        start: int,
        start_offset: int,
        end: int,
        end_offset: int,
        strand: List[str],
        canonical: bool,
    ) -> List[RnaMappablePosition]:
        if canonical:
            transcript_ids = [i for i in transcript_ids if self.is_canonical_transcript(i)]

        def convert(position: int, offset: int):
            result = []

            # If the offset position can be normalized to a non-offset position, do so. Otherwise
            # just return an offset position.
            if offset:
                for dna in self._rna_to_dna(
                    transcript_ids, position, offset, position, offset, strand, canonical
                ):
                    for rna in dna.to_rna(canonical=canonical):
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
            return join_positions(result_start, result_end, TRANSCRIPT_ID)

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
        canonical: bool,
    ) -> List[RnaSmallVariant]:
        result: List[RnaSmallVariant] = []

        for rna in self._rna_to_rna(
            transcript_ids, start, start_offset, end, end_offset, strand, canonical
        ):
            result.extend(RnaSmallVariant.from_rna(rna, refseq, altseq))

        return result

    def _normalize_map_args(
        self,
        start: Optional[int],
        start_offset: Optional[int],
        end: Optional[int],
        end_offset: Optional[int],
        strand: Optional[str],
    ) -> Tuple[Optional[int], int, Optional[int], int, List[str]]:
        start_ = start if start is not None else end
        end_ = end if end is not None else start
        start_offset_ = start_offset or 0
        end_offset_ = end_offset or 0
        strand_ = [strand] if strand is not None else ["+", "-"]

        return start_, start_offset_, end_, end_offset_, strand_

    def _map(
        self,
        idfunc: Callable,
        mapfunc: Callable,
        feature: str,
        start: int,
        start_offset: Optional[int],
        end: Optional[int],
        end_offset: Optional[int],
        strand: Optional[str],
        canonical: bool,
    ) -> List:
        start_, start_offset_, end_, end_offset_, strand_ = self._normalize_map_args(
            start, start_offset, end, end_offset, strand
        )
        feature_ = idfunc(feature)
        result = mapfunc(feature_, start_, start_offset_, end_, end_offset_, strand_, canonical)

        return sorted(set(result))

    def _map_variant(
        self,
        idfunc: Callable,
        mapfunc: Callable,
        feature: str,
        start: int,
        refseq: str,
        altseq: str,
        start_offset: Optional[int],
        end: Optional[int],
        end_offset: Optional[int],
        strand: Optional[str],
        canonical: bool,
    ) -> List:
        start_, start_offset_, end_, end_offset_, strand_ = self._normalize_map_args(
            start, start_offset, end, end_offset, strand
        )
        feature_ = idfunc(feature)
        result = mapfunc(
            feature_, start_, start_offset_, end_, end_offset_, strand_, refseq, altseq, canonical
        )

        return sorted(set(result))

    def _map_exon(
        self,
        mapfunc: Callable,
        feature: str,
        start: Optional[int],
        start_offset: Optional[int],
        end: Optional[int],
        end_offset: Optional[int],
        strand: Optional[str],
        canonical: bool,
    ) -> List:
        start_, start_offset_, end_, end_offset_, strand_ = self._normalize_map_args(
            start, start_offset, end, end_offset, strand
        )
        feature_ = self.transcript_ids(feature)

        if start_ is not None:
            positions = list(product(feature_, [start_], [end_], strand_))
        else:
            positions = []
            for exon_id in self.exon_ids(feature):
                positions.extend(self.exon_index(exon_id, transcript_ids=feature_))

        result = []
        for transcript_id, start__, end__, strand__ in positions:
            result.extend(
                mapfunc(
                    [transcript_id],
                    start__,
                    start_offset_,
                    end__,
                    end_offset_,
                    [strand__],
                    canonical,
                )
            )

        return sorted(set(result))

    # ---------------------------------------------------------------------------------------------
    # Utility functions
    # ---------------------------------------------------------------------------------------------
    def cds_offset(self, transcript_id: str) -> int:
        """Get the integer offset of the CDS from the start of the spliced RNA.

        Args:
            transcript_id (str): transcript ID

        Returns:
            int: Nucleotide offset of the cDNA start from the start of the transcript
        """
        rna = self.cdna_to_rna(transcript_id, 1)
        assert len(rna) == 1, f"Ambiguous transcript start for '{transcript_id}'"
        offset = rna[0].start - 1
        assert offset >= 0

        return offset

    def exon_index(
        self, exon_id: str, transcript_ids: Union[List[str], str] = []
    ) -> List[Tuple[str, int, int, str]]:
        """Get the exon number(s) for an exon ID on each transcript. Optionally, restrict to a
        specific transcript(s).

        Args:
            exon_id (str): exon ID
            transcript_ids (Union[List[str], str], optional): Restrict the search to the given
                transcript IDs. Defaults to [].

        Returns:
            List[Tuple[str, int, int, str]]: List of (exon ID, exon index, exon index, exon strand)
        """
        result = []
        restrict = [transcript_ids] if isinstance(transcript_ids, str) else transcript_ids

        mask = (self.df[EXON_ID] == exon_id) & (self.df["feature"] == "exon")
        for _, exon in self.df[mask].iterrows():
            if restrict and exon.transcript_id not in restrict:
                continue
            else:
                result.append((exon.transcript_id, exon.exon_number, exon.exon_number, exon.strand))

        return sorted(set(result))

    def translate_cds_variant(
        self, transcript_id: str, cdna_start: int, cdna_end: int, cdna_alt: str
    ) -> List[str]:
        """Return the mutated protein sequence, given a cDNA position and alt allele.

        Args:
            transcript_id (str): transcript ID
            cdna_start (int): cDNA start position
            cdna_end (int): cDNA end position
            cdna_alt (str): cDNA alternate allele

        Returns:
            List[str]: protein alternate allele(s)
        """
        pep_altseq_set = set()

        # If no alt, assume we're talking about a deletion variant
        if not cdna_alt:
            return [""]

        # Get the codon sequence
        codon_start_offset = (cdna_start - 1) % 3
        codon_start = cdna_start - codon_start_offset
        codon_end_offset = 2 - ((cdna_end - 1) % 3)
        codon_end = cdna_end + codon_end_offset
        codon_refseq = self.cds_sequence(transcript_id, codon_start, codon_end)

        # Assert that the codon sequence is divisible by 3
        assert len(codon_refseq) % 3 == 0

        # Mutate the codon sequence
        codon_refseq_left = codon_refseq[:codon_start_offset]
        codon_refseq_right = codon_refseq[-codon_end_offset:] if codon_end_offset else ""
        for i in expand_nt(cdna_alt):
            codon_altseq = codon_refseq_left + i + codon_refseq_right
            # Assert that the altered codon sequence is divisible by 3
            assert len(codon_altseq) % 3 == 0, codon_altseq
            pep_altseq = "".join(AMINO_ACID_TABLE[codon] for codon in split_by_codon(codon_altseq))
            pep_altseq_set.add(pep_altseq)

        return sorted(pep_altseq_set)


def join_positions(
    start: List[MappablePositionOrSmallVariant],
    end: List[MappablePositionOrSmallVariant],
    merge_on: str,
) -> List[MappablePositionOrSmallVariant]:
    """Return the combination of two list of position or variant objects - one of start positions
    and one of end positions - into one list by the given key (e.g. 'transcript_id').

    All positions must be of the same class.

    Args:
        start (List[MappablePositionOrSmallVariant]): Start positions
        end (List[MappablePositionOrSmallVariant]): End positions
        merge_on (str): Attribute to merge on (e.g. 'transcript_id')

    Returns:
        List[MappablePositionOrSmallVariant]: New position or variant objects
    """
    result = set()

    for start_pos, end_pos in product(start, end):
        # Make sure that the positions are sorted by position
        if start_pos.start > end_pos.start:
            start_pos, end_pos = end_pos, start_pos
        # Make a new position from the start of the 1st position and end of the 2nd
        start_key = getattr(start_pos, merge_on)
        end_key = getattr(end_pos, merge_on)
        if start_key == end_key:
            assert start_pos.__class__ == end_pos.__class__
            new = start_pos.copy_from(start_pos, end=end_pos.end, end_offset=end_pos.end_offset)
            result.add(new)

    return sorted(result)  # type: ignore
