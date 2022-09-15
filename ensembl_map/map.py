from typing import List, Optional

from .core import CdnaPosition, DnaPosition, ExonPosition, ProteinPosition, RnaPosition, instance


def cdna_to_cdna(
    feature: str,
    start: int,
    end: Optional[int] = None,
    strand: Optional[str] = None,
) -> List[CdnaPosition]:
    """Map a cDNA position to a cDNA position."""
    return instance().cdna_to_cdna(feature, start, end, strand)


def cdna_to_dna(
    feature: str,
    start: int,
    end: Optional[int] = None,
    strand: Optional[str] = None,
) -> List[DnaPosition]:
    """Map a cDNA position to a DNA position."""
    return instance().cdna_to_dna(feature, start, end, strand)


def cdna_to_exon(
    feature: str,
    start: int,
    end: Optional[int] = None,
    strand: Optional[str] = None,
) -> List[ExonPosition]:
    """Map a cDNA position to an exon position."""
    return instance().cdna_to_exon(feature, start, end, strand)


def cdna_to_protein(
    feature: str,
    start: int,
    end: Optional[int] = None,
    strand: Optional[str] = None,
) -> List[ProteinPosition]:
    """Map a cDNA position to a protein position."""
    return instance().cdna_to_protein(feature, start, end, strand)


def cdna_to_rna(
    feature: str,
    start: int,
    end: Optional[int] = None,
    strand: Optional[str] = None,
) -> List[RnaPosition]:
    """Map a cDNA position to a RNA position."""
    return instance().cdna_to_rna(feature, start, end, strand)


def dna_to_cdna(
    feature: str,
    start: int,
    end: Optional[int] = None,
    strand: Optional[str] = None,
) -> List[CdnaPosition]:
    """Map a DNA position to a cDNA position."""
    return instance().dna_to_cdna(feature, start, end, strand)


def dna_to_dna(
    feature: str,
    start: int,
    end: Optional[int] = None,
    strand: Optional[str] = None,
) -> List[DnaPosition]:
    """Map a DNA position to a DNA position."""
    return instance().dna_to_dna(feature, start, end, strand)


def dna_to_exon(
    feature: str,
    start: int,
    end: Optional[int] = None,
    strand: Optional[str] = None,
) -> List[ExonPosition]:
    """Map a DNA position to an exon position."""
    return instance().dna_to_exon(feature, start, end, strand)


def dna_to_protein(
    feature: str,
    start: int,
    end: Optional[int] = None,
    strand: Optional[str] = None,
) -> List[ProteinPosition]:
    """Map a DNA position to a protein position."""
    return instance().dna_to_protein(feature, start, end, strand)


def dna_to_rna(
    feature: str,
    start: int,
    end: Optional[int] = None,
    strand: Optional[str] = None,
) -> List[RnaPosition]:
    """Map a DNA position to a RNA position."""
    return instance().dna_to_rna(feature, start, end, strand)


def exon_to_cdna(
    feature: str,
    start: Optional[int] = None,
    end: Optional[int] = None,
    strand: Optional[str] = None,
) -> List[CdnaPosition]:
    """Map an exon position to a cDNA position."""
    return instance().exon_to_cdna(feature, start, end, strand)


def exon_to_dna(
    feature: str,
    start: Optional[int] = None,
    end: Optional[int] = None,
    strand: Optional[str] = None,
) -> List[DnaPosition]:
    """Map an exon position to a DNA position."""
    return instance().exon_to_dna(feature, start, end, strand)


def exon_to_exon(
    feature: str,
    start: Optional[int] = None,
    end: Optional[int] = None,
    strand: Optional[str] = None,
) -> List[ExonPosition]:
    """Map an exon position to an exon position."""
    return instance().exon_to_exon(feature, start, end, strand)


def exon_to_protein(
    feature: str,
    start: Optional[int] = None,
    end: Optional[int] = None,
    strand: Optional[str] = None,
) -> List[ProteinPosition]:
    """Map an exon position to a protein position."""
    return instance().exon_to_protein(feature, start, end, strand)


def exon_to_rna(
    feature: str,
    start: Optional[int] = None,
    end: Optional[int] = None,
    strand: Optional[str] = None,
) -> List[RnaPosition]:
    """Map an exon position to a RNA position."""
    return instance().exon_to_rna(feature, start, end, strand)


def protein_to_cdna(
    feature: str,
    start: int,
    end: Optional[int] = None,
    strand: Optional[str] = None,
) -> List[CdnaPosition]:
    """Map a protein position to a cDNA position."""
    return instance().protein_to_cdna(feature, start, end, strand)


def protein_to_dna(
    feature: str,
    start: int,
    end: Optional[int] = None,
    strand: Optional[str] = None,
) -> List[DnaPosition]:
    """Map a protein position to a DNA position."""
    return instance().protein_to_dna(feature, start, end, strand)


def protein_to_exon(
    feature: str,
    start: int,
    end: Optional[int] = None,
    strand: Optional[str] = None,
) -> List[ExonPosition]:
    """Map a protein position to an exon position."""
    return instance().protein_to_exon(feature, start, end, strand)


def protein_to_protein(
    feature: str,
    start: int,
    end: Optional[int] = None,
    strand: Optional[str] = None,
) -> List[ProteinPosition]:
    """Map a protein position to a protein position."""
    return instance().protein_to_protein(feature, start, end, strand)


def protein_to_rna(
    feature: str,
    start: int,
    end: Optional[int] = None,
    strand: Optional[str] = None,
) -> List[RnaPosition]:
    """Map a protein position to a RNA position."""
    return instance().protein_to_rna(feature, start, end, strand)


def rna_to_cdna(
    feature: str,
    start: int,
    end: Optional[int] = None,
    strand: Optional[str] = None,
) -> List[CdnaPosition]:
    """Map a RNA position to a cDNA position."""
    return instance().rna_to_cdna(feature, start, end, strand)


def rna_to_dna(
    feature: str,
    start: int,
    end: Optional[int] = None,
    strand: Optional[str] = None,
) -> List[DnaPosition]:
    """Map a RNA position to a DNA position."""
    return instance().rna_to_dna(feature, start, end, strand)


def rna_to_exon(
    feature: str,
    start: int,
    end: Optional[int] = None,
    strand: Optional[str] = None,
) -> List[ExonPosition]:
    """Map a RNA position to an exon position."""
    return instance().rna_to_exon(feature, start, end, strand)


def rna_to_protein(
    feature: str,
    start: int,
    end: Optional[int] = None,
    strand: Optional[str] = None,
) -> List[ProteinPosition]:
    """Map a RNA position to a protein position."""
    return instance().rna_to_protein(feature, start, end, strand)


def rna_to_rna(
    feature: str,
    start: int,
    end: Optional[int] = None,
    strand: Optional[str] = None,
) -> List[RnaPosition]:
    """Map a RNA position to a RNA position."""
    return instance().rna_to_rna(feature, start, end, strand)
