from typing import List, Optional, Tuple

from .core import (
    CdnaMappablePosition,
    DnaMappablePosition,
    ExonMappablePosition,
    ProteinMappablePosition,
    RnaMappablePosition,
)
from .ensembl_release import EnsemblRelease


def all_contig_ids(release: Optional[int] = None, species: Optional[str] = None) -> List[str]:
    """Return a list of all contig (chromosome) IDs."""
    return EnsemblRelease.instance(release=release, species=species).all_contig_ids()


def all_exon_ids(release: Optional[int] = None, species: Optional[str] = None) -> List[str]:
    """Return a list of all exon IDs."""
    return EnsemblRelease.instance(release=release, species=species).all_exon_ids()


def all_gene_ids(release: Optional[int] = None, species: Optional[str] = None) -> List[str]:
    """Return a list of all gene IDs."""
    return EnsemblRelease.instance(release=release, species=species).all_gene_ids()


def all_gene_names(release: Optional[int] = None, species: Optional[str] = None) -> List[str]:
    """Return a list of all gene names."""
    return EnsemblRelease.instance(release=release, species=species).all_gene_names()


def all_protein_ids(release: Optional[int] = None, species: Optional[str] = None) -> List[str]:
    """Return a list of all protein IDs."""
    return EnsemblRelease.instance(release=release, species=species).all_protein_ids()


def all_transcript_ids(release: Optional[int] = None, species: Optional[str] = None) -> List[str]:
    """Return a list of all transcript IDs."""
    return EnsemblRelease.instance(release=release, species=species).all_transcript_ids()


def all_transcript_names(release: Optional[int] = None, species: Optional[str] = None) -> List[str]:
    """Return a list of all transcript names."""
    return EnsemblRelease.instance(release=release, species=species).all_transcript_names()


def contig_ids(
    feature: str,
    feature_type: str = "",
    release: Optional[int] = None,
    species: Optional[str] = None,
) -> List[str]:
    """Given a feature symbol, return the corresponding contig ID(s)."""
    return EnsemblRelease.instance(release=release, species=species).contig_ids(
        feature, feature_type=feature_type
    )


def exon_ids(
    feature: str,
    feature_type: str = "",
    release: Optional[int] = None,
    species: Optional[str] = None,
) -> List[str]:
    """Given a feature symbol, return the corresponding exon ID(s)."""
    return EnsemblRelease.instance(release=release, species=species).exon_ids(
        feature, feature_type=feature_type
    )


def gene_ids(
    feature: str,
    feature_type: str = "",
    release: Optional[int] = None,
    species: Optional[str] = None,
) -> List[str]:
    """Given a feature symbol, return the corresponding gene ID(s)."""
    return EnsemblRelease.instance(release=release, species=species).gene_ids(
        feature, feature_type=feature_type
    )


def gene_names(
    feature: str,
    feature_type: str = "",
    release: Optional[int] = None,
    species: Optional[str] = None,
) -> List[str]:
    """Given a feature symbol, return the corresponding gene names(s)."""
    return EnsemblRelease.instance(release=release, species=species).gene_names(
        feature, feature_type=feature_type
    )


def protein_ids(
    feature: str,
    feature_type: str = "",
    release: Optional[int] = None,
    species: Optional[str] = None,
) -> List[str]:
    """Given a feature symbol, return the corresponding protein ID(s)."""
    return EnsemblRelease.instance(release=release, species=species).protein_ids(
        feature, feature_type=feature_type
    )


def transcript_ids(
    feature: str,
    feature_type: str = "",
    release: Optional[int] = None,
    species: Optional[str] = None,
) -> List[str]:
    """Given a feature symbol, return the corresponding transcript ID(s)."""
    return EnsemblRelease.instance(release=release, species=species).transcript_ids(
        feature, feature_type=feature_type
    )


def transcript_names(
    feature: str,
    feature_type: str = "",
    release: Optional[int] = None,
    species: Optional[str] = None,
) -> List[str]:
    """Given a feature symbol, return the corresponding transcript names(s)."""
    return EnsemblRelease.instance(release=release, species=species).transcript_names(
        feature, feature_type=feature_type
    )


def normalize_feature(
    feature: str,
    feature_type: str = "",
    release: Optional[int] = None,
    species: Optional[str] = None,
) -> List[Tuple[str, str]]:
    """Normalize a feature ID to one or more representations."""
    return EnsemblRelease.instance(release=release, species=species).normalize_feature(
        feature, feature_type=feature_type
    )


def is_contig(feature: str, release: Optional[int] = None, species: Optional[str] = None) -> bool:
    """Return True if the given feature is a contig ID."""
    return EnsemblRelease.instance(release=release, species=species).is_contig(feature)


def is_exon(feature: str, release: Optional[int] = None, species: Optional[str] = None) -> bool:
    """Return True if the given feature is an exon ID."""
    return EnsemblRelease.instance(release=release, species=species).is_exon(feature)


def is_gene(feature: str, release: Optional[int] = None, species: Optional[str] = None) -> bool:
    """Return True if the given feature is a gene ID or name."""
    return EnsemblRelease.instance(release=release, species=species).is_gene(feature)


def is_protein(feature: str, release: Optional[int] = None, species: Optional[str] = None) -> bool:
    """Return True if the given feature is a protein ID."""
    return EnsemblRelease.instance(release=release, species=species).is_protein(feature)


def is_transcript(
    feature: str, release: Optional[int] = None, species: Optional[str] = None
) -> bool:
    """Return True if the given feature is a transcript ID or name."""
    return EnsemblRelease.instance(release=release, species=species).is_transcript(feature)


def cds_sequence(
    transcript_id: str,
    start: Optional[int] = None,
    end: Optional[int] = None,
    release: Optional[int] = None,
    species: Optional[str] = None,
) -> str:
    """Return the nucleotide sequence at the given CDS coordinates."""
    return EnsemblRelease.instance(release=release, species=species).cds_sequence(
        transcript_id, start=start, end=end
    )


def dna_sequence(
    contig_id: str,
    start: Optional[int] = None,
    end: Optional[int] = None,
    strand: str = "+",
    release: Optional[int] = None,
    species: Optional[str] = None,
) -> str:
    """Return the nucleotide sequence at the given contig coordinates."""
    return EnsemblRelease.instance(release=release, species=species).dna_sequence(
        contig_id, start=start, end=end, strand=strand
    )


def peptide_sequence(
    protein_id: str,
    start: Optional[int] = None,
    end: Optional[int] = None,
    release: Optional[int] = None,
    species: Optional[str] = None,
) -> str:
    """Return the amino acid sequence at the given peptide coordinates."""
    return EnsemblRelease.instance(release=release, species=species).peptide_sequence(
        protein_id, start=start, end=end
    )


def rna_sequence(
    transcript_id: str,
    start: Optional[int] = None,
    end: Optional[int] = None,
    release: Optional[int] = None,
    species: Optional[str] = None,
) -> str:
    """Return the nucleotide sequence at the given cDNA or ncRNA coordinates."""
    return EnsemblRelease.instance(release=release, species=species).rna_sequence(
        transcript_id, start=start, end=end
    )


def get_cdna(
    feature: str, release: Optional[int] = None, species: Optional[str] = None
) -> List[CdnaMappablePosition]:
    """Return the cDNA position(s) of the given feature."""
    return EnsemblRelease.instance(release=release, species=species).get_cdna(feature)


def get_dna(
    feature: str, release: Optional[int] = None, species: Optional[str] = None
) -> List[DnaMappablePosition]:
    """Return the DNA position(s) of the given feature."""
    return EnsemblRelease.instance(release=release, species=species).get_dna(feature)


def get_exons(
    feature: str, release: Optional[int] = None, species: Optional[str] = None
) -> List[ExonMappablePosition]:
    """Return the exon position(s) of the given feature."""
    return EnsemblRelease.instance(release=release, species=species).get_exons(feature)


def get_genes(
    feature: str, release: Optional[int] = None, species: Optional[str] = None
) -> List[DnaMappablePosition]:
    """Return the gene position(s) of the given feature."""
    return EnsemblRelease.instance(release=release, species=species).get_genes(feature)


def get_transcripts(
    feature: str, release: Optional[int] = None, species: Optional[str] = None
) -> List[RnaMappablePosition]:
    """Return the transcript position(s) of the given feature."""
    return EnsemblRelease.instance(release=release, species=species).get_transcripts(feature)


def cdna_to_cdna(
    feature: str,
    start: int,
    end: Optional[int] = None,
    strand: Optional[str] = None,
    start_offset: Optional[int] = None,
    end_offset: Optional[int] = None,
    release: Optional[int] = None,
    species: Optional[str] = None,
) -> List[CdnaMappablePosition]:
    """Map a cDNA position to zero or more cDNA positions."""
    return EnsemblRelease.instance(release=release, species=species).cdna_to_cdna(
        feature, start, end=end, strand=strand, start_offset=start_offset, end_offset=end_offset
    )


def cdna_to_dna(
    feature: str,
    start: int,
    end: Optional[int] = None,
    strand: Optional[str] = None,
    start_offset: Optional[int] = None,
    end_offset: Optional[int] = None,
    release: Optional[int] = None,
    species: Optional[str] = None,
) -> List[DnaMappablePosition]:
    """Map a cDNA position to zero or more DNA positions."""
    return EnsemblRelease.instance(release=release, species=species).cdna_to_dna(
        feature, start, end=end, strand=strand, start_offset=start_offset, end_offset=end_offset
    )


def cdna_to_exon(
    feature: str,
    start: int,
    end: Optional[int] = None,
    strand: Optional[str] = None,
    start_offset: Optional[int] = None,
    end_offset: Optional[int] = None,
    release: Optional[int] = None,
    species: Optional[str] = None,
) -> List[ExonMappablePosition]:
    """Map a cDNA position to an exon positions."""
    return EnsemblRelease.instance(release=release, species=species).cdna_to_exon(
        feature, start, end=end, strand=strand, start_offset=start_offset, end_offset=end_offset
    )


def cdna_to_protein(
    feature: str,
    start: int,
    end: Optional[int] = None,
    strand: Optional[str] = None,
    start_offset: Optional[int] = None,
    end_offset: Optional[int] = None,
    release: Optional[int] = None,
    species: Optional[str] = None,
) -> List[ProteinMappablePosition]:
    """Map a cDNA position to zero or more protein positions."""
    return EnsemblRelease.instance(release=release, species=species).cdna_to_protein(
        feature, start, end=end, strand=strand, start_offset=start_offset, end_offset=end_offset
    )


def cdna_to_rna(
    feature: str,
    start: int,
    end: Optional[int] = None,
    strand: Optional[str] = None,
    start_offset: Optional[int] = None,
    end_offset: Optional[int] = None,
    release: Optional[int] = None,
    species: Optional[str] = None,
) -> List[RnaMappablePosition]:
    """Map a cDNA position to zero or more RNA positions."""
    return EnsemblRelease.instance(release=release, species=species).cdna_to_rna(
        feature, start, end=end, strand=strand, start_offset=start_offset, end_offset=end_offset
    )


def dna_to_cdna(
    feature: str,
    start: int,
    end: Optional[int] = None,
    strand: Optional[str] = None,
    start_offset: Optional[int] = None,
    end_offset: Optional[int] = None,
    release: Optional[int] = None,
    species: Optional[str] = None,
) -> List[CdnaMappablePosition]:
    """Map a DNA position to zero or more cDNA positions."""
    return EnsemblRelease.instance(release=release, species=species).dna_to_cdna(
        feature, start, end=end, strand=strand, start_offset=start_offset, end_offset=end_offset
    )


def dna_to_dna(
    feature: str,
    start: int,
    end: Optional[int] = None,
    strand: Optional[str] = None,
    start_offset: Optional[int] = None,
    end_offset: Optional[int] = None,
    release: Optional[int] = None,
    species: Optional[str] = None,
) -> List[DnaMappablePosition]:
    """Map a DNA position to zero or more DNA positions."""
    return EnsemblRelease.instance(release=release, species=species).dna_to_dna(
        feature, start, end=end, strand=strand, start_offset=start_offset, end_offset=end_offset
    )


def dna_to_exon(
    feature: str,
    start: int,
    end: Optional[int] = None,
    strand: Optional[str] = None,
    start_offset: Optional[int] = None,
    end_offset: Optional[int] = None,
    release: Optional[int] = None,
    species: Optional[str] = None,
) -> List[ExonMappablePosition]:
    """Map a DNA position to zero or more exon positions."""
    return EnsemblRelease.instance(release=release, species=species).dna_to_exon(
        feature, start, end=end, strand=strand, start_offset=start_offset, end_offset=end_offset
    )


def dna_to_protein(
    feature: str,
    start: int,
    end: Optional[int] = None,
    strand: Optional[str] = None,
    start_offset: Optional[int] = None,
    end_offset: Optional[int] = None,
    release: Optional[int] = None,
    species: Optional[str] = None,
) -> List[ProteinMappablePosition]:
    """Map a DNA position to zero or more protein positions."""
    return EnsemblRelease.instance(release=release, species=species).dna_to_protein(
        feature, start, end=end, strand=strand, start_offset=start_offset, end_offset=end_offset
    )


def dna_to_rna(
    feature: str,
    start: int,
    end: Optional[int] = None,
    strand: Optional[str] = None,
    start_offset: Optional[int] = None,
    end_offset: Optional[int] = None,
    release: Optional[int] = None,
    species: Optional[str] = None,
) -> List[RnaMappablePosition]:
    """Map a DNA position to zero or more RNA positions."""
    return EnsemblRelease.instance(release=release, species=species).dna_to_rna(
        feature, start, end=end, strand=strand, start_offset=start_offset, end_offset=end_offset
    )


def exon_to_cdna(
    feature: str,
    start: Optional[int] = None,
    end: Optional[int] = None,
    strand: Optional[str] = None,
    start_offset: Optional[int] = None,
    end_offset: Optional[int] = None,
    release: Optional[int] = None,
    species: Optional[str] = None,
) -> List[CdnaMappablePosition]:
    """Map an exon position to zero or more cDNA positions."""
    return EnsemblRelease.instance(release=release, species=species).exon_to_cdna(
        feature,
        start=start,
        end=end,
        strand=strand,
        start_offset=start_offset,
        end_offset=end_offset,
    )


def exon_to_dna(
    feature: str,
    start: Optional[int] = None,
    end: Optional[int] = None,
    strand: Optional[str] = None,
    start_offset: Optional[int] = None,
    end_offset: Optional[int] = None,
    release: Optional[int] = None,
    species: Optional[str] = None,
) -> List[DnaMappablePosition]:
    """Map an exon position to zero or more DNA positions."""
    return EnsemblRelease.instance(release=release, species=species).exon_to_dna(
        feature,
        start=start,
        end=end,
        strand=strand,
        start_offset=start_offset,
        end_offset=end_offset,
    )


def exon_to_exon(
    feature: str,
    start: Optional[int] = None,
    end: Optional[int] = None,
    strand: Optional[str] = None,
    start_offset: Optional[int] = None,
    end_offset: Optional[int] = None,
    release: Optional[int] = None,
    species: Optional[str] = None,
) -> List[ExonMappablePosition]:
    """Map an exon position to an exon positions."""
    return EnsemblRelease.instance(release=release, species=species).exon_to_exon(
        feature,
        start=start,
        end=end,
        strand=strand,
        start_offset=start_offset,
        end_offset=end_offset,
    )


def exon_to_protein(
    feature: str,
    start: Optional[int] = None,
    end: Optional[int] = None,
    strand: Optional[str] = None,
    start_offset: Optional[int] = None,
    end_offset: Optional[int] = None,
    release: Optional[int] = None,
    species: Optional[str] = None,
) -> List[ProteinMappablePosition]:
    """Map an exon position to zero or more protein positions."""
    return EnsemblRelease.instance(release=release, species=species).exon_to_protein(
        feature,
        start=start,
        end=end,
        strand=strand,
        start_offset=start_offset,
        end_offset=end_offset,
    )


def exon_to_rna(
    feature: str,
    start: Optional[int] = None,
    end: Optional[int] = None,
    strand: Optional[str] = None,
    start_offset: Optional[int] = None,
    end_offset: Optional[int] = None,
    release: Optional[int] = None,
    species: Optional[str] = None,
) -> List[RnaMappablePosition]:
    """Map an exon position to zero or more RNA positions."""
    return EnsemblRelease.instance(release=release, species=species).exon_to_rna(
        feature,
        start=start,
        end=end,
        strand=strand,
        start_offset=start_offset,
        end_offset=end_offset,
    )


def protein_to_cdna(
    feature: str,
    start: int,
    end: Optional[int] = None,
    strand: Optional[str] = None,
    start_offset: Optional[int] = None,
    end_offset: Optional[int] = None,
    release: Optional[int] = None,
    species: Optional[str] = None,
) -> List[CdnaMappablePosition]:
    """Map a protein position to zero or more cDNA positions."""
    return EnsemblRelease.instance(release=release, species=species).protein_to_cdna(
        feature, start, end=end, strand=strand, start_offset=start_offset, end_offset=end_offset
    )


def protein_to_dna(
    feature: str,
    start: int,
    end: Optional[int] = None,
    strand: Optional[str] = None,
    start_offset: Optional[int] = None,
    end_offset: Optional[int] = None,
    release: Optional[int] = None,
    species: Optional[str] = None,
) -> List[DnaMappablePosition]:
    """Map a protein position to zero or more DNA positions."""
    return EnsemblRelease.instance(release=release, species=species).protein_to_dna(
        feature, start, end=end, strand=strand, start_offset=start_offset, end_offset=end_offset
    )


def protein_to_exon(
    feature: str,
    start: int,
    end: Optional[int] = None,
    strand: Optional[str] = None,
    start_offset: Optional[int] = None,
    end_offset: Optional[int] = None,
    release: Optional[int] = None,
    species: Optional[str] = None,
) -> List[ExonMappablePosition]:
    """Map a protein position to zero or more exon positions."""
    return EnsemblRelease.instance(release=release, species=species).protein_to_exon(
        feature, start, end=end, strand=strand, start_offset=start_offset, end_offset=end_offset
    )


def protein_to_protein(
    feature: str,
    start: int,
    end: Optional[int] = None,
    strand: Optional[str] = None,
    start_offset: Optional[int] = None,
    end_offset: Optional[int] = None,
    release: Optional[int] = None,
    species: Optional[str] = None,
) -> List[ProteinMappablePosition]:
    """Map a protein position to zero or more protein positions."""
    return EnsemblRelease.instance(release=release, species=species).protein_to_protein(
        feature, start, end=end, strand=strand, start_offset=start_offset, end_offset=end_offset
    )


def protein_to_rna(
    feature: str,
    start: int,
    end: Optional[int] = None,
    strand: Optional[str] = None,
    start_offset: Optional[int] = None,
    end_offset: Optional[int] = None,
    release: Optional[int] = None,
    species: Optional[str] = None,
) -> List[RnaMappablePosition]:
    """Map a protein position to zero or more RNA positions."""
    return EnsemblRelease.instance(release=release, species=species).protein_to_rna(
        feature, start, end=end, strand=strand, start_offset=start_offset, end_offset=end_offset
    )


def rna_to_cdna(
    feature: str,
    start: int,
    end: Optional[int] = None,
    strand: Optional[str] = None,
    start_offset: Optional[int] = None,
    end_offset: Optional[int] = None,
    release: Optional[int] = None,
    species: Optional[str] = None,
) -> List[CdnaMappablePosition]:
    """Map a RNA position to zero or more cDNA positions."""
    return EnsemblRelease.instance(release=release, species=species).rna_to_cdna(
        feature, start, end=end, strand=strand, start_offset=start_offset, end_offset=end_offset
    )


def rna_to_dna(
    feature: str,
    start: int,
    end: Optional[int] = None,
    strand: Optional[str] = None,
    start_offset: Optional[int] = None,
    end_offset: Optional[int] = None,
    release: Optional[int] = None,
    species: Optional[str] = None,
) -> List[DnaMappablePosition]:
    """Map a RNA position to zero or more DNA positions."""
    return EnsemblRelease.instance(release=release, species=species).rna_to_dna(
        feature, start, end=end, strand=strand, start_offset=start_offset, end_offset=end_offset
    )


def rna_to_exon(
    feature: str,
    start: int,
    end: Optional[int] = None,
    strand: Optional[str] = None,
    start_offset: Optional[int] = None,
    end_offset: Optional[int] = None,
    release: Optional[int] = None,
    species: Optional[str] = None,
) -> List[ExonMappablePosition]:
    """Map a RNA position to an exon positions."""
    return EnsemblRelease.instance(release=release, species=species).rna_to_exon(
        feature, start, end=end, strand=strand, start_offset=start_offset, end_offset=end_offset
    )


def rna_to_protein(
    feature: str,
    start: int,
    end: Optional[int] = None,
    strand: Optional[str] = None,
    start_offset: Optional[int] = None,
    end_offset: Optional[int] = None,
    release: Optional[int] = None,
    species: Optional[str] = None,
) -> List[ProteinMappablePosition]:
    """Map a RNA position to zero or more protein positions."""
    return EnsemblRelease.instance(release=release, species=species).rna_to_protein(
        feature, start, end=end, strand=strand, start_offset=start_offset, end_offset=end_offset
    )


def rna_to_rna(
    feature: str,
    start: int,
    end: Optional[int] = None,
    strand: Optional[str] = None,
    start_offset: Optional[int] = None,
    end_offset: Optional[int] = None,
    release: Optional[int] = None,
    species: Optional[str] = None,
) -> List[RnaMappablePosition]:
    """Map a RNA position to zero or more RNA positions."""
    return EnsemblRelease.instance(release=release, species=species).rna_to_rna(
        feature, start, end=end, strand=strand, start_offset=start_offset, end_offset=end_offset
    )
