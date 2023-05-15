import pytest

from variant_map.constants import SUBSTITUTION
from variant_map.positions import (
    CdnaSubstitution,
    DnaSubstitution,
    ProteinSubstitution,
    RnaSubstitution,
)


@pytest.fixture()
def variant():
    return RnaSubstitution(
        contig_id="12",
        start=96,
        start_offset=0,
        end=96,
        end_offset=0,
        strand="-",
        gene_id="ENSG00000133703",
        gene_name="KRAS",
        transcript_id="ENST00000256078",
        transcript_name="KRAS-004",
        refseq="C",
        altseq="A",
    )


def test_str(variant):
    assert str(variant) == "ENST00000256078:r.96C>A"


def test_type(variant):
    assert not variant.is_cdna
    assert not variant.is_dna
    assert not variant.is_exon
    assert not variant.is_protein
    assert variant.is_rna


def test_variant_type(variant):
    assert variant.variant_type == SUBSTITUTION
    assert not variant.is_deletion
    assert not variant.is_delins
    assert not variant.is_duplication
    assert not variant.is_frameshift
    assert not variant.is_fusion
    assert not variant.is_insertion
    assert variant.is_substitution


def test_to_cdna(ensembl69, variant):
    expected = [
        CdnaSubstitution(
            contig_id="12",
            start=32,
            start_offset=0,
            end=32,
            end_offset=0,
            strand="-",
            gene_id="ENSG00000133703",
            gene_name="KRAS",
            transcript_id="ENST00000256078",
            transcript_name="KRAS-004",
            protein_id="ENSP00000256078",
            refseq="C",
            altseq="A",
        )
    ]
    assert ensembl69.to_cdna(variant) == expected


def test_to_dna(ensembl69, variant):
    expected = [
        DnaSubstitution(
            contig_id="12",
            start=25398287,
            start_offset=0,
            end=25398287,
            end_offset=0,
            strand="-",
            refseq="C",
            altseq="A",
        )
    ]
    assert ensembl69.to_dna(variant) == expected


def test_to_protein(ensembl69, variant):
    expected = [
        ProteinSubstitution(
            contig_id="12",
            start=11,
            start_offset=0,
            end=11,
            end_offset=0,
            strand="-",
            gene_id="ENSG00000133703",
            gene_name="KRAS",
            transcript_id="ENST00000256078",
            transcript_name="KRAS-004",
            protein_id="ENSP00000256078",
            refseq="A",
            altseq="D",
        )
    ]
    assert ensembl69.to_protein(variant) == expected


def test_to_rna(ensembl69, variant):
    expected = [
        RnaSubstitution(
            contig_id="12",
            start=96,
            start_offset=0,
            end=96,
            end_offset=0,
            strand="-",
            gene_id="ENSG00000133703",
            gene_name="KRAS",
            transcript_id="ENST00000256078",
            transcript_name="KRAS-004",
            refseq="C",
            altseq="A",
        )
    ]
    assert ensembl69.to_rna(variant) == expected
