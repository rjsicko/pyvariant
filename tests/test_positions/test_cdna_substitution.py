import pytest

from pyvariant.constants import SUBSTITUTION
from pyvariant.positions import (
    CdnaSubstitution,
    DnaSubstitution,
    ProteinSubstitution,
    RnaSubstitution,
)


@pytest.fixture()
def variant(ensembl69):
    return CdnaSubstitution(
        _core=ensembl69,
        contig_id="12",
        start=38,
        start_offset=0,
        end=38,
        end_offset=0,
        strand="-",
        gene_id="ENSG00000133703",
        gene_name="KRAS",
        transcript_id="ENST00000256078",
        transcript_name="KRAS-004",
        protein_id="ENSP00000256078",
        refseq="G",
        altseq="A",
    )


def test_str(variant):
    assert str(variant) == "ENST00000256078:c.38G>A"


def test_type(variant):
    assert variant.is_cdna
    assert not variant.is_dna
    assert not variant.is_exon
    assert not variant.is_protein
    assert not variant.is_rna


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
            _core=ensembl69,
            contig_id="12",
            start=38,
            start_offset=0,
            end=38,
            end_offset=0,
            strand="-",
            gene_id="ENSG00000133703",
            gene_name="KRAS",
            transcript_id="ENST00000256078",
            transcript_name="KRAS-004",
            protein_id="ENSP00000256078",
            refseq="G",
            altseq="A",
        )
    ]
    assert variant.to_cdna() == expected


def test_to_dna(ensembl69, variant):
    expected = [
        DnaSubstitution(
            _core=ensembl69,
            contig_id="12",
            start=25398281,
            start_offset=0,
            end=25398281,
            end_offset=0,
            strand="-",
            refseq="G",
            altseq="A",
        )
    ]
    assert variant.to_dna() == expected


def test_to_protein(ensembl69, variant):
    expected = [
        ProteinSubstitution(
            _core=ensembl69,
            contig_id="12",
            start=13,
            start_offset=0,
            end=13,
            end_offset=0,
            strand="-",
            gene_id="ENSG00000133703",
            gene_name="KRAS",
            transcript_id="ENST00000256078",
            transcript_name="KRAS-004",
            protein_id="ENSP00000256078",
            refseq="G",  # GGC
            altseq="D",  # GAC
        )
    ]
    assert variant.to_protein() == expected


def test_to_rna(ensembl69, variant):
    expected = [
        RnaSubstitution(
            _core=ensembl69,
            contig_id="12",
            start=102,
            start_offset=0,
            end=102,
            end_offset=0,
            strand="-",
            gene_id="ENSG00000133703",
            gene_name="KRAS",
            transcript_id="ENST00000256078",
            transcript_name="KRAS-004",
            refseq="G",
            altseq="A",
        )
    ]
    assert variant.to_rna() == expected
