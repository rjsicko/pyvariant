import pytest

from pyvariant.constants import DELINS
from pyvariant.positions import CdnaDelins, DnaDelins, ProteinDelins, RnaDelins


@pytest.fixture()
def variant():
    return ProteinDelins(
        contig_id="4",
        start=417,
        start_offset=0,
        end=419,
        end_offset=0,
        strand="+",
        gene_id="ENSG00000157404",
        gene_name="KIT",
        transcript_id="ENST00000288135",
        transcript_name="KIT-001",
        protein_id="ENSP00000288135",
        refseq="TYD",
        altseq="I",
    )


def test_str(variant):
    assert str(variant) == "ENSP00000288135:p.T417_D419delinsI"


def test_type(variant):
    assert not variant.is_cdna
    assert not variant.is_dna
    assert not variant.is_exon
    assert variant.is_protein
    assert not variant.is_rna


def test_variant_type(variant):
    assert variant.variant_type == DELINS
    assert not variant.is_deletion
    assert variant.is_delins
    assert not variant.is_duplication
    assert not variant.is_frameshift
    assert not variant.is_fusion
    assert not variant.is_insertion
    assert not variant.is_substitution


def test_to_cdna(ensembl69, variant):
    expected = CdnaDelins(
        contig_id="4",
        start=1250,
        start_offset=0,
        end=1257,
        end_offset=0,
        strand="+",
        gene_id="ENSG00000157404",
        gene_name="KIT",
        transcript_id="ENST00000288135",
        transcript_name="KIT-001",
        protein_id="ENSP00000288135",
        refseq="CTTACGAC",
        altseq="TA",
    )
    result = ensembl69.to_cdna(variant)
    assert expected in result
    assert len(result) == 3


def test_to_dna(ensembl69, variant):
    expected = DnaDelins(
        contig_id="4",
        start=55589768,
        start_offset=0,
        end=55589775,
        end_offset=0,
        strand="+",
        refseq="CTTACGAC",
        altseq="TA",
    )
    result = ensembl69.to_dna(variant)
    assert expected in result
    assert len(result) == 3


def test_to_protein(ensembl69, variant):
    expected = [
        ProteinDelins(
            contig_id="4",
            start=417,
            start_offset=0,
            end=419,
            end_offset=0,
            strand="+",
            gene_id="ENSG00000157404",
            gene_name="KIT",
            transcript_id="ENST00000288135",
            transcript_name="KIT-001",
            protein_id="ENSP00000288135",
            refseq="TYD",
            altseq="I",
        )
    ]
    assert ensembl69.to_protein(variant) == expected


def test_to_rna(ensembl69, variant):
    expected = RnaDelins(
        contig_id="4",
        start=1347,
        start_offset=0,
        end=1354,
        end_offset=0,
        strand="+",
        gene_id="ENSG00000157404",
        gene_name="KIT",
        transcript_id="ENST00000288135",
        transcript_name="KIT-001",
        refseq="CTTACGAC",
        altseq="TA",
    )
    result = ensembl69.to_rna(variant)
    assert expected in result
    assert len(result) == 3
