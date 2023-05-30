import pytest

from pyvariant.constants import SUBSTITUTION
from pyvariant.variants import ExonSmallVariant


@pytest.fixture()
def variant(ensembl69):
    return ExonSmallVariant(
        _core=ensembl69,
        contig_id="12",
        start=1,
        start_offset=0,
        end=1,
        end_offset=0,
        strand="-",
        gene_id="ENSG00000133703",
        gene_name="KRAS",
        transcript_id="ENST00000256078",
        transcript_name="KRAS-004",
        exon_id="ENSE00002446502",
        refseq="C",
        altseq="A",
    )


def test_str(variant):
    assert str(variant) == "ENST00000256078:e.1mut"


def test_to_string_gene_name(variant):
    assert variant.to_string(reference="gene_name") == "KRAS:e.1mut"


def test_type(variant):
    assert not variant.is_cdna
    assert not variant.is_dna
    assert variant.is_exon
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


# TODO: Is there a reasonable way to map a exon small variant to DNA/cDNA/etc coordinates?
def test_to_cdna(ensembl69, variant):
    expected = []
    assert variant.to_cdna() == expected


def test_to_dna(ensembl69, variant):
    expected = []
    assert variant.to_dna() == expected


def test_to_protein(ensembl69, variant):
    expected = []
    assert variant.to_protein() == expected


def test_to_rna(ensembl69, variant):
    expected = []
    assert variant.to_rna() == expected
