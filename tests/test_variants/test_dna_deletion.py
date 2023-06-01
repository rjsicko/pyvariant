import pytest

from pyvariant.constants import DELETION
from pyvariant.variants import CdnaDeletion, DnaDeletion, ProteinDeletion, RnaDeletion


@pytest.fixture()
def variant(ensembl69):
    return DnaDeletion(
        _core=ensembl69,
        contig_id="17",
        start=7579471,
        start_offset=0,
        end=7579473,
        end_offset=0,
        strand="-",
        refseq="CCC",
        altseq="",
    )


def test_str(variant):
    assert str(variant) == "17:g.7579471_7579473del"


def test_type(variant):
    assert not variant.is_cdna
    assert variant.is_dna
    assert not variant.is_exon
    assert not variant.is_protein
    assert not variant.is_rna


def test_variant_type(variant):
    assert variant.variant_type == DELETION
    assert variant.is_deletion
    assert not variant.is_delins
    assert not variant.is_duplication
    assert not variant.is_frameshift
    assert not variant.is_fusion
    assert not variant.is_insertion
    assert not variant.is_substitution


def test_to_cdna(ensembl69, variant):
    expected = CdnaDeletion(
        _core=ensembl69,
        contig_id="17",
        start=214,
        start_offset=0,
        end=216,
        end_offset=0,
        gene_id="ENSG00000141510",
        gene_name="TP53",
        transcript_id="ENST00000269305",
        transcript_name="TP53-001",
        protein_id="ENSP00000269305",
        strand="-",
        refseq="CCC",
        altseq="",
    )
    result = variant.to_cdna()
    assert expected in result
    assert len(result) == 8


def test_to_dna(ensembl69, variant):
    expected = DnaDeletion(
        _core=ensembl69,
        contig_id="17",
        start=7579471,
        start_offset=0,
        end=7579473,
        end_offset=0,
        strand="-",
        refseq="CCC",
        altseq="",
    )
    result = variant.to_dna()
    assert expected in result
    assert len(result) == 1


def test_to_protein(ensembl69, variant):
    expected = ProteinDeletion(
        _core=ensembl69,
        contig_id="17",
        start=72,
        start_offset=0,
        end=72,
        end_offset=0,
        strand="-",
        gene_id="ENSG00000141510",
        gene_name="TP53",
        transcript_id="ENST00000269305",
        transcript_name="TP53-001",
        protein_id="ENSP00000269305",
        refseq="P",
        altseq="",
    )
    result = variant.to_protein()
    assert expected in result
    assert len(result) == 8


def test_to_rna(ensembl69, variant):
    expected = RnaDeletion(
        _core=ensembl69,
        contig_id="17",
        start=404,
        start_offset=0,
        end=406,
        end_offset=0,
        strand="-",
        gene_id="ENSG00000141510",
        gene_name="TP53",
        transcript_id="ENST00000269305",
        transcript_name="TP53-001",
        refseq="CCC",
        altseq="",
    )
    result = variant.to_rna()
    assert expected in result
    assert len(result) == 9
