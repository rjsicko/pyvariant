import pytest

from pyvariant.constants import DELINS
from pyvariant.variants import CdnaDelins, DnaDelins, ProteinDelins, RnaDelins


@pytest.fixture()
def variant(ensembl69):
    return DnaDelins(
        _core=ensembl69,
        contig_id="17",
        start=7577058,
        start_offset=0,
        end=7577060,
        end_offset=0,
        strand="-",
        refseq="GGG",
        altseq="TTT",
    )


def test_str(variant):
    assert str(variant) == "17:g.7577058_7577060delinsTTT"


def test_type(variant):
    assert not variant.is_cdna
    assert variant.is_dna
    assert not variant.is_exon
    assert not variant.is_protein
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
        _core=ensembl69,
        contig_id="17",
        start=878,
        start_offset=0,
        end=880,
        end_offset=0,
        strand="-",
        gene_id="ENSG00000141510",
        gene_name="TP53",
        transcript_id="ENST00000269305",
        transcript_name="TP53-001",
        protein_id="ENSP00000269305",
        refseq="GGG",
        altseq="TTT",
    )
    result = variant.to_cdna()
    assert expected in result
    assert len(result) == 6


def test_to_dna(ensembl69, variant):
    expected = DnaDelins(
        _core=ensembl69,
        contig_id="17",
        start=7577058,
        start_offset=0,
        end=7577060,
        end_offset=0,
        strand="-",
        refseq="GGG",
        altseq="TTT",
    )
    result = variant.to_dna()
    assert expected in result
    assert len(result) == 1


def test_to_protein(ensembl69, variant):
    expected = ProteinDelins(
        _core=ensembl69,
        contig_id="17",
        start=293,
        start_offset=0,
        end=294,
        end_offset=0,
        strand="-",
        gene_id="ENSG00000141510",
        gene_name="TP53",
        transcript_id="ENST00000269305",
        transcript_name="TP53-001",
        protein_id="ENSP00000269305",
        refseq="GE",
        altseq="V*",
    )
    result = variant.to_protein()
    assert expected in result
    assert len(result) == 6


def test_to_rna(ensembl69, variant):
    expected = RnaDelins(
        _core=ensembl69,
        contig_id="17",
        start=1068,
        start_offset=0,
        end=1070,
        end_offset=0,
        strand="-",
        gene_id="ENSG00000141510",
        gene_name="TP53",
        transcript_id="ENST00000269305",
        transcript_name="TP53-001",
        refseq="GGG",
        altseq="TTT",
    )
    result = variant.to_rna()
    assert expected in result
    assert len(result) == 9
