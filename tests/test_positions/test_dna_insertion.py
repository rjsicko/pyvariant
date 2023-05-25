import pytest

from pyvariant.constants import INSERTION
from pyvariant.positions import CdnaInsertion, DnaInsertion, ProteinInsertion, RnaInsertion


@pytest.fixture()
def variant(ensembl69):
    return DnaInsertion(
        _core=ensembl69,
        contig_id="17",
        start=7579470,
        start_offset=0,
        end=7579471,
        end_offset=0,
        strand="-",
        refseq="CG",
        altseq="CAAAG",
    )


def test_str(variant):
    assert str(variant) == "17:g.7579470_7579471insAAA"


def test_type(variant):
    assert not variant.is_cdna
    assert variant.is_dna
    assert not variant.is_exon
    assert not variant.is_protein
    assert not variant.is_rna


def test_variant_type(variant):
    assert variant.variant_type == INSERTION
    assert not variant.is_deletion
    assert not variant.is_delins
    assert not variant.is_duplication
    assert not variant.is_frameshift
    assert not variant.is_fusion
    assert variant.is_insertion
    assert not variant.is_substitution


def test_to_cdna(ensembl69, variant):
    expected = CdnaInsertion(
        _core=ensembl69,
        contig_id="17",
        start=216,
        start_offset=0,
        end=217,
        end_offset=0,
        strand="-",
        gene_id="ENSG00000141510",
        gene_name="TP53",
        transcript_id="ENST00000269305",
        transcript_name="TP53-001",
        protein_id="ENSP00000269305",
        refseq="CG",
        altseq="CAAAG",
    )
    result = ensembl69.to_cdna(variant)
    assert expected in result
    assert len(result) == 8


def test_to_dna(ensembl69, variant):
    expected = [
        DnaInsertion(
            _core=ensembl69,
            contig_id="17",
            start=7579470,
            start_offset=0,
            end=7579471,
            end_offset=0,
            strand="-",
            refseq="CG",
            altseq="CAAAG",
        )
    ]
    assert ensembl69.to_dna(variant) == expected


def test_to_protein(ensembl69, variant):
    expected = ProteinInsertion(
        _core=ensembl69,
        contig_id="17",
        start=72,
        start_offset=0,
        end=73,
        end_offset=0,
        strand="-",
        gene_id="ENSG00000141510",
        gene_name="TP53",
        transcript_id="ENST00000269305",
        transcript_name="TP53-001",
        protein_id="ENSP00000269305",
        refseq="PV",
        altseq="PKV",
    )
    result = ensembl69.to_protein(variant)
    assert expected in result
    assert len(result) == 8


def test_to_rna(ensembl69, variant):
    expected = RnaInsertion(
        _core=ensembl69,
        contig_id="17",
        start=406,
        start_offset=0,
        end=407,
        end_offset=0,
        strand="-",
        gene_id="ENSG00000141510",
        gene_name="TP53",
        transcript_id="ENST00000269305",
        transcript_name="TP53-001",
        refseq="CG",
        altseq="CAAAG",
    )
    result = ensembl69.to_rna(variant)
    assert expected in result
    assert len(result) == 9
