import pytest

from variant_map.constants import INSERTION
from variant_map.positions import CdnaInsertion, DnaInsertion, ProteinInsertion, RnaInsertion


@pytest.fixture()
def variant():
    return ProteinInsertion(
        contig_id="12",
        start=11,
        start_offset=0,
        end=12,
        end_offset=0,
        strand="-",
        gene_id="ENSG00000133703",
        gene_name="KRAS",
        transcript_id="ENST00000256078",
        transcript_name="KRAS-004",
        protein_id="ENSP00000256078",
        refseq="AG",
        altseq="ATG",
    )


def test_str(variant):
    assert str(variant) == "ENSP00000256078:p.A11_G12insT"


def test_type(variant):
    assert not variant.is_cdna
    assert not variant.is_dna
    assert not variant.is_exon
    assert variant.is_protein
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
        contig_id="12",
        start=33,
        start_offset=0,
        end=34,
        end_offset=0,
        strand="-",
        gene_id="ENSG00000133703",
        gene_name="KRAS",
        transcript_id="ENST00000256078",
        transcript_name="KRAS-004",
        protein_id="ENSP00000256078",
        refseq="TG",
        altseq="TACAG",
    )
    result = ensembl69.to_cdna(variant)
    assert expected in result
    assert len(result) == 7


def test_to_dna(ensembl69, variant):
    expected = DnaInsertion(
        contig_id="12",
        start=25398285,
        start_offset=0,
        end=25398286,
        end_offset=0,
        strand="-",
        refseq="TG",
        altseq="TACAG",
    )
    result = ensembl69.to_dna(variant)
    assert expected in result
    assert len(result) == 7


def test_to_protein(ensembl69, variant):
    expected = [
        ProteinInsertion(
            contig_id="12",
            start=11,
            start_offset=0,
            end=12,
            end_offset=0,
            strand="-",
            gene_id="ENSG00000133703",
            gene_name="KRAS",
            transcript_id="ENST00000256078",
            transcript_name="KRAS-004",
            protein_id="ENSP00000256078",
            refseq="AG",
            altseq="ATG",
        )
    ]
    assert ensembl69.to_protein(variant) == expected


def test_to_rna(ensembl69, variant):
    expected = RnaInsertion(
        contig_id="12",
        start=97,
        start_offset=0,
        end=98,
        end_offset=0,
        strand="-",
        gene_id="ENSG00000133703",
        gene_name="KRAS",
        transcript_id="ENST00000256078",
        transcript_name="KRAS-004",
        refseq="TG",
        altseq="TACAG",
    )
    result = ensembl69.to_rna(variant)
    assert expected in result
    assert len(result) == 7
