import pytest

from ensembl_map.constants import DELETION
from ensembl_map.core import CdnaDeletion, DnaDeletion, ProteinDeletion, RnaDeletion


@pytest.fixture()
def variant(ensembl69):
    return RnaDeletion(
        _data=ensembl69,
        contig_id="3",
        start=1318,
        start_offset=0,
        end=1320,
        end_offset=0,
        gene_id="ENSG00000134086",
        gene_name="VHL",
        transcript_id="ENST00000256474",
        transcript_name="VHL-001",
        strand="+",
        refseq="GAG",
        altseq="",
    )


def test_str(variant):
    assert str(variant) == "ENST00000256474:r.1318_1320del"


def test_type(variant):
    assert not variant.is_cdna
    assert not variant.is_dna
    assert not variant.is_exon
    assert not variant.is_protein
    assert variant.is_rna


def test_variant_type(variant):
    assert variant.type == DELETION
    assert variant.is_deletion
    assert not variant.is_delins
    assert not variant.is_duplication
    assert not variant.is_frameshift
    assert not variant.is_fusion
    assert not variant.is_insertion
    assert not variant.is_substitution


def test_to_cdna(ensembl69, variant):
    expected = [
        CdnaDeletion(
            _data=ensembl69,
            contig_id="3",
            start=478,
            start_offset=0,
            end=480,
            end_offset=0,
            gene_id="ENSG00000134086",
            gene_name="VHL",
            transcript_id="ENST00000256474",
            transcript_name="VHL-001",
            protein_id="ENSP00000256474",
            strand="+",
            refseq="GAG",
            altseq="",
        )
    ]
    assert variant.to_cdna() == expected


def test_to_dna(ensembl69, variant):
    expected = [
        DnaDeletion(
            _data=ensembl69,
            contig_id="3",
            start=10191485,
            start_offset=0,
            end=10191487,
            end_offset=0,
            strand="+",
            refseq="GAG",
            altseq="",
        )
    ]
    assert variant.to_dna() == expected


def test_to_protein(ensembl69, variant):
    expected = [
        ProteinDeletion(
            _data=ensembl69,
            contig_id="3",
            start=160,
            start_offset=0,
            end=160,
            end_offset=0,
            strand="+",
            gene_id="ENSG00000134086",
            gene_name="VHL",
            transcript_id="ENST00000256474",
            transcript_name="VHL-001",
            protein_id="ENSP00000256474",
            refseq="E",
            altseq="",
        )
    ]
    assert variant.to_protein() == expected


def test_to_rna(ensembl69, variant):
    expected = [
        RnaDeletion(
            _data=ensembl69,
            contig_id="3",
            start=1318,
            start_offset=0,
            end=1320,
            end_offset=0,
            strand="+",
            gene_id="ENSG00000134086",
            gene_name="VHL",
            transcript_id="ENST00000256474",
            transcript_name="VHL-001",
            refseq="GAG",
            altseq="",
        )
    ]
    assert variant.to_rna() == expected
