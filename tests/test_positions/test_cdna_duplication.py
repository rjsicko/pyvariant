import pytest

from variant_map.constants import DUPLICATION
from variant_map.positions import (
    CdnaDuplication,
    DnaDuplication,
    ProteinDuplication,
    RnaDuplication,
)


@pytest.fixture()
def variant():
    return CdnaDuplication(
        contig_id="17",
        start=286,
        start_offset=0,
        end=288,
        end_offset=0,
        strand="-",
        gene_id="ENSG00000141510",
        gene_name="TP53",
        transcript_id="ENST00000269305",
        transcript_name="TP53-001",
        protein_id="ENSP00000269305",
        refseq="TCT",
        altseq="TCTTCT",
    )


def test_str(variant):
    assert str(variant) == "ENST00000269305:c.286_288dup"


def test_type(variant):
    assert variant.is_cdna
    assert not variant.is_dna
    assert not variant.is_exon
    assert not variant.is_protein
    assert not variant.is_rna


def test_variant_type(variant):
    assert variant.variant_type == DUPLICATION
    assert not variant.is_deletion
    assert not variant.is_delins
    assert variant.is_duplication
    assert not variant.is_frameshift
    assert not variant.is_fusion
    assert not variant.is_insertion
    assert not variant.is_substitution


def test_to_cdna(ensembl69, variant):
    expected = [
        CdnaDuplication(
            contig_id="17",
            start=286,
            start_offset=0,
            end=288,
            end_offset=0,
            strand="-",
            gene_id="ENSG00000141510",
            gene_name="TP53",
            transcript_id="ENST00000269305",
            transcript_name="TP53-001",
            protein_id="ENSP00000269305",
            refseq="TCT",
            altseq="TCTTCT",
        )
    ]
    assert ensembl69.to_cdna(variant) == expected


def test_to_dna(ensembl69, variant):
    expected = [
        DnaDuplication(
            contig_id="17",
            start=7579399,
            start_offset=0,
            end=7579401,
            end_offset=0,
            strand="-",
            refseq="TCT",
            altseq="TCTTCT",
        )
    ]
    assert ensembl69.to_dna(variant) == expected


def test_to_protein(ensembl69, variant):
    expected = [
        ProteinDuplication(
            contig_id="17",
            start=96,
            start_offset=0,
            end=96,
            end_offset=0,
            strand="-",
            gene_id="ENSG00000141510",
            gene_name="TP53",
            transcript_id="ENST00000269305",
            transcript_name="TP53-001",
            protein_id="ENSP00000269305",
            refseq="S",
            altseq="SS",
        )
    ]
    assert ensembl69.to_protein(variant) == expected


def test_to_rna(ensembl69, variant):
    expected = [
        RnaDuplication(
            contig_id="17",
            start=476,
            start_offset=0,
            end=478,
            end_offset=0,
            strand="-",
            gene_id="ENSG00000141510",
            gene_name="TP53",
            transcript_id="ENST00000269305",
            transcript_name="TP53-001",
            refseq="TCT",
            altseq="TCTTCT",
        )
    ]
    assert ensembl69.to_rna(variant) == expected
