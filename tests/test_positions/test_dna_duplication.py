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
    return DnaDuplication(
        contig_id="4",
        start=55570013,
        start_offset=0,
        end=55570015,
        end_offset=0,
        strand="+",
        refseq="AAT",
        altseq="AATAAT",
    )


def test_str(variant):
    assert str(variant) == "4:g.55570013_55570015dup"


def test_type(variant):
    assert not variant.is_cdna
    assert variant.is_dna
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
            contig_id="4",
            start=880,
            start_offset=0,
            end=882,
            end_offset=0,
            strand="+",
            gene_id="ENSG00000157404",
            gene_name="KIT",
            transcript_id="ENST00000288135",
            transcript_name="KIT-001",
            protein_id="ENSP00000288135",
            refseq="AAT",
            altseq="AATAAT",
        ),
        CdnaDuplication(
            contig_id="4",
            start=880,
            start_offset=0,
            end=882,
            end_offset=0,
            strand="+",
            gene_id="ENSG00000157404",
            gene_name="KIT",
            transcript_id="ENST00000412167",
            transcript_name="KIT-002",
            protein_id="ENSP00000390987",
            refseq="AAT",
            altseq="AATAAT",
        ),
    ]
    assert ensembl69.to_cdna(variant) == expected


def test_to_dna(ensembl69, variant):
    expected = [
        DnaDuplication(
            contig_id="4",
            start=55570013,
            start_offset=0,
            end=55570015,
            end_offset=0,
            strand="+",
            refseq="AAT",
            altseq="AATAAT",
        )
    ]
    assert ensembl69.to_dna(variant) == expected


def test_to_protein(ensembl69, variant):
    expected = [
        ProteinDuplication(
            contig_id="4",
            start=294,
            start_offset=0,
            end=294,
            end_offset=0,
            strand="+",
            gene_id="ENSG00000157404",
            gene_name="KIT",
            transcript_id="ENST00000288135",
            transcript_name="KIT-001",
            protein_id="ENSP00000288135",
            refseq="N",
            altseq="NN",
        ),
        ProteinDuplication(
            contig_id="4",
            start=294,
            start_offset=0,
            end=294,
            end_offset=0,
            strand="+",
            gene_id="ENSG00000157404",
            gene_name="KIT",
            transcript_id="ENST00000412167",
            transcript_name="KIT-002",
            protein_id="ENSP00000390987",
            refseq="N",
            altseq="NN",
        ),
    ]
    assert ensembl69.to_protein(variant) == expected


def test_to_rna(ensembl69, variant):
    expected = [
        RnaDuplication(
            contig_id="4",
            start=977,
            start_offset=0,
            end=979,
            end_offset=0,
            strand="+",
            gene_id="ENSG00000157404",
            gene_name="KIT",
            transcript_id="ENST00000288135",
            transcript_name="KIT-001",
            refseq="AAT",
            altseq="AATAAT",
        ),
        RnaDuplication(
            contig_id="4",
            start=977,
            start_offset=0,
            end=979,
            end_offset=0,
            strand="+",
            gene_id="ENSG00000157404",
            gene_name="KIT",
            transcript_id="ENST00000412167",
            transcript_name="KIT-002",
            refseq="AAT",
            altseq="AATAAT",
        ),
    ]
    assert ensembl69.to_rna(variant) == expected
