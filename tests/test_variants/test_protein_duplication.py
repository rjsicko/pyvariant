import pytest

from pyvariant.constants import DUPLICATION
from pyvariant.variants import (
    CdnaDuplication,
    DnaDuplication,
    ExonSmallVariant,
    ProteinDuplication,
    RnaDuplication,
)


@pytest.fixture()
def variant(ensembl69):
    return ProteinDuplication(
        _core=ensembl69,
        contig_id="4",
        start=501,
        start_offset=0,
        end=501,
        end_offset=0,
        strand="+",
        gene_id="ENSG00000157404",
        gene_name="KIT",
        transcript_id="ENST00000288135",
        transcript_name="KIT-001",
        protein_id="ENSP00000288135",
        refseq="S",
        altseq="SS",
    )


def test_str(variant):
    assert str(variant) == "ENSP00000288135:p.S501dup"


def test_to_string_gene_name(variant):
    assert variant.to_string(reference="gene_name") == "KIT:p.S501dup"


def test_type(variant):
    assert not variant.is_cdna
    assert not variant.is_dna
    assert not variant.is_exon
    assert variant.is_protein
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
    expected = CdnaDuplication(
        _core=ensembl69,
        contig_id="4",
        start=1501,
        start_offset=0,
        end=1503,
        end_offset=0,
        strand="+",
        gene_id="ENSG00000157404",
        gene_name="KIT",
        transcript_id="ENST00000288135",
        transcript_name="KIT-001",
        protein_id="ENSP00000288135",
        refseq="TCT",
        altseq="TCTTCT",
    )
    result = variant.to_cdna()
    assert expected in result
    assert len(result) == 6


def test_to_dna(ensembl69, variant):
    expected = DnaDuplication(
        _core=ensembl69,
        contig_id="4",
        start=55592177,
        start_offset=0,
        end=55592179,
        end_offset=0,
        strand="+",
        refseq="TCT",
        altseq="TCTTCT",
    )
    result = variant.to_dna()
    assert expected in result
    assert len(result) == 6


# TODO: Reduce runtime for protein duplication mapping
def test_to_exon(ensembl69, variant):
    expected = ExonSmallVariant(
        _core=ensembl69,
        refseq="TCT",
        altseq="TCTTCT",
        contig_id="4",
        start=9,
        start_offset=0,
        end=9,
        end_offset=0,
        strand="+",
        gene_id="ENSG00000157404",
        gene_name="KIT",
        transcript_id="ENST00000288135",
        transcript_name="KIT-001",
        exon_id="ENSE00001074423",
    )
    result = variant.to_exon()
    assert expected in result
    assert len(result) == 6


def test_to_protein(ensembl69, variant):
    expected = [
        ProteinDuplication(
            _core=ensembl69,
            contig_id="4",
            start=501,
            start_offset=0,
            end=501,
            end_offset=0,
            strand="+",
            gene_id="ENSG00000157404",
            gene_name="KIT",
            transcript_id="ENST00000288135",
            transcript_name="KIT-001",
            protein_id="ENSP00000288135",
            refseq="S",
            altseq="SS",
        )
    ]
    assert variant.to_protein() == expected


def test_to_rna(ensembl69, variant):
    expected = RnaDuplication(
        _core=ensembl69,
        contig_id="4",
        start=1598,
        start_offset=0,
        end=1600,
        end_offset=0,
        strand="+",
        gene_id="ENSG00000157404",
        gene_name="KIT",
        transcript_id="ENST00000288135",
        transcript_name="KIT-001",
        refseq="TCT",
        altseq="TCTTCT",
    )
    result = variant.to_rna()
    assert expected in result
    assert len(result) == 6
