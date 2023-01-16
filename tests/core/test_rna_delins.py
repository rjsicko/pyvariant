import pytest

from variant_map.constants import DELINS
from variant_map.core import CdnaDelins, DnaDelins, ProteinDelins, RnaDelins


@pytest.fixture()
def variant(ensembl69):
    return RnaDelins(
        _data=ensembl69,
        contig_id="4",
        start=976,
        start_offset=0,
        end=977,
        end_offset=0,
        strand="+",
        gene_id="ENSG00000157404",
        gene_name="KIT",
        transcript_id="ENST00000288135",
        transcript_name="KIT-001",
        refseq="TA",
        altseq="AG",
    )


def test_str(variant):
    assert str(variant) == "ENST00000288135:r.976_977delinsAG"


def test_type(variant):
    assert not variant.is_cdna
    assert not variant.is_dna
    assert not variant.is_exon
    assert not variant.is_protein
    assert variant.is_rna


def test_variant_type(variant):
    assert variant.type == DELINS
    assert not variant.is_deletion
    assert variant.is_delins
    assert not variant.is_duplication
    assert not variant.is_frameshift
    assert not variant.is_fusion
    assert not variant.is_insertion
    assert not variant.is_substitution


def test_to_cdna(ensembl69, variant):
    expected = [
        CdnaDelins(
            _data=ensembl69,
            contig_id="4",
            start=879,
            start_offset=0,
            end=880,
            end_offset=0,
            strand="+",
            gene_id="ENSG00000157404",
            gene_name="KIT",
            transcript_id="ENST00000288135",
            transcript_name="KIT-001",
            protein_id="ENSP00000288135",
            refseq="TA",
            altseq="AG",
        )
    ]
    assert variant.to_cdna() == expected


def test_to_dna(ensembl69, variant):
    expected = [
        DnaDelins(
            _data=ensembl69,
            contig_id="4",
            start=55570012,
            start_offset=0,
            end=55570013,
            end_offset=0,
            strand="+",
            refseq="TA",
            altseq="AG",
        )
    ]
    assert variant.to_dna() == expected


def test_to_protein(ensembl69, variant):
    expected = [
        ProteinDelins(
            _data=ensembl69,
            contig_id="4",
            start=293,
            start_offset=0,
            end=294,
            end_offset=0,
            strand="+",
            gene_id="ENSG00000157404",
            gene_name="KIT",
            transcript_id="ENST00000288135",
            transcript_name="KIT-001",
            protein_id="ENSP00000288135",
            refseq="NN",
            altseq="KD",
        )
    ]
    assert variant.to_protein() == expected


def test_to_rna(ensembl69, variant):
    expected = [
        RnaDelins(
            _data=ensembl69,
            contig_id="4",
            start=976,
            start_offset=0,
            end=977,
            end_offset=0,
            strand="+",
            gene_id="ENSG00000157404",
            gene_name="KIT",
            transcript_id="ENST00000288135",
            transcript_name="KIT-001",
            refseq="TA",
            altseq="AG",
        )
    ]
    assert variant.to_rna() == expected
