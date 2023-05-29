import pytest

from pyvariant.constants import DUPLICATION
from pyvariant.positions import CdnaDuplication, DnaDuplication, ProteinDuplication, RnaDuplication


@pytest.fixture()
def variant(ensembl69):
    return CdnaDuplication(
        _core=ensembl69,
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


def test_to_string_gene_name(variant):
    assert variant.to_string(reference="gene_name") == "TP53:c.286_288dup"


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
            _core=ensembl69,
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
    assert variant.to_cdna() == expected


def test_to_dna(ensembl69, variant):
    expected = [
        DnaDuplication(
            _core=ensembl69,
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
    assert variant.to_dna() == expected


def test_to_protein(ensembl69, variant):
    expected = [
        ProteinDuplication(
            _core=ensembl69,
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
    assert variant.to_protein() == expected


def test_to_rna(ensembl69, variant):
    expected = [
        RnaDuplication(
            _core=ensembl69,
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
    assert variant.to_rna() == expected
