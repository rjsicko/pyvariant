import pytest

from variant_map.constants import FUSION
from variant_map.positions import (
    CdnaFusion,
    CdnaPosition,
    DnaFusion,
    DnaPosition,
    ExonFusion,
    ExonPosition,
    ProteinFusion,
    ProteinPosition,
    RnaFusion,
    RnaPosition,
)


@pytest.fixture()
def variant(ensembl69):
    breakpoint1 = DnaPosition(
        _data=ensembl69,
        contig_id="X",
        start=48895958,
        start_offset=0,
        end=48896632,
        end_offset=0,
        strand="-",
    )
    breakpoint2 = DnaPosition(
        _data=ensembl69,
        contig_id="1",
        start=156752074,
        start_offset=0,
        end=156752084,
        end_offset=0,
        strand="+",
    )
    return DnaFusion(breakpoint1, breakpoint2)


def test_str(variant):
    assert str(variant) == "X:g.48895958_48896632::1:g.156752074_156752084"


def test_type(variant):
    assert not variant.is_cdna
    assert variant.is_dna
    assert not variant.is_exon
    assert not variant.is_protein
    assert not variant.is_rna


def test_variant_type(variant):
    assert variant.type == FUSION
    assert not variant.is_deletion
    assert not variant.is_delins
    assert not variant.is_duplication
    assert not variant.is_frameshift
    assert variant.is_fusion
    assert not variant.is_insertion
    assert not variant.is_substitution


def test_to_cdna(ensembl69, variant):
    breakpoint1 = CdnaPosition(
        _data=ensembl69,
        contig_id="X",
        start=534,
        start_offset=0,
        end=544,
        end_offset=0,
        strand="-",
        gene_id="ENSG00000068323",
        gene_name="TFE3",
        transcript_id="ENST00000315869",
        transcript_name="TFE3-001",
        protein_id="ENSP00000314129",
    )
    breakpoint2 = CdnaPosition(
        _data=ensembl69,
        contig_id="1",
        start=469,
        start_offset=0,
        end=479,
        end_offset=0,
        strand="+",
        gene_id="ENSG00000143294",
        gene_name="PRCC",
        transcript_id="ENST00000271526",
        transcript_name="PRCC-001",
        protein_id="ENSP00000271526",
    )
    expected = CdnaFusion(breakpoint1, breakpoint2)
    result = variant.to_cdna()
    assert expected in result
    assert len(result) == 2


def test_to_dna(ensembl69, variant):
    breakpoint1 = DnaPosition(
        _data=ensembl69,
        contig_id="X",
        start=48895958,
        start_offset=0,
        end=48896632,
        end_offset=0,
        strand="-",
    )
    breakpoint2 = DnaPosition(
        _data=ensembl69,
        contig_id="1",
        start=156752074,
        start_offset=0,
        end=156752084,
        end_offset=0,
        strand="+",
    )
    expected = DnaFusion(breakpoint1, breakpoint2)
    result = variant.to_dna()
    assert expected in result
    assert len(result) == 1


def test_to_exon(ensembl69, variant):
    breakpoint1 = ExonPosition(
        _data=ensembl69,
        contig_id="X",
        start=3,
        start_offset=0,
        end=4,
        end_offset=0,
        strand="-",
        gene_id="ENSG00000068323",
        gene_name="TFE3",
        transcript_id="ENST00000315869",
        transcript_name="TFE3-001",
        exon_id="ENSE00002939607",
    )
    breakpoint2 = ExonPosition(
        _data=ensembl69,
        contig_id="1",
        start=2,
        start_offset=0,
        end=2,
        end_offset=0,
        strand="+",
        gene_id="ENSG00000143294",
        gene_name="PRCC",
        transcript_id="ENST00000271526",
        transcript_name="PRCC-001",
        exon_id="ENSE00002871917",
    )
    expected = ExonFusion(breakpoint1, breakpoint2)
    result = variant.to_exon()
    assert expected in result
    assert len(result) == 16


def test_to_protein(ensembl69, variant):
    breakpoint1 = ProteinPosition(
        _data=ensembl69,
        contig_id="X",
        start=178,
        start_offset=0,
        end=182,
        end_offset=0,
        strand="-",
        gene_id="ENSG00000068323",
        gene_name="TFE3",
        transcript_id="ENST00000315869",
        transcript_name="TFE3-001",
        protein_id="ENSP00000314129",
    )
    breakpoint2 = ProteinPosition(
        _data=ensembl69,
        contig_id="1",
        start=157,
        start_offset=0,
        end=160,
        end_offset=0,
        strand="+",
        gene_id="ENSG00000143294",
        gene_name="PRCC",
        transcript_id="ENST00000271526",
        transcript_name="PRCC-001",
        protein_id="ENSP00000271526",
    )
    expected = ProteinFusion(breakpoint1, breakpoint2)
    result = variant.to_protein()
    assert expected in result
    assert len(result) == 2


def test_to_rna(ensembl69, variant):
    breakpoint1 = RnaPosition(
        _data=ensembl69,
        contig_id="X",
        start=794,
        start_offset=0,
        end=804,
        end_offset=0,
        strand="-",
        gene_id="ENSG00000068323",
        gene_name="TFE3",
        transcript_id="ENST00000315869",
        transcript_name="TFE3-001",
    )
    breakpoint2 = RnaPosition(
        _data=ensembl69,
        contig_id="1",
        start=741,
        start_offset=0,
        end=751,
        end_offset=0,
        strand="+",
        gene_id="ENSG00000143294",
        gene_name="PRCC",
        transcript_id="ENST00000271526",
        transcript_name="PRCC-001",
    )
    expected = RnaFusion(breakpoint1, breakpoint2)
    result = variant.to_rna()
    assert expected in result
    assert len(result) == 16
