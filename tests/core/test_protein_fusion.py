import pytest

from variant_map.constants import FUSION
from variant_map.core import (
    CdnaFusion,
    CdnaMappablePosition,
    DnaFusion,
    DnaMappablePosition,
    ExonFusion,
    ExonMappablePosition,
    ProteinFusion,
    ProteinMappablePosition,
    RnaFusion,
    RnaMappablePosition,
)


@pytest.fixture()
def variant(ensembl69):
    breakpoint1 = ProteinMappablePosition(
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
    breakpoint2 = ProteinMappablePosition(
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
    return ProteinFusion(breakpoint1, breakpoint2)


def test_str(variant):
    assert str(variant) == "ENSP00000314129:p.178_182::ENSP00000271526:p.157_160"


def test_type(variant):
    assert not variant.is_cdna
    assert not variant.is_dna
    assert not variant.is_exon
    assert variant.is_protein
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
    breakpoint1 = CdnaMappablePosition(
        _data=ensembl69,
        contig_id="X",
        start=532,
        start_offset=0,
        end=546,
        end_offset=0,
        strand="-",
        gene_id="ENSG00000068323",
        gene_name="TFE3",
        transcript_id="ENST00000315869",
        transcript_name="TFE3-001",
        protein_id="ENSP00000314129",
    )
    breakpoint2 = CdnaMappablePosition(
        _data=ensembl69,
        contig_id="1",
        start=469,
        start_offset=0,
        end=480,
        end_offset=0,
        strand="+",
        gene_id="ENSG00000143294",
        gene_name="PRCC",
        transcript_id="ENST00000271526",
        transcript_name="PRCC-001",
        protein_id="ENSP00000271526",
    )
    expected = [CdnaFusion(breakpoint1, breakpoint2)]
    assert variant.to_cdna() == expected


def test_to_dna(ensembl69, variant):
    breakpoint1 = DnaMappablePosition(
        _data=ensembl69,
        contig_id="X",
        start=48895956,
        start_offset=0,
        end=48896634,
        end_offset=0,
        strand="-",
    )
    breakpoint2 = DnaMappablePosition(
        _data=ensembl69,
        contig_id="1",
        start=156752074,
        start_offset=0,
        end=156752085,
        end_offset=0,
        strand="+",
    )
    expected = [DnaFusion(breakpoint1, breakpoint2)]
    assert variant.to_dna() == expected


def test_to_exon(ensembl69, variant):
    breakpoint1 = ExonMappablePosition(
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
    breakpoint2 = ExonMappablePosition(
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
    expected = [ExonFusion(breakpoint1, breakpoint2)]
    assert variant.to_exon() == expected


def test_to_protein(ensembl69, variant):
    breakpoint1 = ProteinMappablePosition(
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
    breakpoint2 = ProteinMappablePosition(
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
    expected = [ProteinFusion(breakpoint1, breakpoint2)]
    assert variant.to_protein() == expected


def test_to_rna(ensembl69, variant):
    breakpoint1 = RnaMappablePosition(
        _data=ensembl69,
        contig_id="X",
        start=792,
        start_offset=0,
        end=806,
        end_offset=0,
        strand="-",
        gene_id="ENSG00000068323",
        gene_name="TFE3",
        transcript_id="ENST00000315869",
        transcript_name="TFE3-001",
    )
    breakpoint2 = RnaMappablePosition(
        _data=ensembl69,
        contig_id="1",
        start=741,
        start_offset=0,
        end=752,
        end_offset=0,
        strand="+",
        gene_id="ENSG00000143294",
        gene_name="PRCC",
        transcript_id="ENST00000271526",
        transcript_name="PRCC-001",
    )
    expected = [RnaFusion(breakpoint1, breakpoint2)]
    assert variant.to_rna() == expected
