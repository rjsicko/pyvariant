import pytest

from ensembl_map.core import ExonMappablePosition


def test_negative_strand(ensembl100):
    position = ExonMappablePosition(
        _data=ensembl100,
        contig_id="5",
        start=2,
        start_offset=0,
        end=2,
        end_offset=0,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
        exon_id="ENSE00001197112",
    )
    assert position.to_exon() == [position]


def test_negative_strand_cdna_protein_end(ensembl100):
    position = ExonMappablePosition(
        _data=ensembl100,
        contig_id="5",
        start=16,
        start_offset=0,
        end=16,
        end_offset=0,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
        exon_id="ENSE00001863787",
    )
    assert position.to_exon() == [position]


def test_negative_strand_cdna_protein_start(ensembl100):
    position = ExonMappablePosition(
        _data=ensembl100,
        contig_id="5",
        start=1,
        start_offset=0,
        end=1,
        end_offset=0,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
        exon_id="ENSE00003896691",
    )
    assert position.to_exon() == [position]


def test_negative_strand_overlapping_genes_different_strands(ensembl100):
    position = ExonMappablePosition(
        _data=ensembl100,
        contig_id="6",
        start=1,
        start_offset=0,
        end=1,
        end_offset=0,
        strand="-",
        gene_id="ENSG00000204531",
        gene_name="POU5F1",
        transcript_id="ENST00000471529",
        transcript_name="POU5F1-204",
        exon_id="ENSE00002568331",
    )
    assert position.to_exon() == [position]


def test_negative_strand_transcript_end(ensembl100):
    position = ExonMappablePosition(
        _data=ensembl100,
        contig_id="5",
        start=16,
        start_offset=0,
        end=16,
        end_offset=0,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
        exon_id="ENSE00001863787",
    )
    assert position.to_exon() == [position]


def test_negative_strand_transcript_start(ensembl100):
    position = ExonMappablePosition(
        _data=ensembl100,
        contig_id="5",
        start=1,
        start_offset=0,
        end=1,
        end_offset=0,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
        exon_id="ENSE00003896691",
    )
    assert position.to_exon() == [position]


def test_positive_strand(ensembl100):
    position = ExonMappablePosition(
        _data=ensembl100,
        contig_id="13",
        start=20,
        start_offset=0,
        end=20,
        end_offset=0,
        strand="+",
        gene_id="ENSG00000139618",
        gene_name="BRCA2",
        transcript_id="ENST00000380152",
        transcript_name="BRCA2-201",
        exon_id="ENSE00000939180",
    )
    assert position.to_exon() == [position]


def test_positive_strand_cdna_protein_end(ensembl100):
    position = ExonMappablePosition(
        _data=ensembl100,
        contig_id="13",
        start=27,
        start_offset=0,
        end=27,
        end_offset=0,
        strand="+",
        gene_id="ENSG00000139618",
        gene_name="BRCA2",
        transcript_id="ENST00000380152",
        transcript_name="BRCA2-201",
        exon_id="ENSE00003717596",
    )
    assert position.to_exon() == [position]


def test_positive_strand_cdna_protein_start(ensembl100):
    position = ExonMappablePosition(
        _data=ensembl100,
        contig_id="13",
        start=2,
        start_offset=0,
        end=2,
        end_offset=0,
        strand="+",
        gene_id="ENSG00000139618",
        gene_name="BRCA2",
        transcript_id="ENST00000380152",
        transcript_name="BRCA2-201",
        exon_id="ENSE00001484009",
    )
    assert position.to_exon() == [position]


def test_positive_strand_overlapping_genes_different_strands(ensembl100):
    position = ExonMappablePosition(
        _data=ensembl100,
        contig_id="6",
        start=2,
        start_offset=0,
        end=2,
        end_offset=0,
        strand="+",
        gene_id="ENSG00000137310",
        gene_name="TCF19",
        transcript_id="ENST00000542218",
        transcript_name="TCF19-204",
        exon_id="ENSE00002283659",
    )
    assert position.to_exon() == [position]


def test_positive_strand_transcript_end(ensembl100):
    position = ExonMappablePosition(
        _data=ensembl100,
        contig_id="13",
        start=27,
        start_offset=0,
        end=27,
        end_offset=0,
        strand="+",
        gene_id="ENSG00000139618",
        gene_name="BRCA2",
        transcript_id="ENST00000380152",
        transcript_name="BRCA2-201",
        exon_id="ENSE00003717596",
    )
    assert position.to_exon() == [position]


def test_positive_strand_transcript_start(ensembl100):
    position = ExonMappablePosition(
        _data=ensembl100,
        contig_id="13",
        start=1,
        start_offset=0,
        end=1,
        end_offset=0,
        strand="+",
        gene_id="ENSG00000139618",
        gene_name="BRCA2",
        transcript_id="ENST00000380152",
        transcript_name="BRCA2-201",
        exon_id="ENSE00001184784",
    )
    assert position.to_exon() == [position]


def test_offset_error(ensembl100):
    position = ExonMappablePosition(
        _data=ensembl100,
        contig_id="5",
        start=1,
        start_offset=1,
        end=1,
        end_offset=1,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
        exon_id="ENSE00001197112",
    )
    with pytest.raises(AssertionError):
        position.to_exon()
