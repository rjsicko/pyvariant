import pytest

from variant_map.core import ProteinMappablePosition, RnaMappablePosition


def test_negative_strand(ensembl100):
    position = ProteinMappablePosition(
        _data=ensembl100,
        contig_id="5",
        start=42,
        start_offset=0,
        end=42,
        end_offset=0,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
        protein_id="ENSP00000309572",
    )
    expected = RnaMappablePosition(
        _data=ensembl100,
        contig_id="5",
        start=203,
        start_offset=0,
        end=205,
        end_offset=0,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
    )
    assert position.to_rna() == [expected]


def test_negative_strand_across_exon_boundary(ensembl100):
    position = ProteinMappablePosition(
        _data=ensembl100,
        contig_id="5",
        start=525,
        start_offset=0,
        end=525,
        end_offset=0,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
        protein_id="ENSP00000309572",
    )
    expected = RnaMappablePosition(
        _data=ensembl100,
        contig_id="5",
        start=1652,
        start_offset=0,
        end=1654,
        end_offset=0,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
    )
    assert position.to_rna() == [expected]


def test_negative_strand_cdna_protein_end(ensembl100):
    position = ProteinMappablePosition(
        _data=ensembl100,
        contig_id="5",
        start=1132,
        start_offset=0,
        end=1132,
        end_offset=0,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
        protein_id="ENSP00000309572",
    )
    expected = RnaMappablePosition(
        _data=ensembl100,
        contig_id="5",
        start=3473,
        start_offset=0,
        end=3475,
        end_offset=0,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
    )
    assert position.to_rna() == [expected]


def test_negative_strand_cdna_protein_start(ensembl100):
    position = ProteinMappablePosition(
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
        protein_id="ENSP00000309572",
    )
    expected = RnaMappablePosition(
        _data=ensembl100,
        contig_id="5",
        start=80,
        start_offset=0,
        end=82,
        end_offset=0,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
    )
    assert position.to_rna() == [expected]


def test_positive_strand(ensembl100):
    position = ProteinMappablePosition(
        _data=ensembl100,
        contig_id="13",
        start=2856,
        start_offset=0,
        end=2856,
        end_offset=0,
        strand="+",
        gene_id="ENSG00000139618",
        gene_name="BRCA2",
        transcript_id="ENST00000380152",
        transcript_name="BRCA2-201",
        protein_id="ENSP00000369497",
    )
    expected = RnaMappablePosition(
        _data=ensembl100,
        contig_id="13",
        start=8799,
        start_offset=0,
        end=8801,
        end_offset=0,
        strand="+",
        gene_id="ENSG00000139618",
        gene_name="BRCA2",
        transcript_id="ENST00000380152",
        transcript_name="BRCA2-201",
    )
    assert position.to_rna() == [expected]


def test_positive_strand_across_exon_boundary(ensembl100):
    position = ProteinMappablePosition(
        _data=ensembl100,
        contig_id="4",
        start=23,
        start_offset=0,
        end=23,
        end_offset=0,
        strand="+",
        gene_id="ENSG00000157404",
        gene_name="KIT",
        transcript_id="ENST00000288135",
        transcript_name="KIT-201",
        protein_id="ENSP00000288135",
    )
    expected = RnaMappablePosition(
        _data=ensembl100,
        contig_id="4",
        start=125,
        start_offset=0,
        end=127,
        end_offset=0,
        strand="+",
        gene_id="ENSG00000157404",
        gene_name="KIT",
        transcript_id="ENST00000288135",
        transcript_name="KIT-201",
    )
    assert position.to_rna() == [expected]


def test_positive_strand_cdna_protein_end(ensembl100):
    position = ProteinMappablePosition(
        _data=ensembl100,
        contig_id="13",
        start=3418,
        start_offset=0,
        end=3418,
        end_offset=0,
        strand="+",
        gene_id="ENSG00000139618",
        gene_name="BRCA2",
        transcript_id="ENST00000380152",
        transcript_name="BRCA2-201",
        protein_id="ENSP00000369497",
    )
    expected = RnaMappablePosition(
        _data=ensembl100,
        contig_id="13",
        start=10485,
        start_offset=0,
        end=10487,
        end_offset=0,
        strand="+",
        gene_id="ENSG00000139618",
        gene_name="BRCA2",
        transcript_id="ENST00000380152",
        transcript_name="BRCA2-201",
    )
    assert position.to_rna() == [expected]


def test_positive_strand_cdna_protein_start(ensembl100):
    position = ProteinMappablePosition(
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
        protein_id="ENSP00000369497",
    )
    expected = RnaMappablePosition(
        _data=ensembl100,
        contig_id="13",
        start=234,
        start_offset=0,
        end=236,
        end_offset=0,
        strand="+",
        gene_id="ENSG00000139618",
        gene_name="BRCA2",
        transcript_id="ENST00000380152",
        transcript_name="BRCA2-201",
    )
    assert position.to_rna() == [expected]


def test_positive_strand_overlapping_genes_different_strands(ensembl100):
    position = ProteinMappablePosition(
        _data=ensembl100,
        contig_id="6",
        start=194,
        start_offset=0,
        end=194,
        end_offset=0,
        strand="+",
        gene_id="ENSG00000137310",
        gene_name="TCF19",
        transcript_id="ENST00000542218",
        transcript_name="TCF19-204",
        protein_id="ENSP00000439397",
    )
    expected = RnaMappablePosition(
        _data=ensembl100,
        contig_id="6",
        start=580,
        start_offset=0,
        end=582,
        end_offset=0,
        strand="+",
        gene_id="ENSG00000137310",
        gene_name="TCF19",
        transcript_id="ENST00000542218",
        transcript_name="TCF19-204",
    )
    assert position.to_rna() == [expected]


def test_offset_error(ensembl100):
    position = ProteinMappablePosition(
        _data=ensembl100,
        contig_id="4",
        start=21,
        start_offset=1,
        end=21,
        end_offset=1,
        strand="+",
        gene_id="ENSG00000157404",
        gene_name="KIT",
        transcript_id="ENST00000288135",
        transcript_name="KIT-201",
        protein_id="ENSP00000288135",
    )
    with pytest.raises(AssertionError):
        position.to_rna()
