from ensembl_map.core import RnaMappablePosition


def test_negative_strand(ensembl100):
    position = RnaMappablePosition(
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
    assert position.to_rna() == [position]


def test_negative_strand_across_exon_boundary(ensembl100):
    position = RnaMappablePosition(
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
    assert position.to_rna() == [position]


def test_negative_strand_cdna_protein_end(ensembl100):
    position = RnaMappablePosition(
        _data=ensembl100,
        contig_id="5",
        start=3476,
        start_offset=0,
        end=3478,
        end_offset=0,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
    )
    assert position.to_rna() == [position]


def test_negative_strand_cdna_protein_start(ensembl100):
    position = RnaMappablePosition(
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
    assert position.to_rna() == [position]


def test_negative_strand_overlapping_genes_different_strands(ensembl100):
    position = RnaMappablePosition(
        _data=ensembl100,
        contig_id="6",
        start=535,
        start_offset=0,
        end=537,
        end_offset=0,
        strand="-",
        gene_id="ENSG00000204531",
        gene_name="POU5F1",
        transcript_id="ENST00000471529",
        transcript_name="POU5F1-204",
    )
    assert position.to_rna() == [position]


def test_negative_strand_transcript_end(ensembl100):
    position = RnaMappablePosition(
        _data=ensembl100,
        contig_id="5",
        start=4037,
        start_offset=0,
        end=4039,
        end_offset=0,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
    )
    assert position.to_rna() == [position]


def test_negative_strand_transcript_start(ensembl100):
    position = RnaMappablePosition(
        _data=ensembl100,
        contig_id="5",
        start=1,
        start_offset=0,
        end=3,
        end_offset=0,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
    )
    assert position.to_rna() == [position]


def test_positive_strand(ensembl100):
    position = RnaMappablePosition(
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
    assert position.to_rna() == [position]


def test_positive_strand_across_exon_boundary(ensembl100):
    position = RnaMappablePosition(
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
    assert position.to_rna() == [position]


def test_positive_strand_cdna_protein_end(ensembl100):
    position = RnaMappablePosition(
        _data=ensembl100,
        contig_id="13",
        start=10488,
        start_offset=0,
        end=10490,
        end_offset=0,
        strand="+",
        gene_id="ENSG00000139618",
        gene_name="BRCA2",
        transcript_id="ENST00000380152",
        transcript_name="BRCA2-201",
    )
    assert position.to_rna() == [position]


def test_positive_strand_cdna_protein_start(ensembl100):
    position = RnaMappablePosition(
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
    assert position.to_rna() == [position]


def test_positive_strand_overlapping_genes_different_strands(ensembl100):
    position = RnaMappablePosition(
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
    assert position.to_rna() == [position]


def test_positive_strand_transcript_end(ensembl100):
    position = RnaMappablePosition(
        _data=ensembl100,
        contig_id="13",
        start=11984,
        start_offset=0,
        end=11986,
        end_offset=0,
        strand="+",
        gene_id="ENSG00000139618",
        gene_name="BRCA2",
        transcript_id="ENST00000380152",
        transcript_name="BRCA2-201",
    )
    assert position.to_rna() == [position]


def test_positive_strand_transcript_start(ensembl100):
    position = RnaMappablePosition(
        _data=ensembl100,
        contig_id="13",
        start=1,
        start_offset=0,
        end=3,
        end_offset=0,
        strand="+",
        gene_id="ENSG00000139618",
        gene_name="BRCA2",
        transcript_id="ENST00000380152",
        transcript_name="BRCA2-201",
    )
    assert position.to_rna() == [position]


def test_offset(ensembl100):
    position = RnaMappablePosition(
        _data=ensembl100,
        contig_id="5",
        start=298,
        start_offset=1,
        end=298,
        end_offset=1,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
    )
    assert position.to_rna() == [position]


def test_offset_normalized(ensembl100):
    position = RnaMappablePosition(
        _data=ensembl100,
        contig_id="5",
        start=100,
        start_offset=-1,
        end=100,
        end_offset=-1,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
    )
    expected = RnaMappablePosition(
        _data=ensembl100,
        contig_id="5",
        start=99,
        start_offset=0,
        end=99,
        end_offset=0,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
    )
    assert position.to_rna() == [expected]
