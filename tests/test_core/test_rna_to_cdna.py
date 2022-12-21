from ensembl_map.core import CdnaMappablePosition, RnaMappablePosition


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
    expected = CdnaMappablePosition(
        _data=ensembl100,
        contig_id="5",
        start=124,
        start_offset=0,
        end=126,
        end_offset=0,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
        protein_id="ENSP00000309572",
    )
    assert position.to_cdna() == [expected]


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
    expected = CdnaMappablePosition(
        _data=ensembl100,
        contig_id="5",
        start=1573,
        start_offset=0,
        end=1575,
        end_offset=0,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
        protein_id="ENSP00000309572",
    )
    assert position.to_cdna() == [expected]


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
    expected = CdnaMappablePosition(
        _data=ensembl100,
        contig_id="5",
        start=3397,
        start_offset=0,
        end=3399,
        end_offset=0,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
        protein_id="ENSP00000309572",
    )
    assert position.to_cdna() == [expected]


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
    expected = CdnaMappablePosition(
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
        protein_id="ENSP00000309572",
    )
    assert position.to_cdna() == [expected]


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
    expected = CdnaMappablePosition(
        _data=ensembl100,
        contig_id="13",
        start=8566,
        start_offset=0,
        end=8568,
        end_offset=0,
        strand="+",
        gene_id="ENSG00000139618",
        gene_name="BRCA2",
        transcript_id="ENST00000380152",
        transcript_name="BRCA2-201",
        protein_id="ENSP00000369497",
    )
    assert position.to_cdna() == [expected]


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
    expected = CdnaMappablePosition(
        _data=ensembl100,
        contig_id="4",
        start=67,
        start_offset=0,
        end=69,
        end_offset=0,
        strand="+",
        gene_id="ENSG00000157404",
        gene_name="KIT",
        transcript_id="ENST00000288135",
        transcript_name="KIT-201",
        protein_id="ENSP00000288135",
    )
    assert position.to_cdna() == [expected]


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
    expected = CdnaMappablePosition(
        _data=ensembl100,
        contig_id="13",
        start=10255,
        start_offset=0,
        end=10257,
        end_offset=0,
        strand="+",
        gene_id="ENSG00000139618",
        gene_name="BRCA2",
        transcript_id="ENST00000380152",
        transcript_name="BRCA2-201",
        protein_id="ENSP00000369497",
    )
    assert position.to_cdna() == [expected]


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
    expected = CdnaMappablePosition(
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
        protein_id="ENSP00000369497",
    )
    assert position.to_cdna() == [expected]


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
    expected = CdnaMappablePosition(
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
        protein_id="ENSP00000439397",
    )
    assert position.to_cdna() == [expected]


def test_offset(ensembl100):
    position = RnaMappablePosition(
        _data=ensembl100,
        contig_id="5",
        start=299,
        start_offset=-1,
        end=299,
        end_offset=-1,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
    )
    expected = CdnaMappablePosition(
        _data=ensembl100,
        contig_id="5",
        start=220,
        start_offset=-1,
        end=220,
        end_offset=-1,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
        protein_id="ENSP00000309572",
    )
    assert position.to_cdna() == [expected]


def test_offset_normalized(ensembl100):
    position = RnaMappablePosition(
        _data=ensembl100,
        contig_id="5",
        start=297,
        start_offset=-1,
        end=297,
        end_offset=-1,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
    )
    expected = CdnaMappablePosition(
        _data=ensembl100,
        contig_id="5",
        start=217,
        start_offset=0,
        end=217,
        end_offset=0,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
        protein_id="ENSP00000309572",
    )
    assert position.to_cdna() == [expected]
