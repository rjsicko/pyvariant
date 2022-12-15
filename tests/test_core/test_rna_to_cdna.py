from ensembl_map.core import CdnaPosition, RnaPosition


def test_negative_strand(ensembl100):
    position = RnaPosition(
        _data=ensembl100,
        contig_id="5",
        start=203,
        end=205,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
    )
    expected = CdnaPosition(
        _data=ensembl100,
        contig_id="5",
        start=124,
        end=126,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
        protein_id="ENSP00000309572",
    )
    assert position.to_cdna() == [expected]


def test_negative_strand_across_exon_boundary(ensembl100):
    position = RnaPosition(
        _data=ensembl100,
        contig_id="5",
        start=1652,
        end=1654,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
    )
    expected = CdnaPosition(
        _data=ensembl100,
        contig_id="5",
        start=1573,
        end=1575,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
        protein_id="ENSP00000309572",
    )
    assert position.to_cdna() == [expected]


def test_negative_strand_cdna_protein_end(ensembl100):
    position = RnaPosition(
        _data=ensembl100,
        contig_id="5",
        start=3476,
        end=3478,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
    )
    expected = CdnaPosition(
        _data=ensembl100,
        contig_id="5",
        start=3397,
        end=3399,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
        protein_id="ENSP00000309572",
    )
    assert position.to_cdna() == [expected]


def test_negative_strand_cdna_protein_start(ensembl100):
    position = RnaPosition(
        _data=ensembl100,
        contig_id="5",
        start=80,
        end=82,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
    )
    expected = CdnaPosition(
        _data=ensembl100,
        contig_id="5",
        start=1,
        end=3,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
        protein_id="ENSP00000309572",
    )
    assert position.to_cdna() == [expected]


def test_positive_strand(ensembl100):
    position = RnaPosition(
        _data=ensembl100,
        contig_id="13",
        start=8799,
        end=8801,
        strand="+",
        gene_id="ENSG00000139618",
        gene_name="BRCA2",
        transcript_id="ENST00000380152",
        transcript_name="BRCA2-201",
    )
    expected = CdnaPosition(
        _data=ensembl100,
        contig_id="13",
        start=8566,
        end=8568,
        strand="+",
        gene_id="ENSG00000139618",
        gene_name="BRCA2",
        transcript_id="ENST00000380152",
        transcript_name="BRCA2-201",
        protein_id="ENSP00000369497",
    )
    assert position.to_cdna() == [expected]


def test_positive_strand_across_exon_boundary(ensembl100):
    position = RnaPosition(
        _data=ensembl100,
        contig_id="4",
        start=125,
        end=127,
        strand="+",
        gene_id="ENSG00000157404",
        gene_name="KIT",
        transcript_id="ENST00000288135",
        transcript_name="KIT-201",
    )
    expected = CdnaPosition(
        _data=ensembl100,
        contig_id="4",
        start=67,
        end=69,
        strand="+",
        gene_id="ENSG00000157404",
        gene_name="KIT",
        transcript_id="ENST00000288135",
        transcript_name="KIT-201",
        protein_id="ENSP00000288135",
    )
    assert position.to_cdna() == [expected]


def test_positive_strand_cdna_protein_end(ensembl100):
    position = RnaPosition(
        _data=ensembl100,
        contig_id="13",
        start=10488,
        end=10490,
        strand="+",
        gene_id="ENSG00000139618",
        gene_name="BRCA2",
        transcript_id="ENST00000380152",
        transcript_name="BRCA2-201",
    )
    expected = CdnaPosition(
        _data=ensembl100,
        contig_id="13",
        start=10255,
        end=10257,
        strand="+",
        gene_id="ENSG00000139618",
        gene_name="BRCA2",
        transcript_id="ENST00000380152",
        transcript_name="BRCA2-201",
        protein_id="ENSP00000369497",
    )
    assert position.to_cdna() == [expected]


def test_positive_strand_cdna_protein_start(ensembl100):
    position = RnaPosition(
        _data=ensembl100,
        contig_id="13",
        start=234,
        end=236,
        strand="+",
        gene_id="ENSG00000139618",
        gene_name="BRCA2",
        transcript_id="ENST00000380152",
        transcript_name="BRCA2-201",
    )
    expected = CdnaPosition(
        _data=ensembl100,
        contig_id="13",
        start=1,
        end=3,
        strand="+",
        gene_id="ENSG00000139618",
        gene_name="BRCA2",
        transcript_id="ENST00000380152",
        transcript_name="BRCA2-201",
        protein_id="ENSP00000369497",
    )
    assert position.to_cdna() == [expected]


def test_positive_strand_overlapping_genes_different_strands(ensembl100):
    position = RnaPosition(
        _data=ensembl100,
        contig_id="6",
        start=580,
        end=582,
        strand="+",
        gene_id="ENSG00000137310",
        gene_name="TCF19",
        transcript_id="ENST00000542218",
        transcript_name="TCF19-204",
    )
    expected = CdnaPosition(
        _data=ensembl100,
        contig_id="6",
        start=580,
        end=582,
        strand="+",
        gene_id="ENSG00000137310",
        gene_name="TCF19",
        transcript_id="ENST00000542218",
        transcript_name="TCF19-204",
        protein_id="ENSP00000439397",
    )
    assert position.to_cdna() == [expected]
