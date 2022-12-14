from ensembl_map.core import ExonPosition, RnaPosition


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
    expected = ExonPosition(
        _data=ensembl100,
        contig_id="5",
        start=1,
        end=1,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
        exon_id="ENSE00003896691",
    )
    assert position.to_exon() == [expected]


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
    expected = ExonPosition(
        _data=ensembl100,
        contig_id="5",
        start=16,
        end=16,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
        exon_id="ENSE00001863787",
    )
    assert position.to_exon() == [expected]


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
    expected = ExonPosition(
        _data=ensembl100,
        contig_id="5",
        start=1,
        end=1,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
        exon_id="ENSE00003896691",
    )
    assert position.to_exon() == [expected]


def test_negative_strand_overlapping_genes_different_strands(ensembl100):
    position = RnaPosition(
        _data=ensembl100,
        contig_id="6",
        start=535,
        end=537,
        strand="-",
        gene_id="ENSG00000204531",
        gene_name="POU5F1",
        transcript_id="ENST00000471529",
        transcript_name="POU5F1-204",
    )
    expected = ExonPosition(
        _data=ensembl100,
        contig_id="6",
        start=1,
        end=1,
        strand="-",
        gene_id="ENSG00000204531",
        gene_name="POU5F1",
        transcript_id="ENST00000471529",
        transcript_name="POU5F1-204",
        exon_id="ENSE00002568331",
    )
    assert position.to_exon() == [expected]


def test_negative_strand_transcript_end(ensembl100):
    position = RnaPosition(
        _data=ensembl100,
        contig_id="5",
        start=4037,
        end=4039,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
    )
    expected = ExonPosition(
        _data=ensembl100,
        contig_id="5",
        start=16,
        end=16,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
        exon_id="ENSE00001863787",
    )
    assert position.to_exon() == [expected]


def test_negative_strand_transcript_start(ensembl100):
    position = RnaPosition(
        _data=ensembl100,
        contig_id="5",
        start=1,
        end=3,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
    )
    expected = ExonPosition(
        _data=ensembl100,
        contig_id="5",
        start=1,
        end=1,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
        exon_id="ENSE00003896691",
    )
    assert position.to_exon() == [expected]


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
    expected = ExonPosition(
        _data=ensembl100,
        contig_id="13",
        start=20,
        end=20,
        strand="+",
        gene_id="ENSG00000139618",
        gene_name="BRCA2",
        transcript_id="ENST00000380152",
        transcript_name="BRCA2-201",
        exon_id="ENSE00000939180",
    )
    assert position.to_exon() == [expected]


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
    expected = ExonPosition(
        _data=ensembl100,
        contig_id="13",
        start=27,
        end=27,
        strand="+",
        gene_id="ENSG00000139618",
        gene_name="BRCA2",
        transcript_id="ENST00000380152",
        transcript_name="BRCA2-201",
        exon_id="ENSE00003717596",
    )
    assert position.to_exon() == [expected]


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
    expected = ExonPosition(
        _data=ensembl100,
        contig_id="13",
        start=2,
        end=2,
        strand="+",
        gene_id="ENSG00000139618",
        gene_name="BRCA2",
        transcript_id="ENST00000380152",
        transcript_name="BRCA2-201",
        exon_id="ENSE00001484009",
    )
    assert position.to_exon() == [expected]


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
    expected = ExonPosition(
        _data=ensembl100,
        contig_id="6",
        start=2,
        end=2,
        strand="+",
        gene_id="ENSG00000137310",
        gene_name="TCF19",
        transcript_id="ENST00000542218",
        transcript_name="TCF19-204",
        exon_id="ENSE00002283659",
    )
    assert position.to_exon() == [expected]


def test_positive_strand_transcript_end(ensembl100):
    position = RnaPosition(
        _data=ensembl100,
        contig_id="13",
        start=11984,
        end=11986,
        strand="+",
        gene_id="ENSG00000139618",
        gene_name="BRCA2",
        transcript_id="ENST00000380152",
        transcript_name="BRCA2-201",
    )
    expected = ExonPosition(
        _data=ensembl100,
        contig_id="13",
        start=27,
        end=27,
        strand="+",
        gene_id="ENSG00000139618",
        gene_name="BRCA2",
        transcript_id="ENST00000380152",
        transcript_name="BRCA2-201",
        exon_id="ENSE00003717596",
    )
    assert position.to_exon() == [expected]


def test_positive_strand_transcript_start(ensembl100):
    position = RnaPosition(
        _data=ensembl100,
        contig_id="13",
        start=1,
        end=3,
        strand="+",
        gene_id="ENSG00000139618",
        gene_name="BRCA2",
        transcript_id="ENST00000380152",
        transcript_name="BRCA2-201",
    )
    expected = ExonPosition(
        _data=ensembl100,
        contig_id="13",
        start=1,
        end=1,
        strand="+",
        gene_id="ENSG00000139618",
        gene_name="BRCA2",
        transcript_id="ENST00000380152",
        transcript_name="BRCA2-201",
        exon_id="ENSE00001184784",
    )
    assert position.to_exon() == [expected]
