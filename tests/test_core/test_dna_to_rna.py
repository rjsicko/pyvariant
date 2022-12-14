from ensembl_map.core import DnaPosition, RnaPosition


def test_negative_strand(ensembl100):
    position = DnaPosition(_data=ensembl100, contig_id="5", start=1294864, end=1294866, strand="-")
    expected = RnaPosition(
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
    assert expected in position.to_rna()


def test_negative_strand_across_exon_boundary(ensembl100):
    position = DnaPosition(_data=ensembl100, contig_id="5", start=1282623, end=1293313, strand="-")
    expected = RnaPosition(
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
    assert expected in position.to_rna()


def test_negative_strand_cdna_protein_end(ensembl100):
    position = DnaPosition(_data=ensembl100, contig_id="5", start=1253728, end=1253730, strand="-")
    expected = RnaPosition(
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
    assert expected in position.to_rna()


def test_negative_strand_cdna_protein_start(ensembl100):
    position = DnaPosition(_data=ensembl100, contig_id="5", start=1294987, end=1294989, strand="-")
    expected = RnaPosition(
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
    assert expected in position.to_rna()


def test_negative_strand_overlapping_genes_different_strands(ensembl100):
    position = DnaPosition(
        _data=ensembl100, contig_id="6", start=31166302, end=31166304, strand="-"
    )
    expected = [
        RnaPosition(
            _data=ensembl100,
            contig_id="6",
            start=535,
            end=537,
            strand="-",
            gene_id="ENSG00000204531",
            gene_name="POU5F1",
            transcript_id="ENST00000471529",
            transcript_name="POU5F1-204",
        ),
        RnaPosition(
            _data=ensembl100,
            contig_id="6",
            start=867,
            end=869,
            strand="-",
            gene_id="ENSG00000204531",
            gene_name="POU5F1",
            transcript_id="ENST00000513407",
            transcript_name="POU5F1-206",
        ),
    ]
    assert position.to_rna() == expected


def test_negative_strand_transcript_end(ensembl100):
    position = DnaPosition(_data=ensembl100, contig_id="5", start=1253167, end=1253169, strand="-")
    expected = RnaPosition(
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
    assert expected in position.to_rna()


def test_negative_strand_transcript_start(ensembl100):
    position = DnaPosition(_data=ensembl100, contig_id="5", start=1295066, end=1295068, strand="-")
    expected = RnaPosition(
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
    assert expected in position.to_rna()


def test_positive_strand(ensembl100):
    position = DnaPosition(
        _data=ensembl100, contig_id="13", start=32371034, end=32371036, strand="+"
    )
    expected = RnaPosition(
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
    assert expected in position.to_rna()


def test_positive_strand_across_exon_boundary(ensembl100):
    position = DnaPosition(
        _data=ensembl100, contig_id="4", start=54658081, end=54695513, strand="+"
    )
    expected = RnaPosition(
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
    assert expected in position.to_rna()


def test_positive_strand_cdna_protein_end(ensembl100):
    position = DnaPosition(
        _data=ensembl100, contig_id="13", start=32398768, end=32398770, strand="+"
    )
    expected = RnaPosition(
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
    assert expected in position.to_rna()


def test_positive_strand_cdna_protein_start(ensembl100):
    position = DnaPosition(
        _data=ensembl100, contig_id="13", start=32316461, end=32316463, strand="+"
    )
    expected = RnaPosition(
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
    assert expected in position.to_rna()


def test_positive_strand_overlapping_genes_different_strands(ensembl100):
    position = DnaPosition(
        _data=ensembl100, contig_id="6", start=31166302, end=31166304, strand="+"
    )
    expected = RnaPosition(
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
    assert position.to_rna() == [expected]


def test_positive_strand_transcript_end(ensembl100):
    position = DnaPosition(
        _data=ensembl100, contig_id="13", start=32400264, end=32400266, strand="+"
    )
    expected = RnaPosition(
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
    assert expected in position.to_rna()


def test_positive_strand_transcript_start(ensembl100):
    position = DnaPosition(
        _data=ensembl100, contig_id="13", start=32315474, end=32315476, strand="+"
    )
    expected = RnaPosition(
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
    assert expected in position.to_rna()
