from ensembl_map.core import CdnaPosition, DnaPosition


def test_negative_strand(ensembl100):
    position = DnaPosition(_data=ensembl100, contig_id="5", start=1294864, end=1294866, strand="-")
    expected = CdnaPosition(
        _data=ensembl100,
        contig_id="5",
        start=124,
        end=126,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000334602",
        transcript_name="TERT-202",
        protein_id="ENSP00000334346",
    )
    assert expected in position.to_cdna()


def test_negative_strand_across_exon_boundary(ensembl100):
    position = DnaPosition(_data=ensembl100, contig_id="5", start=1282623, end=1293313, strand="-")
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
    assert expected in position.to_cdna()


def test_negative_strand_cdna_protein_end(ensembl100):
    position = DnaPosition(_data=ensembl100, contig_id="5", start=1253728, end=1253730, strand="-")
    expected = CdnaPosition(
        _data=ensembl100,
        contig_id="5",
        start=3208,
        end=3210,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000334602",
        transcript_name="TERT-202",
        protein_id="ENSP00000334346",
    )
    assert expected in position.to_cdna()


def test_negative_strand_cdna_protein_start(ensembl100):
    position = DnaPosition(_data=ensembl100, contig_id="5", start=1294987, end=1294989, strand="-")
    expected = CdnaPosition(
        _data=ensembl100,
        contig_id="5",
        start=1,
        end=3,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000334602",
        transcript_name="TERT-202",
        protein_id="ENSP00000334346",
    )
    assert expected in position.to_cdna()


def test_positive_strand(ensembl100):
    position = DnaPosition(
        _data=ensembl100, contig_id="13", start=32371034, end=32371036, strand="+"
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
    assert expected in position.to_cdna()


def test_positive_strand_across_exon_boundary(ensembl100):
    position = DnaPosition(
        _data=ensembl100, contig_id="4", start=54658081, end=54695513, strand="+"
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
    assert expected in position.to_cdna()


def test_positive_strand_cdna_protein_end(ensembl100):
    position = DnaPosition(
        _data=ensembl100, contig_id="13", start=32398768, end=32398770, strand="+"
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
    assert expected in position.to_cdna()


def test_positive_strand_cdna_protein_start(ensembl100):
    position = DnaPosition(
        _data=ensembl100, contig_id="13", start=32316461, end=32316463, strand="+"
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
    assert expected in position.to_cdna()


def test_positive_strand_overlapping_genes_different_strands(ensembl100):
    position = DnaPosition(
        _data=ensembl100, contig_id="6", start=31166302, end=31166304, strand="+"
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
