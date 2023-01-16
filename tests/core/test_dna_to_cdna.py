from variant_map.core import CdnaMappablePosition, DnaMappablePosition


def test_negative_strand(ensembl100):
    position = DnaMappablePosition(
        _data=ensembl100,
        contig_id="5",
        start=1294864,
        start_offset=0,
        end=1294866,
        end_offset=0,
        strand="-",
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
        transcript_id="ENST00000334602",
        transcript_name="TERT-202",
        protein_id="ENSP00000334346",
    )
    assert expected in position.to_cdna()


def test_negative_strand_across_exon_boundary(ensembl100):
    position = DnaMappablePosition(
        _data=ensembl100,
        contig_id="5",
        start=1282623,
        start_offset=0,
        end=1293313,
        end_offset=0,
        strand="-",
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
    assert expected in position.to_cdna()


def test_negative_strand_cdna_protein_end(ensembl100):
    position = DnaMappablePosition(
        _data=ensembl100,
        contig_id="5",
        start=1253728,
        start_offset=0,
        end=1253730,
        end_offset=0,
        strand="-",
    )
    expected = CdnaMappablePosition(
        _data=ensembl100,
        contig_id="5",
        start=3208,
        start_offset=0,
        end=3210,
        end_offset=0,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000334602",
        transcript_name="TERT-202",
        protein_id="ENSP00000334346",
    )
    assert expected in position.to_cdna()


def test_negative_strand_cdna_protein_start(ensembl100):
    position = DnaMappablePosition(
        _data=ensembl100,
        contig_id="5",
        start=1294987,
        start_offset=0,
        end=1294989,
        end_offset=0,
        strand="-",
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
        transcript_id="ENST00000334602",
        transcript_name="TERT-202",
        protein_id="ENSP00000334346",
    )
    assert expected in position.to_cdna()


def test_positive_strand(ensembl100):
    position = DnaMappablePosition(
        _data=ensembl100,
        contig_id="13",
        start=32371034,
        start_offset=0,
        end=32371036,
        end_offset=0,
        strand="+",
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
    assert expected in position.to_cdna()


def test_positive_strand_across_exon_boundary(ensembl100):
    position = DnaMappablePosition(
        _data=ensembl100,
        contig_id="4",
        start=54658081,
        start_offset=0,
        end=54695513,
        end_offset=0,
        strand="+",
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
    assert expected in position.to_cdna()


def test_positive_strand_cdna_protein_end(ensembl100):
    position = DnaMappablePosition(
        _data=ensembl100,
        contig_id="13",
        start=32398768,
        start_offset=0,
        end=32398770,
        end_offset=0,
        strand="+",
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
    assert expected in position.to_cdna()


def test_positive_strand_cdna_protein_start(ensembl100):
    position = DnaMappablePosition(
        _data=ensembl100,
        contig_id="13",
        start=32316461,
        start_offset=0,
        end=32316463,
        end_offset=0,
        strand="+",
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
    assert expected in position.to_cdna()


def test_positive_strand_overlapping_genes_different_strands(ensembl100):
    position = DnaMappablePosition(
        _data=ensembl100,
        contig_id="6",
        start=31166302,
        start_offset=0,
        end=31166304,
        end_offset=0,
        strand="+",
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


def test_offset_into_intron(ensembl100):
    position = DnaMappablePosition(
        _data=ensembl100,
        contig_id="5",
        start=1294666,
        start_offset=-1,
        end=1294666,
        end_offset=-1,
        strand="-",
    )
    assert position.to_cdna() == []


def test_offset_out_of_intron(ensembl100):
    position = DnaMappablePosition(
        _data=ensembl100,
        contig_id="5",
        start=1294667,
        start_offset=1,
        end=1294667,
        end_offset=1,
        strand="-",
    )
    expected = CdnaMappablePosition(
        _data=ensembl100,
        contig_id="5",
        start=220,
        start_offset=0,
        end=220,
        end_offset=0,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
        protein_id="ENSP00000309572",
    )
    assert expected in position.to_cdna()
