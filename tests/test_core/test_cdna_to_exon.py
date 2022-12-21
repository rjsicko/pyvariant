from ensembl_map.core import CdnaMappablePosition, ExonMappablePosition


def test_negative_strand(ensembl100):
    position = CdnaMappablePosition(
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
    expected = ExonMappablePosition(
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
    assert position.to_exon() == [expected]


def test_negative_strand_cdna_protein_end(ensembl100):
    position = CdnaMappablePosition(
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
    expected = ExonMappablePosition(
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
    assert position.to_exon() == [expected]


def test_negative_strand_cdna_protein_start(ensembl100):
    position = CdnaMappablePosition(
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
    expected = ExonMappablePosition(
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
    assert position.to_exon() == [expected]


def test_positive_strand(ensembl100):
    position = CdnaMappablePosition(
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
    expected = ExonMappablePosition(
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
    assert position.to_exon() == [expected]


def test_positive_strand_cdna_protein_end(ensembl100):
    position = CdnaMappablePosition(
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
    expected = ExonMappablePosition(
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
    assert position.to_exon() == [expected]


def test_positive_strand_cdna_protein_start(ensembl100):
    position = CdnaMappablePosition(
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
    expected = ExonMappablePosition(
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
    assert position.to_exon() == [expected]


def test_positive_strand_overlapping_genes_different_strands(ensembl100):
    position = CdnaMappablePosition(
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
    expected = ExonMappablePosition(
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
    assert position.to_exon() == [expected]


def test_negative_strand_negative_offset(ensembl100):
    position = CdnaMappablePosition(
        _data=ensembl100,
        contig_id="5",
        start=3,
        start_offset=-1,
        end=4,
        end_offset=-1,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
        protein_id="ENSP00000309572",
    )
    expected = ExonMappablePosition(
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
    assert position.to_exon() == [expected]


def test_negative_strand_positive_offset(ensembl100):
    position = CdnaMappablePosition(
        _data=ensembl100,
        contig_id="5",
        start=3,
        start_offset=1,
        end=4,
        end_offset=1,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
        protein_id="ENSP00000309572",
    )
    expected = ExonMappablePosition(
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
    assert position.to_exon() == [expected]


def test_positive_strand_positive_offset(ensembl100):
    position = CdnaMappablePosition(
        _data=ensembl100,
        contig_id="4",
        start=67,
        start_offset=37432,
        end=67,
        end_offset=37432,
        strand="+",
        gene_id="ENSG00000157404",
        gene_name="KIT",
        transcript_id="ENST00000288135",
        transcript_name="KIT-201",
        protein_id="ENSP00000288135",
    )
    expected = ExonMappablePosition(
        _data=ensembl100,
        contig_id="4",
        start=2,
        start_offset=0,
        end=2,
        end_offset=0,
        strand="+",
        gene_id="ENSG00000157404",
        gene_name="KIT",
        transcript_id="ENST00000288135",
        transcript_name="KIT-201",
        exon_id="ENSE00001032350",
    )
    assert position.to_exon() == [expected]


def test_positive_strand_negative_offset(ensembl100):
    position = CdnaMappablePosition(
        _data=ensembl100,
        contig_id="4",
        start=68,
        start_offset=-37432,
        end=68,
        end_offset=-37432,
        strand="+",
        gene_id="ENSG00000157404",
        gene_name="KIT",
        transcript_id="ENST00000288135",
        transcript_name="KIT-201",
        protein_id="ENSP00000288135",
    )
    expected = ExonMappablePosition(
        _data=ensembl100,
        contig_id="4",
        start=1,
        start_offset=0,
        end=1,
        end_offset=0,
        strand="+",
        gene_id="ENSG00000157404",
        gene_name="KIT",
        transcript_id="ENST00000288135",
        transcript_name="KIT-201",
        exon_id="ENSE00000000233",
    )
    assert position.to_exon() == [expected]


def test_positive_strand_positive_offset_into_intron(ensembl100):
    position = CdnaMappablePosition(
        _data=ensembl100,
        contig_id="4",
        start=67,
        start_offset=1,
        end=67,
        end_offset=1,
        strand="+",
        gene_id="ENSG00000157404",
        gene_name="KIT",
        transcript_id="ENST00000288135",
        transcript_name="KIT-201",
        protein_id="ENSP00000288135",
    )
    assert position.to_exon() == []


def test_positive_strand_negative_offset_into_intron(ensembl100):
    position = CdnaMappablePosition(
        _data=ensembl100,
        contig_id="4",
        start=68,
        start_offset=-1,
        end=68,
        end_offset=-1,
        strand="+",
        gene_id="ENSG00000157404",
        gene_name="KIT",
        transcript_id="ENST00000288135",
        transcript_name="KIT-201",
        protein_id="ENSP00000288135",
    )
    assert position.to_exon() == []
