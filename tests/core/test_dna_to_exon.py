from variant_map.core import DnaMappablePosition, ExonMappablePosition


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
    assert expected in position.to_exon()


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
    assert expected in position.to_exon()


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
    assert expected in position.to_exon()


def test_negative_strand_overlapping_genes_different_strands(ensembl100):
    position = DnaMappablePosition(
        _data=ensembl100,
        contig_id="6",
        start=31166302,
        start_offset=0,
        end=31166304,
        end_offset=0,
        strand="-",
    )
    expected = [
        ExonMappablePosition(
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
        ),
        ExonMappablePosition(
            _data=ensembl100,
            contig_id="6",
            start=1,
            start_offset=0,
            end=1,
            end_offset=0,
            strand="-",
            gene_id="ENSG00000204531",
            gene_name="POU5F1",
            transcript_id="ENST00000513407",
            transcript_name="POU5F1-206",
            exon_id="ENSE00002033137",
        ),
    ]
    assert position.to_exon() == expected


def test_negative_strand_transcript_end(ensembl100):
    position = DnaMappablePosition(
        _data=ensembl100,
        contig_id="5",
        start=1253167,
        start_offset=0,
        end=1253169,
        end_offset=0,
        strand="-",
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
    assert expected in position.to_exon()


def test_negative_strand_transcript_start(ensembl100):
    position = DnaMappablePosition(
        _data=ensembl100,
        contig_id="5",
        start=1295066,
        start_offset=0,
        end=1295068,
        end_offset=0,
        strand="-",
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
    assert expected in position.to_exon()


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
    assert expected in position.to_exon()


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
    assert expected in position.to_exon()


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
    assert expected in position.to_exon()


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


def test_positive_strand_transcript_end(ensembl100):
    position = DnaMappablePosition(
        _data=ensembl100,
        contig_id="13",
        start=32400264,
        start_offset=0,
        end=32400266,
        end_offset=0,
        strand="+",
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
    assert expected in position.to_exon()


def test_positive_strand_transcript_start(ensembl100):
    position = DnaMappablePosition(
        _data=ensembl100,
        contig_id="13",
        start=32315474,
        start_offset=0,
        end=32315476,
        end_offset=0,
        strand="+",
    )
    expected = ExonMappablePosition(
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
    assert expected in position.to_exon()


def test_negative_strand_negative_offset(ensembl100):
    position = DnaMappablePosition(
        _data=ensembl100,
        contig_id="5",
        start=1294770,
        start_offset=-1,
        end=1294770,
        end_offset=-1,
        strand="-",
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
    assert expected in position.to_exon()


def test_negative_strand_positive_offset(ensembl100):
    position = DnaMappablePosition(
        _data=ensembl100,
        contig_id="5",
        start=1294667,
        start_offset=1,
        end=1294667,
        end_offset=1,
        strand="-",
    )
    expected = ExonMappablePosition(
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
    assert expected in position.to_exon()


def test_positive_strand_negative_offset(ensembl100):
    position = DnaMappablePosition(
        _data=ensembl100,
        contig_id="4",
        start=54658082,
        start_offset=-1,
        end=54658082,
        end_offset=-1,
        strand="+",
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
    assert expected in position.to_exon()


def test_positive_strand_positive_offset(ensembl100):
    position = DnaMappablePosition(
        _data=ensembl100,
        contig_id="4",
        start=54695511,
        start_offset=1,
        end=54695511,
        end_offset=1,
        strand="+",
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
    assert expected in position.to_exon()


def test_positive_strand_positive_offset_into_intron(ensembl100):
    position = DnaMappablePosition(
        _data=ensembl100,
        contig_id="4",
        start=54658081,
        start_offset=1,
        end=54658081,
        end_offset=1,
        strand="+",
    )
    assert position.to_exon() == []


def test_positive_strand_negative_offset_into_intron(ensembl100):
    position = DnaMappablePosition(
        _data=ensembl100,
        contig_id="4",
        start=54695512,
        start_offset=-1,
        end=54695512,
        end_offset=-1,
        strand="+",
    )
    assert position.to_exon() == []
