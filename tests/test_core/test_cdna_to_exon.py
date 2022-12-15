from ensembl_map.core import CdnaPosition, ExonPosition


def test_negative_strand(ensembl100):
    position = CdnaPosition(
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
    position = CdnaPosition(
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
    position = CdnaPosition(
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
    position = CdnaPosition(
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
    position = CdnaPosition(
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
    position = CdnaPosition(
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
    position = CdnaPosition(
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
