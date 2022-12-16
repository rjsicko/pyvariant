from ensembl_map.core import CdnaPosition, DnaPosition


def test_negative_strand(ensembl100):
    position = CdnaPosition(
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
    expected = DnaPosition(
        _data=ensembl100,
        contig_id="5",
        start=1294864,
        start_offset=0,
        end=1294866,
        end_offset=0,
        strand="-",
    )
    assert position.to_dna() == [expected]


def test_negative_strand_across_exon_boundary(ensembl100):
    position = CdnaPosition(
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
    expected = DnaPosition(
        _data=ensembl100,
        contig_id="5",
        start=1282623,
        start_offset=0,
        end=1293313,
        end_offset=0,
        strand="-",
    )
    assert position.to_dna() == [expected]


def test_negative_strand_cdna_protein_end(ensembl100):
    position = CdnaPosition(
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
    expected = DnaPosition(
        _data=ensembl100,
        contig_id="5",
        start=1253728,
        start_offset=0,
        end=1253730,
        end_offset=0,
        strand="-",
    )
    assert position.to_dna() == [expected]


def test_negative_strand_cdna_protein_start(ensembl100):
    position = CdnaPosition(
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
    expected = DnaPosition(
        _data=ensembl100,
        contig_id="5",
        start=1294987,
        start_offset=0,
        end=1294989,
        end_offset=0,
        strand="-",
    )
    assert position.to_dna() == [expected]


def test_positive_strand(ensembl100):
    position = CdnaPosition(
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
    expected = DnaPosition(
        _data=ensembl100,
        contig_id="13",
        start=32371034,
        start_offset=0,
        end=32371036,
        end_offset=0,
        strand="+",
    )
    assert position.to_dna() == [expected]


def test_positive_strand_across_exon_boundary(ensembl100):
    position = CdnaPosition(
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
    expected = DnaPosition(
        _data=ensembl100,
        contig_id="4",
        start=54658081,
        start_offset=0,
        end=54695513,
        end_offset=0,
        strand="+",
    )
    assert position.to_dna() == [expected]


def test_positive_strand_cdna_protein_end(ensembl100):
    position = CdnaPosition(
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
    expected = DnaPosition(
        _data=ensembl100,
        contig_id="13",
        start=32398768,
        start_offset=0,
        end=32398770,
        end_offset=0,
        strand="+",
    )
    assert position.to_dna() == [expected]


def test_positive_strand_cdna_protein_start(ensembl100):
    position = CdnaPosition(
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
    expected = DnaPosition(
        _data=ensembl100,
        contig_id="13",
        start=32316461,
        start_offset=0,
        end=32316463,
        end_offset=0,
        strand="+",
    )
    assert position.to_dna() == [expected]


def test_positive_strand_overlapping_genes_different_strands(ensembl100):
    position = CdnaPosition(
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
    expected = DnaPosition(
        _data=ensembl100,
        contig_id="6",
        start=31166302,
        start_offset=0,
        end=31166304,
        end_offset=0,
        strand="+",
    )
    assert position.to_dna() == [expected]


def test_negative_strand_negative_offset(ensembl100):
    position = CdnaPosition(
        _data=ensembl100,
        contig_id="5",
        start=220,
        start_offset=-2,
        end=221,
        end_offset=-4,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
        protein_id="ENSP00000309572",
    )
    expected = DnaPosition(
        _data=ensembl100,
        contig_id="5",
        start=1294668,
        start_offset=0,
        end=1294669,
        end_offset=0,
        strand="-",
    )
    assert position.to_dna() == [expected]


def test_negative_strand_positive_offset(ensembl100):
    position = CdnaPosition(
        _data=ensembl100,
        contig_id="5",
        start=218,
        start_offset=4,
        end=219,
        end_offset=6,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
        protein_id="ENSP00000309572",
    )
    expected = DnaPosition(
        _data=ensembl100,
        contig_id="5",
        start=1294765,
        start_offset=0,
        end=1294768,
        end_offset=0,
        strand="-",
    )
    assert position.to_dna() == [expected]


def test_positive_strand_negative_offset(ensembl100):
    position = CdnaPosition(
        _data=ensembl100,
        contig_id="4",
        start=68,
        start_offset=-4,
        end=69,
        end_offset=-2,
        strand="+",
        gene_id="ENSG00000157404",
        gene_name="KIT",
        transcript_id="ENST00000288135",
        transcript_name="KIT-201",
        protein_id="ENSP00000288135",
    )
    expected = DnaPosition(
        _data=ensembl100,
        contig_id="4",
        start=54695508,
        start_offset=0,
        end=54695511,
        end_offset=0,
        strand="+",
    )
    assert position.to_dna() == [expected]


def test_positive_strand_positive_offset(ensembl100):
    position = CdnaPosition(
        _data=ensembl100,
        contig_id="4",
        start=66,
        start_offset=2,
        end=67,
        end_offset=4,
        strand="+",
        gene_id="ENSG00000157404",
        gene_name="KIT",
        transcript_id="ENST00000288135",
        transcript_name="KIT-201",
        protein_id="ENSP00000288135",
    )
    expected = DnaPosition(
        _data=ensembl100,
        contig_id="4",
        start=54658082,
        start_offset=0,
        end=54658085,
        end_offset=0,
        strand="+",
    )
    assert position.to_dna() == [expected]
