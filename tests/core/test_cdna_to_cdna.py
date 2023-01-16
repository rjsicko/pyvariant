from variant_map.core import CdnaMappablePosition


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
    assert position.to_cdna() == [position]


def test_negative_strand_across_exon_boundary(ensembl100):
    position = CdnaMappablePosition(
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
    assert position.to_cdna() == [position]


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
    assert position.to_cdna() == [position]


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
    assert position.to_cdna() == [position]


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
    assert position.to_cdna() == [position]


def test_positive_strand_across_exon_boundary(ensembl100):
    position = CdnaMappablePosition(
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
    assert position.to_cdna() == [position]


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
    assert position.to_cdna() == [position]


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
    assert position.to_cdna() == [position]


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
    assert position.to_cdna() == [position]


def test_offset(ensembl100):
    position = CdnaMappablePosition(
        _data=ensembl100,
        contig_id="5",
        start=219,
        start_offset=1,
        end=219,
        end_offset=1,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
        protein_id="ENSP00000309572",
    )
    assert position.to_cdna() == [position]


def test_offset_normalized(ensembl100):
    position = CdnaMappablePosition(
        _data=ensembl100,
        contig_id="5",
        start=100,
        start_offset=-1,
        end=100,
        end_offset=-1,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
        protein_id="ENSP00000309572",
    )
    expected = CdnaMappablePosition(
        _data=ensembl100,
        contig_id="5",
        start=99,
        start_offset=0,
        end=99,
        end_offset=0,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
        protein_id="ENSP00000309572",
    )
    assert position.to_cdna() == [expected]
