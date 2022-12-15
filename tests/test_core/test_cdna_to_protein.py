from ensembl_map.core import CdnaPosition, ProteinPosition


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
    expected = ProteinPosition(
        _data=ensembl100,
        contig_id="5",
        start=42,
        end=42,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
        protein_id="ENSP00000309572",
    )
    assert position.to_protein() == [expected]


def test_negative_strand_across_exon_boundary(ensembl100):
    position = CdnaPosition(
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
    expected = ProteinPosition(
        _data=ensembl100,
        contig_id="5",
        start=525,
        end=525,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
        protein_id="ENSP00000309572",
    )
    assert position.to_protein() == [expected]


def test_negative_strand_cdna_protein_end(ensembl100):
    position = CdnaPosition(
        _data=ensembl100,
        contig_id="5",
        start=3394,
        end=3396,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
        protein_id="ENSP00000309572",
    )
    expected = ProteinPosition(
        _data=ensembl100,
        contig_id="5",
        start=1132,
        end=1132,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
        protein_id="ENSP00000309572",
    )
    assert position.to_protein() == [expected]


def test_negative_strand_stop_codon(ensembl100):
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
    assert position.to_protein() == []


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
    expected = ProteinPosition(
        _data=ensembl100,
        contig_id="5",
        start=1,
        end=1,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
        protein_id="ENSP00000309572",
    )
    assert position.to_protein() == [expected]


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
    expected = ProteinPosition(
        _data=ensembl100,
        contig_id="13",
        start=2856,
        end=2856,
        strand="+",
        gene_id="ENSG00000139618",
        gene_name="BRCA2",
        transcript_id="ENST00000380152",
        transcript_name="BRCA2-201",
        protein_id="ENSP00000369497",
    )
    assert position.to_protein() == [expected]


def test_positive_strand_across_exon_boundary(ensembl100):
    position = CdnaPosition(
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
    expected = ProteinPosition(
        _data=ensembl100,
        contig_id="4",
        start=23,
        end=23,
        strand="+",
        gene_id="ENSG00000157404",
        gene_name="KIT",
        transcript_id="ENST00000288135",
        transcript_name="KIT-201",
        protein_id="ENSP00000288135",
    )
    assert position.to_protein() == [expected]


def test_positive_strand_cdna_protein_end(ensembl100):
    position = CdnaPosition(
        _data=ensembl100,
        contig_id="13",
        start=10252,
        end=10254,
        strand="+",
        gene_id="ENSG00000139618",
        gene_name="BRCA2",
        transcript_id="ENST00000380152",
        transcript_name="BRCA2-201",
        protein_id="ENSP00000369497",
    )
    expected = ProteinPosition(
        _data=ensembl100,
        contig_id="13",
        start=3418,
        end=3418,
        strand="+",
        gene_id="ENSG00000139618",
        gene_name="BRCA2",
        transcript_id="ENST00000380152",
        transcript_name="BRCA2-201",
        protein_id="ENSP00000369497",
    )
    assert position.to_protein() == [expected]


def test_positive_strand_stop_codon(ensembl100):
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
    assert position.to_protein() == []


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
    expected = ProteinPosition(
        _data=ensembl100,
        contig_id="13",
        start=1,
        end=1,
        strand="+",
        gene_id="ENSG00000139618",
        gene_name="BRCA2",
        transcript_id="ENST00000380152",
        transcript_name="BRCA2-201",
        protein_id="ENSP00000369497",
    )
    assert position.to_protein() == [expected]


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
    expected = ProteinPosition(
        _data=ensembl100,
        contig_id="6",
        start=194,
        end=194,
        strand="+",
        gene_id="ENSG00000137310",
        gene_name="TCF19",
        transcript_id="ENST00000542218",
        transcript_name="TCF19-204",
        protein_id="ENSP00000439397",
    )
    assert position.to_protein() == [expected]
