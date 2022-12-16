import pytest

from ensembl_map.core import DnaPosition, ProteinPosition


def test_negative_strand(ensembl100):
    position = ProteinPosition(
        _data=ensembl100,
        contig_id="5",
        start=42,
        start_offset=0,
        end=42,
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
    position = ProteinPosition(
        _data=ensembl100,
        contig_id="5",
        start=525,
        start_offset=0,
        end=525,
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
    position = ProteinPosition(
        _data=ensembl100,
        contig_id="5",
        start=1132,
        start_offset=0,
        end=1132,
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
        start=1253731,
        start_offset=0,
        end=1253733,
        end_offset=0,
        strand="-",
    )
    assert position.to_dna() == [expected]


def test_negative_strand_cdna_protein_start(ensembl100):
    position = ProteinPosition(
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
    position = ProteinPosition(
        _data=ensembl100,
        contig_id="13",
        start=2856,
        start_offset=0,
        end=2856,
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
    position = ProteinPosition(
        _data=ensembl100,
        contig_id="4",
        start=23,
        start_offset=0,
        end=23,
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
    position = ProteinPosition(
        _data=ensembl100,
        contig_id="13",
        start=3418,
        start_offset=0,
        end=3418,
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
        start=32398765,
        start_offset=0,
        end=32398767,
        end_offset=0,
        strand="+",
    )
    assert position.to_dna() == [expected]


def test_positive_strand_cdna_protein_start(ensembl100):
    position = ProteinPosition(
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
    position = ProteinPosition(
        _data=ensembl100,
        contig_id="6",
        start=194,
        start_offset=0,
        end=194,
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


def test_offset_error(ensembl100):
    position = ProteinPosition(
        _data=ensembl100,
        contig_id="4",
        start=21,
        start_offset=1,
        end=21,
        end_offset=1,
        strand="+",
        gene_id="ENSG00000157404",
        gene_name="KIT",
        transcript_id="ENST00000288135",
        transcript_name="KIT-201",
        protein_id="ENSP00000288135",
    )
    with pytest.raises(AssertionError):
        position.to_dna()
