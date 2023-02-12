import pytest

from variant_map.positions import CdnaPosition


def test_is_cdna(ensembl100):
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
    assert position.is_cdna is True


def test_is_dna(ensembl100):
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
    assert position.is_dna is False


def test_is_exon(ensembl100):
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
    assert position.is_exon is False


def test_is_protein(ensembl100):
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
    assert position.is_protein is False


def test_is_rna(ensembl100):
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
    assert position.is_rna is False


def test_is_on_negative_strand(ensembl100):
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
    assert position.on_negative_strand is True


def test_is_on_positive_strand(ensembl100):
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
    assert position.on_positive_strand is False


def test_sequence(ensembl100):
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
    assert position.sequence() == "GGG"


def test_sequence_offset(ensembl100):
    position = CdnaPosition(
        _data=ensembl100,
        contig_id="5",
        start=1573,
        start_offset=1,
        end=1575,
        end_offset=1,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
        protein_id="ENSP00000309572",
    )
    with pytest.raises(ValueError):
        position.sequence()


def test_str_mutli_position(ensembl100):
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
    assert str(position) == "ENST00000310581:c.1573_1575"


def test_str_one_position(ensembl100):
    position = CdnaPosition(
        _data=ensembl100,
        contig_id="5",
        start=1573,
        start_offset=0,
        end=1573,
        end_offset=0,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
        protein_id="ENSP00000309572",
    )
    assert str(position) == "ENST00000310581:c.1573"
