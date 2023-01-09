import pytest

from ensembl_map.core import ExonMappablePosition, ProteinMappablePosition


def test_negative_strand(ensembl100):
    position = ExonMappablePosition(
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
    expected = ProteinMappablePosition(
        _data=ensembl100,
        contig_id="5",
        start=74,
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
    assert position.to_protein() == [expected]


def test_positive_strand(ensembl100):
    position = ExonMappablePosition(
        _data=ensembl100,
        contig_id="13",
        start=3,
        start_offset=0,
        end=3,
        end_offset=0,
        strand="+",
        gene_id="ENSG00000139618",
        gene_name="BRCA2",
        transcript_id="ENST00000380152",
        transcript_name="BRCA2-201",
        exon_id="ENSE00003666217",
    )
    expected = ProteinMappablePosition(
        _data=ensembl100,
        contig_id="13",
        start=23,
        start_offset=0,
        end=106,
        end_offset=0,
        strand="+",
        gene_id="ENSG00000139618",
        gene_name="BRCA2",
        transcript_id="ENST00000380152",
        transcript_name="BRCA2-201",
        protein_id="ENSP00000369497",
    )
    assert position.to_protein() == [expected]


def test_offset_error(ensembl100):
    position = ExonMappablePosition(
        _data=ensembl100,
        contig_id="5",
        start=1,
        start_offset=1,
        end=1,
        end_offset=1,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
        exon_id="ENSE00001197112",
    )
    with pytest.raises(AssertionError):
        position.to_protein()
