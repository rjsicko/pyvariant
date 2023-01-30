import pytest

from variant_map.core import (
    CdnaMappablePosition,
    DnaMappablePosition,
    ExonMappablePosition,
    ProteinMappablePosition,
    RnaMappablePosition,
)


def test_is_cdna(ensembl100):
    position = DnaMappablePosition(
        _data=ensembl100,
        contig_id="5",
        start=1282623,
        start_offset=0,
        end=1293313,
        end_offset=0,
        strand="-",
    )
    assert position.is_cdna is False


def test_is_dna(ensembl100):
    position = DnaMappablePosition(
        _data=ensembl100,
        contig_id="5",
        start=1282623,
        start_offset=0,
        end=1293313,
        end_offset=0,
        strand="-",
    )
    assert position.is_dna is True


def test_is_exon(ensembl100):
    position = DnaMappablePosition(
        _data=ensembl100,
        contig_id="5",
        start=1282623,
        start_offset=0,
        end=1293313,
        end_offset=0,
        strand="-",
    )
    assert position.is_exon is False


def test_is_protein(ensembl100):
    position = DnaMappablePosition(
        _data=ensembl100,
        contig_id="5",
        start=1282623,
        start_offset=0,
        end=1293313,
        end_offset=0,
        strand="-",
    )
    assert position.is_protein is False


def test_is_rna(ensembl100):
    position = DnaMappablePosition(
        _data=ensembl100,
        contig_id="5",
        start=1282623,
        start_offset=0,
        end=1293313,
        end_offset=0,
        strand="-",
    )
    assert position.is_rna is False


def test_is_on_negative_strand(ensembl100):
    position = DnaMappablePosition(
        _data=ensembl100,
        contig_id="5",
        start=1282623,
        start_offset=0,
        end=1293313,
        end_offset=0,
        strand="-",
    )
    assert position.on_negative_strand is True


def test_is_on_positive_strand(ensembl100):
    position = DnaMappablePosition(
        _data=ensembl100,
        contig_id="5",
        start=1282623,
        start_offset=0,
        end=1293313,
        end_offset=0,
        strand="-",
    )
    assert position.on_positive_strand is False


def test_sequence(ensembl100):
    position = DnaMappablePosition(
        _data=ensembl100,
        contig_id="5",
        start=1282623,
        start_offset=0,
        end=1293313,
        end_offset=0,
        strand="-",
    )
    sequence = position.sequence()
    assert (sequence[:3], sequence[-3:]) == ("GGT", "GGG")


def test_sequence_offset(ensembl100):
    position = DnaMappablePosition(
        _data=ensembl100,
        contig_id="5",
        start=1282623,
        start_offset=1,
        end=1293313,
        end_offset=1,
        strand="-",
    )
    with pytest.raises(ValueError):
        position.sequence()


def test_str_mutli_position(ensembl100):
    position = DnaMappablePosition(
        _data=ensembl100,
        contig_id="5",
        start=1282623,
        start_offset=0,
        end=1293313,
        end_offset=0,
        strand="-",
    )
    assert str(position) == "5:g.1282623_1293313"


def test_str_one_position(ensembl100):
    position = DnaMappablePosition(
        _data=ensembl100,
        contig_id="5",
        start=1282623,
        start_offset=0,
        end=1282623,
        end_offset=0,
        strand="-",
    )
    assert str(position) == "5:g.1282623"


def test_to_cdna_dna_end_minus_strand(ensembl100):
    from_position = DnaMappablePosition(
        _data=ensembl100,
        contig_id="5",
        start=1253147,
        start_offset=0,
        end=1253149,
        end_offset=0,
        strand="-",
    )
    to_positions = []
    assert from_position.to_cdna() == to_positions


def test_to_cdna_dna_end_plus_strand(ensembl100):
    from_position = DnaMappablePosition(
        _data=ensembl100,
        contig_id="13",
        start=32400266,
        start_offset=0,
        end=32400268,
        end_offset=0,
        strand="+",
    )
    to_positions = []
    assert from_position.to_cdna() == to_positions


def test_to_cdna_dna_exon_boundary_minus_strand(ensembl100):
    from_position = DnaMappablePosition(
        _data=ensembl100,
        contig_id="5",
        start=1282623,
        start_offset=0,
        end=1293313,
        end_offset=0,
        strand="-",
    )
    to_positions = [
        CdnaMappablePosition(
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
        ),
        CdnaMappablePosition(
            _data=ensembl100,
            contig_id="5",
            start=1573,
            start_offset=0,
            end=1575,
            end_offset=0,
            strand="-",
            gene_id="ENSG00000164362",
            gene_name="TERT",
            transcript_id="ENST00000334602",
            transcript_name="TERT-202",
            protein_id="ENSP00000334346",
        ),
        CdnaMappablePosition(
            _data=ensembl100,
            contig_id="5",
            start=1573,
            start_offset=0,
            end=1575,
            end_offset=0,
            strand="-",
            gene_id="ENSG00000164362",
            gene_name="TERT",
            transcript_id="ENST00000460137",
            transcript_name="TERT-203",
            protein_id="ENSP00000425003",
        ),
    ]
    assert from_position.to_cdna() == to_positions


def test_to_cdna_dna_exon_boundary_plus_strand(ensembl100):
    from_position = DnaMappablePosition(
        _data=ensembl100,
        contig_id="13",
        start=32316527,
        start_offset=0,
        end=32319078,
        end_offset=0,
        strand="+",
    )
    to_positions = [
        CdnaMappablePosition(
            _data=ensembl100,
            contig_id="13",
            start=67,
            start_offset=0,
            end=69,
            end_offset=0,
            strand="+",
            gene_id="ENSG00000139618",
            gene_name="BRCA2",
            transcript_id="ENST00000380152",
            transcript_name="BRCA2-201",
            protein_id="ENSP00000369497",
        ),
        CdnaMappablePosition(
            _data=ensembl100,
            contig_id="13",
            start=67,
            start_offset=0,
            end=69,
            end_offset=0,
            strand="+",
            gene_id="ENSG00000139618",
            gene_name="BRCA2",
            transcript_id="ENST00000544455",
            transcript_name="BRCA2-206",
            protein_id="ENSP00000439902",
        ),
    ]
    assert from_position.to_cdna() == to_positions


def test_to_cdna_dna_negative_offset_into_intron_minus_strand(ensembl100):
    from_position = DnaMappablePosition(
        _data=ensembl100,
        contig_id="5",
        start=1294666,
        start_offset=-1,
        end=1294666,
        end_offset=-1,
        strand="-",
    )
    to_positions = []
    assert from_position.to_cdna() == to_positions


def test_to_cdna_dna_overlapping_genes_different_strands_minus_strand(ensembl100):
    from_position = DnaMappablePosition(
        _data=ensembl100,
        contig_id="6",
        start=31166302,
        start_offset=0,
        end=31166304,
        end_offset=0,
        strand="-",
    )
    to_positions = []
    assert from_position.to_cdna() == to_positions


def test_to_cdna_dna_overlapping_genes_different_strands_plus_strand(ensembl100):
    from_position = DnaMappablePosition(
        _data=ensembl100,
        contig_id="6",
        start=31166302,
        start_offset=0,
        end=31166304,
        end_offset=0,
        strand="+",
    )
    to_positions = [
        CdnaMappablePosition(
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
    ]
    assert from_position.to_cdna() == to_positions


def test_to_cdna_dna_positive_offset_normalized_plus_strand(ensembl100):
    from_position = DnaMappablePosition(
        _data=ensembl100,
        contig_id="5",
        start=1294771,
        start_offset=1,
        end=1294771,
        end_offset=2,
        strand="+",
    )
    to_positions = []
    assert from_position.to_cdna() == to_positions


def test_to_cdna_dna_positive_offset_out_of_intron_minus_strand(ensembl100):
    from_position = DnaMappablePosition(
        _data=ensembl100,
        contig_id="5",
        start=1294667,
        start_offset=1,
        end=1294667,
        end_offset=1,
        strand="-",
    )
    to_positions = [
        CdnaMappablePosition(
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
        ),
        CdnaMappablePosition(
            _data=ensembl100,
            contig_id="5",
            start=220,
            start_offset=0,
            end=220,
            end_offset=0,
            strand="-",
            gene_id="ENSG00000164362",
            gene_name="TERT",
            transcript_id="ENST00000334602",
            transcript_name="TERT-202",
            protein_id="ENSP00000334346",
        ),
        CdnaMappablePosition(
            _data=ensembl100,
            contig_id="5",
            start=220,
            start_offset=0,
            end=220,
            end_offset=0,
            strand="-",
            gene_id="ENSG00000164362",
            gene_name="TERT",
            transcript_id="ENST00000460137",
            transcript_name="TERT-203",
            protein_id="ENSP00000425003",
        ),
        CdnaMappablePosition(
            _data=ensembl100,
            contig_id="5",
            start=220,
            start_offset=0,
            end=220,
            end_offset=0,
            strand="-",
            gene_id="ENSG00000164362",
            gene_name="TERT",
            transcript_id="ENST00000656021",
            transcript_name="TERT-206",
            protein_id="ENSP00000499759",
        ),
    ]
    assert from_position.to_cdna() == to_positions


def test_to_cdna_dna_start_minus_strand(ensembl100):
    from_position = DnaMappablePosition(
        _data=ensembl100,
        contig_id="5",
        start=1295066,
        start_offset=0,
        end=1295068,
        end_offset=0,
        strand="-",
    )
    to_positions = []
    assert from_position.to_cdna() == to_positions


def test_to_cdna_dna_start_plus_strand(ensembl100):
    from_position = DnaMappablePosition(
        _data=ensembl100,
        contig_id="13",
        start=32315508,
        start_offset=0,
        end=32315510,
        end_offset=0,
        strand="+",
    )
    to_positions = []
    assert from_position.to_cdna() == to_positions


def test_to_dna_dna_end_minus_strand(ensembl100):
    from_position = DnaMappablePosition(
        _data=ensembl100,
        contig_id="5",
        start=1253147,
        start_offset=0,
        end=1253149,
        end_offset=0,
        strand="-",
    )
    to_positions = [
        DnaMappablePosition(
            _data=ensembl100,
            contig_id="5",
            start=1253147,
            start_offset=0,
            end=1253149,
            end_offset=0,
            strand="-",
        )
    ]
    assert from_position.to_dna() == to_positions


def test_to_dna_dna_end_plus_strand(ensembl100):
    from_position = DnaMappablePosition(
        _data=ensembl100,
        contig_id="13",
        start=32400266,
        start_offset=0,
        end=32400268,
        end_offset=0,
        strand="+",
    )
    to_positions = [
        DnaMappablePosition(
            _data=ensembl100,
            contig_id="13",
            start=32400266,
            start_offset=0,
            end=32400268,
            end_offset=0,
            strand="+",
        )
    ]
    assert from_position.to_dna() == to_positions


def test_to_dna_dna_exon_boundary_minus_strand(ensembl100):
    from_position = DnaMappablePosition(
        _data=ensembl100,
        contig_id="5",
        start=1282623,
        start_offset=0,
        end=1293313,
        end_offset=0,
        strand="-",
    )
    to_positions = [
        DnaMappablePosition(
            _data=ensembl100,
            contig_id="5",
            start=1282623,
            start_offset=0,
            end=1293313,
            end_offset=0,
            strand="-",
        )
    ]
    assert from_position.to_dna() == to_positions


def test_to_dna_dna_exon_boundary_plus_strand(ensembl100):
    from_position = DnaMappablePosition(
        _data=ensembl100,
        contig_id="13",
        start=32316527,
        start_offset=0,
        end=32319078,
        end_offset=0,
        strand="+",
    )
    to_positions = [
        DnaMappablePosition(
            _data=ensembl100,
            contig_id="13",
            start=32316527,
            start_offset=0,
            end=32319078,
            end_offset=0,
            strand="+",
        )
    ]
    assert from_position.to_dna() == to_positions


def test_to_dna_dna_negative_offset_into_intron_minus_strand(ensembl100):
    from_position = DnaMappablePosition(
        _data=ensembl100,
        contig_id="5",
        start=1294666,
        start_offset=-1,
        end=1294666,
        end_offset=-1,
        strand="-",
    )
    to_positions = [
        DnaMappablePosition(
            _data=ensembl100,
            contig_id="5",
            start=1294667,
            start_offset=0,
            end=1294667,
            end_offset=0,
            strand="-",
        )
    ]
    assert from_position.to_dna() == to_positions


def test_to_dna_dna_overlapping_genes_different_strands_minus_strand(ensembl100):
    from_position = DnaMappablePosition(
        _data=ensembl100,
        contig_id="6",
        start=31166302,
        start_offset=0,
        end=31166304,
        end_offset=0,
        strand="-",
    )
    to_positions = [
        DnaMappablePosition(
            _data=ensembl100,
            contig_id="6",
            start=31166302,
            start_offset=0,
            end=31166304,
            end_offset=0,
            strand="-",
        )
    ]
    assert from_position.to_dna() == to_positions


def test_to_dna_dna_overlapping_genes_different_strands_plus_strand(ensembl100):
    from_position = DnaMappablePosition(
        _data=ensembl100,
        contig_id="6",
        start=31166302,
        start_offset=0,
        end=31166304,
        end_offset=0,
        strand="+",
    )
    to_positions = [
        DnaMappablePosition(
            _data=ensembl100,
            contig_id="6",
            start=31166302,
            start_offset=0,
            end=31166304,
            end_offset=0,
            strand="+",
        )
    ]
    assert from_position.to_dna() == to_positions


def test_to_dna_dna_positive_offset_normalized_plus_strand(ensembl100):
    from_position = DnaMappablePosition(
        _data=ensembl100,
        contig_id="5",
        start=1294771,
        start_offset=1,
        end=1294771,
        end_offset=2,
        strand="+",
    )
    to_positions = [
        DnaMappablePosition(
            _data=ensembl100,
            contig_id="5",
            start=1294772,
            start_offset=0,
            end=1294773,
            end_offset=0,
            strand="+",
        )
    ]
    assert from_position.to_dna() == to_positions


def test_to_dna_dna_positive_offset_out_of_intron_minus_strand(ensembl100):
    from_position = DnaMappablePosition(
        _data=ensembl100,
        contig_id="5",
        start=1294667,
        start_offset=1,
        end=1294667,
        end_offset=1,
        strand="-",
    )
    to_positions = [
        DnaMappablePosition(
            _data=ensembl100,
            contig_id="5",
            start=1294666,
            start_offset=0,
            end=1294666,
            end_offset=0,
            strand="-",
        )
    ]
    assert from_position.to_dna() == to_positions


def test_to_dna_dna_start_minus_strand(ensembl100):
    from_position = DnaMappablePosition(
        _data=ensembl100,
        contig_id="5",
        start=1295066,
        start_offset=0,
        end=1295068,
        end_offset=0,
        strand="-",
    )
    to_positions = [
        DnaMappablePosition(
            _data=ensembl100,
            contig_id="5",
            start=1295066,
            start_offset=0,
            end=1295068,
            end_offset=0,
            strand="-",
        )
    ]
    assert from_position.to_dna() == to_positions


def test_to_dna_dna_start_plus_strand(ensembl100):
    from_position = DnaMappablePosition(
        _data=ensembl100,
        contig_id="13",
        start=32315508,
        start_offset=0,
        end=32315510,
        end_offset=0,
        strand="+",
    )
    to_positions = [
        DnaMappablePosition(
            _data=ensembl100,
            contig_id="13",
            start=32315508,
            start_offset=0,
            end=32315510,
            end_offset=0,
            strand="+",
        )
    ]
    assert from_position.to_dna() == to_positions


def test_to_exon_dna_end_minus_strand(ensembl100):
    from_position = DnaMappablePosition(
        _data=ensembl100,
        contig_id="5",
        start=1253147,
        start_offset=0,
        end=1253149,
        end_offset=0,
        strand="-",
    )
    to_positions = [
        ExonMappablePosition(
            _data=ensembl100,
            contig_id="5",
            start=11,
            start_offset=0,
            end=11,
            end_offset=0,
            strand="-",
            gene_id="ENSG00000164362",
            gene_name="TERT",
            transcript_id="ENST00000484238",
            transcript_name="TERT-204",
            exon_id="ENSE00001871170",
        )
    ]
    assert from_position.to_exon() == to_positions


def test_to_exon_dna_end_plus_strand(ensembl100):
    from_position = DnaMappablePosition(
        _data=ensembl100,
        contig_id="13",
        start=32400266,
        start_offset=0,
        end=32400268,
        end_offset=0,
        strand="+",
    )
    to_positions = []
    assert from_position.to_exon() == to_positions


def test_to_exon_dna_exon_boundary_minus_strand(ensembl100):
    from_position = DnaMappablePosition(
        _data=ensembl100,
        contig_id="5",
        start=1282623,
        start_offset=0,
        end=1293313,
        end_offset=0,
        strand="-",
    )
    to_positions = [
        ExonMappablePosition(
            _data=ensembl100,
            contig_id="5",
            start=2,
            start_offset=0,
            end=3,
            end_offset=0,
            strand="-",
            gene_id="ENSG00000164362",
            gene_name="TERT",
            transcript_id="ENST00000310581",
            transcript_name="TERT-201",
            exon_id="ENSE00001197112",
        ),
        ExonMappablePosition(
            _data=ensembl100,
            contig_id="5",
            start=2,
            start_offset=0,
            end=3,
            end_offset=0,
            strand="-",
            gene_id="ENSG00000164362",
            gene_name="TERT",
            transcript_id="ENST00000334602",
            transcript_name="TERT-202",
            exon_id="ENSE00001197112",
        ),
        ExonMappablePosition(
            _data=ensembl100,
            contig_id="5",
            start=2,
            start_offset=0,
            end=3,
            end_offset=0,
            strand="-",
            gene_id="ENSG00000164362",
            gene_name="TERT",
            transcript_id="ENST00000460137",
            transcript_name="TERT-203",
            exon_id="ENSE00001197112",
        ),
        ExonMappablePosition(
            _data=ensembl100,
            contig_id="5",
            start=2,
            start_offset=0,
            end=4,
            end_offset=0,
            strand="-",
            gene_id="ENSG00000164362",
            gene_name="TERT",
            transcript_id="ENST00000656021",
            transcript_name="TERT-206",
            exon_id="ENSE00003863765",
        ),
    ]
    assert from_position.to_exon() == to_positions


def test_to_exon_dna_exon_boundary_plus_strand(ensembl100):
    from_position = DnaMappablePosition(
        _data=ensembl100,
        contig_id="13",
        start=32316527,
        start_offset=0,
        end=32319078,
        end_offset=0,
        strand="+",
    )
    to_positions = [
        ExonMappablePosition(
            _data=ensembl100,
            contig_id="13",
            start=2,
            start_offset=0,
            end=3,
            end_offset=0,
            strand="+",
            gene_id="ENSG00000139618",
            gene_name="BRCA2",
            transcript_id="ENST00000380152",
            transcript_name="BRCA2-201",
            exon_id="ENSE00001484009",
        ),
        ExonMappablePosition(
            _data=ensembl100,
            contig_id="13",
            start=2,
            start_offset=0,
            end=3,
            end_offset=0,
            strand="+",
            gene_id="ENSG00000139618",
            gene_name="BRCA2",
            transcript_id="ENST00000544455",
            transcript_name="BRCA2-206",
            exon_id="ENSE00001484009",
        ),
        ExonMappablePosition(
            _data=ensembl100,
            contig_id="13",
            start=1,
            start_offset=0,
            end=2,
            end_offset=0,
            strand="+",
            gene_id="ENSG00000139618",
            gene_name="BRCA2",
            transcript_id="ENST00000614259",
            transcript_name="BRCA2-207",
            exon_id="ENSE00003728983",
        ),
    ]
    assert from_position.to_exon() == to_positions


def test_to_exon_dna_negative_offset_into_intron_minus_strand(ensembl100):
    from_position = DnaMappablePosition(
        _data=ensembl100,
        contig_id="5",
        start=1294666,
        start_offset=-1,
        end=1294666,
        end_offset=-1,
        strand="-",
    )
    to_positions = []
    assert from_position.to_exon() == to_positions


def test_to_exon_dna_overlapping_genes_different_strands_minus_strand(ensembl100):
    from_position = DnaMappablePosition(
        _data=ensembl100,
        contig_id="6",
        start=31166302,
        start_offset=0,
        end=31166304,
        end_offset=0,
        strand="-",
    )
    to_positions = [
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
    assert from_position.to_exon() == to_positions


def test_to_exon_dna_overlapping_genes_different_strands_plus_strand(ensembl100):
    from_position = DnaMappablePosition(
        _data=ensembl100,
        contig_id="6",
        start=31166302,
        start_offset=0,
        end=31166304,
        end_offset=0,
        strand="+",
    )
    to_positions = [
        ExonMappablePosition(
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
    ]
    assert from_position.to_exon() == to_positions


def test_to_exon_dna_positive_offset_normalized_plus_strand(ensembl100):
    from_position = DnaMappablePosition(
        _data=ensembl100,
        contig_id="5",
        start=1294771,
        start_offset=1,
        end=1294771,
        end_offset=2,
        strand="+",
    )
    to_positions = []
    assert from_position.to_exon() == to_positions


def test_to_exon_dna_positive_offset_out_of_intron_minus_strand(ensembl100):
    from_position = DnaMappablePosition(
        _data=ensembl100,
        contig_id="5",
        start=1294667,
        start_offset=1,
        end=1294667,
        end_offset=1,
        strand="-",
    )
    to_positions = [
        ExonMappablePosition(
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
        ),
        ExonMappablePosition(
            _data=ensembl100,
            contig_id="5",
            start=2,
            start_offset=0,
            end=2,
            end_offset=0,
            strand="-",
            gene_id="ENSG00000164362",
            gene_name="TERT",
            transcript_id="ENST00000334602",
            transcript_name="TERT-202",
            exon_id="ENSE00001197112",
        ),
        ExonMappablePosition(
            _data=ensembl100,
            contig_id="5",
            start=2,
            start_offset=0,
            end=2,
            end_offset=0,
            strand="-",
            gene_id="ENSG00000164362",
            gene_name="TERT",
            transcript_id="ENST00000460137",
            transcript_name="TERT-203",
            exon_id="ENSE00001197112",
        ),
        ExonMappablePosition(
            _data=ensembl100,
            contig_id="5",
            start=2,
            start_offset=0,
            end=2,
            end_offset=0,
            strand="-",
            gene_id="ENSG00000164362",
            gene_name="TERT",
            transcript_id="ENST00000656021",
            transcript_name="TERT-206",
            exon_id="ENSE00003863765",
        ),
    ]
    assert from_position.to_exon() == to_positions


def test_to_exon_dna_start_minus_strand(ensembl100):
    from_position = DnaMappablePosition(
        _data=ensembl100,
        contig_id="5",
        start=1295066,
        start_offset=0,
        end=1295068,
        end_offset=0,
        strand="-",
    )
    to_positions = [
        ExonMappablePosition(
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
    ]
    assert from_position.to_exon() == to_positions


def test_to_exon_dna_start_plus_strand(ensembl100):
    from_position = DnaMappablePosition(
        _data=ensembl100,
        contig_id="13",
        start=32315508,
        start_offset=0,
        end=32315510,
        end_offset=0,
        strand="+",
    )
    to_positions = [
        ExonMappablePosition(
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
        ),
        ExonMappablePosition(
            _data=ensembl100,
            contig_id="13",
            start=1,
            start_offset=0,
            end=1,
            end_offset=0,
            strand="+",
            gene_id="ENSG00000139618",
            gene_name="BRCA2",
            transcript_id="ENST00000530893",
            transcript_name="BRCA2-204",
            exon_id="ENSE00002143308",
        ),
        ExonMappablePosition(
            _data=ensembl100,
            contig_id="13",
            start=1,
            start_offset=0,
            end=1,
            end_offset=0,
            strand="+",
            gene_id="ENSG00000139618",
            gene_name="BRCA2",
            transcript_id="ENST00000544455",
            transcript_name="BRCA2-206",
            exon_id="ENSE00002209788",
        ),
    ]
    assert from_position.to_exon() == to_positions


def test_to_protein_dna_end_minus_strand(ensembl100):
    from_position = DnaMappablePosition(
        _data=ensembl100,
        contig_id="5",
        start=1253147,
        start_offset=0,
        end=1253149,
        end_offset=0,
        strand="-",
    )
    to_positions = []
    assert from_position.to_protein() == to_positions


def test_to_protein_dna_end_plus_strand(ensembl100):
    from_position = DnaMappablePosition(
        _data=ensembl100,
        contig_id="13",
        start=32400266,
        start_offset=0,
        end=32400268,
        end_offset=0,
        strand="+",
    )
    to_positions = []
    assert from_position.to_protein() == to_positions


def test_to_protein_dna_exon_boundary_minus_strand(ensembl100):
    from_position = DnaMappablePosition(
        _data=ensembl100,
        contig_id="5",
        start=1282623,
        start_offset=0,
        end=1293313,
        end_offset=0,
        strand="-",
    )
    to_positions = [
        ProteinMappablePosition(
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
        ),
        ProteinMappablePosition(
            _data=ensembl100,
            contig_id="5",
            start=525,
            start_offset=0,
            end=525,
            end_offset=0,
            strand="-",
            gene_id="ENSG00000164362",
            gene_name="TERT",
            transcript_id="ENST00000334602",
            transcript_name="TERT-202",
            protein_id="ENSP00000334346",
        ),
        ProteinMappablePosition(
            _data=ensembl100,
            contig_id="5",
            start=525,
            start_offset=0,
            end=525,
            end_offset=0,
            strand="-",
            gene_id="ENSG00000164362",
            gene_name="TERT",
            transcript_id="ENST00000460137",
            transcript_name="TERT-203",
            protein_id="ENSP00000425003",
        ),
    ]
    assert from_position.to_protein() == to_positions


def test_to_protein_dna_exon_boundary_plus_strand(ensembl100):
    from_position = DnaMappablePosition(
        _data=ensembl100,
        contig_id="13",
        start=32316527,
        start_offset=0,
        end=32319078,
        end_offset=0,
        strand="+",
    )
    to_positions = [
        ProteinMappablePosition(
            _data=ensembl100,
            contig_id="13",
            start=23,
            start_offset=0,
            end=23,
            end_offset=0,
            strand="+",
            gene_id="ENSG00000139618",
            gene_name="BRCA2",
            transcript_id="ENST00000380152",
            transcript_name="BRCA2-201",
            protein_id="ENSP00000369497",
        ),
        ProteinMappablePosition(
            _data=ensembl100,
            contig_id="13",
            start=23,
            start_offset=0,
            end=23,
            end_offset=0,
            strand="+",
            gene_id="ENSG00000139618",
            gene_name="BRCA2",
            transcript_id="ENST00000544455",
            transcript_name="BRCA2-206",
            protein_id="ENSP00000439902",
        ),
    ]
    assert from_position.to_protein() == to_positions


def test_to_protein_dna_negative_offset_into_intron_minus_strand(ensembl100):
    from_position = DnaMappablePosition(
        _data=ensembl100,
        contig_id="5",
        start=1294666,
        start_offset=-1,
        end=1294666,
        end_offset=-1,
        strand="-",
    )
    to_positions = []
    assert from_position.to_protein() == to_positions


def test_to_protein_dna_overlapping_genes_different_strands_minus_strand(ensembl100):
    from_position = DnaMappablePosition(
        _data=ensembl100,
        contig_id="6",
        start=31166302,
        start_offset=0,
        end=31166304,
        end_offset=0,
        strand="-",
    )
    to_positions = []
    assert from_position.to_protein() == to_positions


def test_to_protein_dna_overlapping_genes_different_strands_plus_strand(ensembl100):
    from_position = DnaMappablePosition(
        _data=ensembl100,
        contig_id="6",
        start=31166302,
        start_offset=0,
        end=31166304,
        end_offset=0,
        strand="+",
    )
    to_positions = [
        ProteinMappablePosition(
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
    ]
    assert from_position.to_protein() == to_positions


def test_to_protein_dna_positive_offset_normalized_plus_strand(ensembl100):
    from_position = DnaMappablePosition(
        _data=ensembl100,
        contig_id="5",
        start=1294771,
        start_offset=1,
        end=1294771,
        end_offset=2,
        strand="+",
    )
    to_positions = []
    assert from_position.to_protein() == to_positions


def test_to_protein_dna_positive_offset_out_of_intron_minus_strand(ensembl100):
    from_position = DnaMappablePosition(
        _data=ensembl100,
        contig_id="5",
        start=1294667,
        start_offset=1,
        end=1294667,
        end_offset=1,
        strand="-",
    )
    to_positions = [
        ProteinMappablePosition(
            _data=ensembl100,
            contig_id="5",
            start=74,
            start_offset=0,
            end=74,
            end_offset=0,
            strand="-",
            gene_id="ENSG00000164362",
            gene_name="TERT",
            transcript_id="ENST00000310581",
            transcript_name="TERT-201",
            protein_id="ENSP00000309572",
        ),
        ProteinMappablePosition(
            _data=ensembl100,
            contig_id="5",
            start=74,
            start_offset=0,
            end=74,
            end_offset=0,
            strand="-",
            gene_id="ENSG00000164362",
            gene_name="TERT",
            transcript_id="ENST00000334602",
            transcript_name="TERT-202",
            protein_id="ENSP00000334346",
        ),
        ProteinMappablePosition(
            _data=ensembl100,
            contig_id="5",
            start=74,
            start_offset=0,
            end=74,
            end_offset=0,
            strand="-",
            gene_id="ENSG00000164362",
            gene_name="TERT",
            transcript_id="ENST00000460137",
            transcript_name="TERT-203",
            protein_id="ENSP00000425003",
        ),
        ProteinMappablePosition(
            _data=ensembl100,
            contig_id="5",
            start=74,
            start_offset=0,
            end=74,
            end_offset=0,
            strand="-",
            gene_id="ENSG00000164362",
            gene_name="TERT",
            transcript_id="ENST00000656021",
            transcript_name="TERT-206",
            protein_id="ENSP00000499759",
        ),
    ]
    assert from_position.to_protein() == to_positions


def test_to_protein_dna_start_minus_strand(ensembl100):
    from_position = DnaMappablePosition(
        _data=ensembl100,
        contig_id="5",
        start=1295066,
        start_offset=0,
        end=1295068,
        end_offset=0,
        strand="-",
    )
    to_positions = []
    assert from_position.to_protein() == to_positions


def test_to_protein_dna_start_plus_strand(ensembl100):
    from_position = DnaMappablePosition(
        _data=ensembl100,
        contig_id="13",
        start=32315508,
        start_offset=0,
        end=32315510,
        end_offset=0,
        strand="+",
    )
    to_positions = []
    assert from_position.to_protein() == to_positions


def test_to_rna_dna_end_minus_strand(ensembl100):
    from_position = DnaMappablePosition(
        _data=ensembl100,
        contig_id="5",
        start=1253147,
        start_offset=0,
        end=1253149,
        end_offset=0,
        strand="-",
    )
    to_positions = [
        RnaMappablePosition(
            _data=ensembl100,
            contig_id="5",
            start=2420,
            start_offset=0,
            end=2422,
            end_offset=0,
            strand="-",
            gene_id="ENSG00000164362",
            gene_name="TERT",
            transcript_id="ENST00000484238",
            transcript_name="TERT-204",
        )
    ]
    assert from_position.to_rna() == to_positions


def test_to_rna_dna_end_plus_strand(ensembl100):
    from_position = DnaMappablePosition(
        _data=ensembl100,
        contig_id="13",
        start=32400266,
        start_offset=0,
        end=32400268,
        end_offset=0,
        strand="+",
    )
    to_positions = []
    assert from_position.to_rna() == to_positions


def test_to_rna_dna_exon_boundary_minus_strand(ensembl100):
    from_position = DnaMappablePosition(
        _data=ensembl100,
        contig_id="5",
        start=1282623,
        start_offset=0,
        end=1293313,
        end_offset=0,
        strand="-",
    )
    to_positions = [
        RnaMappablePosition(
            _data=ensembl100,
            contig_id="5",
            start=1652,
            start_offset=0,
            end=1654,
            end_offset=0,
            strand="-",
            gene_id="ENSG00000164362",
            gene_name="TERT",
            transcript_id="ENST00000310581",
            transcript_name="TERT-201",
        ),
        RnaMappablePosition(
            _data=ensembl100,
            contig_id="5",
            start=1573,
            start_offset=0,
            end=1575,
            end_offset=0,
            strand="-",
            gene_id="ENSG00000164362",
            gene_name="TERT",
            transcript_id="ENST00000334602",
            transcript_name="TERT-202",
        ),
        RnaMappablePosition(
            _data=ensembl100,
            contig_id="5",
            start=1573,
            start_offset=0,
            end=1575,
            end_offset=0,
            strand="-",
            gene_id="ENSG00000164362",
            gene_name="TERT",
            transcript_id="ENST00000460137",
            transcript_name="TERT-203",
        ),
        RnaMappablePosition(
            _data=ensembl100,
            contig_id="5",
            start=1631,
            start_offset=0,
            end=2904,
            end_offset=0,
            strand="-",
            gene_id="ENSG00000164362",
            gene_name="TERT",
            transcript_id="ENST00000656021",
            transcript_name="TERT-206",
        ),
    ]
    assert from_position.to_rna() == to_positions


def test_to_rna_dna_exon_boundary_plus_strand(ensembl100):
    from_position = DnaMappablePosition(
        _data=ensembl100,
        contig_id="13",
        start=32316527,
        start_offset=0,
        end=32319078,
        end_offset=0,
        strand="+",
    )
    to_positions = [
        RnaMappablePosition(
            _data=ensembl100,
            contig_id="13",
            start=300,
            start_offset=0,
            end=302,
            end_offset=0,
            strand="+",
            gene_id="ENSG00000139618",
            gene_name="BRCA2",
            transcript_id="ENST00000380152",
            transcript_name="BRCA2-201",
        ),
        RnaMappablePosition(
            _data=ensembl100,
            contig_id="13",
            start=294,
            start_offset=0,
            end=296,
            end_offset=0,
            strand="+",
            gene_id="ENSG00000139618",
            gene_name="BRCA2",
            transcript_id="ENST00000544455",
            transcript_name="BRCA2-206",
        ),
        RnaMappablePosition(
            _data=ensembl100,
            contig_id="13",
            start=67,
            start_offset=0,
            end=69,
            end_offset=0,
            strand="+",
            gene_id="ENSG00000139618",
            gene_name="BRCA2",
            transcript_id="ENST00000614259",
            transcript_name="BRCA2-207",
        ),
    ]
    assert from_position.to_rna() == to_positions


def test_to_rna_dna_negative_offset_into_intron_minus_strand(ensembl100):
    from_position = DnaMappablePosition(
        _data=ensembl100,
        contig_id="5",
        start=1294666,
        start_offset=-1,
        end=1294666,
        end_offset=-1,
        strand="-",
    )
    to_positions = []
    assert from_position.to_rna() == to_positions


def test_to_rna_dna_overlapping_genes_different_strands_minus_strand(ensembl100):
    from_position = DnaMappablePosition(
        _data=ensembl100,
        contig_id="6",
        start=31166302,
        start_offset=0,
        end=31166304,
        end_offset=0,
        strand="-",
    )
    to_positions = [
        RnaMappablePosition(
            _data=ensembl100,
            contig_id="6",
            start=535,
            start_offset=0,
            end=537,
            end_offset=0,
            strand="-",
            gene_id="ENSG00000204531",
            gene_name="POU5F1",
            transcript_id="ENST00000471529",
            transcript_name="POU5F1-204",
        ),
        RnaMappablePosition(
            _data=ensembl100,
            contig_id="6",
            start=867,
            start_offset=0,
            end=869,
            end_offset=0,
            strand="-",
            gene_id="ENSG00000204531",
            gene_name="POU5F1",
            transcript_id="ENST00000513407",
            transcript_name="POU5F1-206",
        ),
    ]
    assert from_position.to_rna() == to_positions


def test_to_rna_dna_overlapping_genes_different_strands_plus_strand(ensembl100):
    from_position = DnaMappablePosition(
        _data=ensembl100,
        contig_id="6",
        start=31166302,
        start_offset=0,
        end=31166304,
        end_offset=0,
        strand="+",
    )
    to_positions = [
        RnaMappablePosition(
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
        )
    ]
    assert from_position.to_rna() == to_positions


def test_to_rna_dna_positive_offset_normalized_plus_strand(ensembl100):
    from_position = DnaMappablePosition(
        _data=ensembl100,
        contig_id="5",
        start=1294771,
        start_offset=1,
        end=1294771,
        end_offset=2,
        strand="+",
    )
    to_positions = []
    assert from_position.to_rna() == to_positions


def test_to_rna_dna_positive_offset_out_of_intron_minus_strand(ensembl100):
    from_position = DnaMappablePosition(
        _data=ensembl100,
        contig_id="5",
        start=1294667,
        start_offset=1,
        end=1294667,
        end_offset=1,
        strand="-",
    )
    to_positions = [
        RnaMappablePosition(
            _data=ensembl100,
            contig_id="5",
            start=299,
            start_offset=0,
            end=299,
            end_offset=0,
            strand="-",
            gene_id="ENSG00000164362",
            gene_name="TERT",
            transcript_id="ENST00000310581",
            transcript_name="TERT-201",
        ),
        RnaMappablePosition(
            _data=ensembl100,
            contig_id="5",
            start=220,
            start_offset=0,
            end=220,
            end_offset=0,
            strand="-",
            gene_id="ENSG00000164362",
            gene_name="TERT",
            transcript_id="ENST00000334602",
            transcript_name="TERT-202",
        ),
        RnaMappablePosition(
            _data=ensembl100,
            contig_id="5",
            start=220,
            start_offset=0,
            end=220,
            end_offset=0,
            strand="-",
            gene_id="ENSG00000164362",
            gene_name="TERT",
            transcript_id="ENST00000460137",
            transcript_name="TERT-203",
        ),
        RnaMappablePosition(
            _data=ensembl100,
            contig_id="5",
            start=278,
            start_offset=0,
            end=278,
            end_offset=0,
            strand="-",
            gene_id="ENSG00000164362",
            gene_name="TERT",
            transcript_id="ENST00000656021",
            transcript_name="TERT-206",
        ),
    ]
    assert from_position.to_rna() == to_positions


def test_to_rna_dna_start_minus_strand(ensembl100):
    from_position = DnaMappablePosition(
        _data=ensembl100,
        contig_id="5",
        start=1295066,
        start_offset=0,
        end=1295068,
        end_offset=0,
        strand="-",
    )
    to_positions = [
        RnaMappablePosition(
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
        )
    ]
    assert from_position.to_rna() == to_positions


def test_to_rna_dna_start_plus_strand(ensembl100):
    from_position = DnaMappablePosition(
        _data=ensembl100,
        contig_id="13",
        start=32315508,
        start_offset=0,
        end=32315510,
        end_offset=0,
        strand="+",
    )
    to_positions = [
        RnaMappablePosition(
            _data=ensembl100,
            contig_id="13",
            start=35,
            start_offset=0,
            end=37,
            end_offset=0,
            strand="+",
            gene_id="ENSG00000139618",
            gene_name="BRCA2",
            transcript_id="ENST00000380152",
            transcript_name="BRCA2-201",
        ),
        RnaMappablePosition(
            _data=ensembl100,
            contig_id="13",
            start=4,
            start_offset=0,
            end=6,
            end_offset=0,
            strand="+",
            gene_id="ENSG00000139618",
            gene_name="BRCA2",
            transcript_id="ENST00000530893",
            transcript_name="BRCA2-204",
        ),
        RnaMappablePosition(
            _data=ensembl100,
            contig_id="13",
            start=29,
            start_offset=0,
            end=31,
            end_offset=0,
            strand="+",
            gene_id="ENSG00000139618",
            gene_name="BRCA2",
            transcript_id="ENST00000544455",
            transcript_name="BRCA2-206",
        ),
    ]
    assert from_position.to_rna() == to_positions
