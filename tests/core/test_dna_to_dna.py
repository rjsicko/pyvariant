from ensembl_map.core import DnaMappablePosition


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
    assert position.to_dna() == [position]


def test_negative_strand_across_exon_boundary(ensembl100):
    position = DnaMappablePosition(
        _data=ensembl100,
        contig_id="5",
        start=1282623,
        start_offset=0,
        end=1293313,
        end_offset=0,
        strand="-",
    )
    assert position.to_dna() == [position]


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
    assert position.to_dna() == [position]


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
    assert position.to_dna() == [position]


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
    assert position.to_dna() == [position]


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
    assert position.to_dna() == [position]


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
    assert position.to_dna() == [position]


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
    assert position.to_dna() == [position]


def test_positive_strand_across_exon_boundary(ensembl100):
    position = DnaMappablePosition(
        _data=ensembl100,
        contig_id="4",
        start=54658081,
        start_offset=0,
        end=54695513,
        end_offset=0,
        strand="+",
    )
    assert position.to_dna() == [position]


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
    assert position.to_dna() == [position]


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
    assert position.to_dna() == [position]


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
    assert position.to_dna() == [position]


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
    assert position.to_dna() == [position]


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
    assert position.to_dna() == [position]


def test_offset_normalized(ensembl100):
    position = DnaMappablePosition(
        _data=ensembl100,
        contig_id="5",
        start=1294771,
        start_offset=1,
        end=1294771,
        end_offset=2,
        strand="+",
    )
    expected = DnaMappablePosition(
        _data=ensembl100,
        contig_id="5",
        start=1294772,
        start_offset=0,
        end=1294773,
        end_offset=0,
        strand="+",
    )
    assert position.to_dna() == [expected]
