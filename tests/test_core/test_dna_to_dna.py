from ensembl_map.core import DnaPosition


def test_negative_strand(ensembl100):
    position = DnaPosition(_data=ensembl100, contig_id="5", start=1294864, end=1294866, strand="-")
    assert position.to_dna() == [position]


def test_negative_strand_across_exon_boundary(ensembl100):
    position = DnaPosition(_data=ensembl100, contig_id="5", start=1282623, end=1293313, strand="-")
    assert position.to_dna() == [position]


def test_negative_strand_cdna_protein_end(ensembl100):
    position = DnaPosition(_data=ensembl100, contig_id="5", start=1253728, end=1253730, strand="-")
    assert position.to_dna() == [position]


def test_negative_strand_cdna_protein_start(ensembl100):
    position = DnaPosition(_data=ensembl100, contig_id="5", start=1294987, end=1294989, strand="-")
    assert position.to_dna() == [position]


def test_negative_strand_overlapping_genes_different_strands(ensembl100):
    position = DnaPosition(
        _data=ensembl100, contig_id="6", start=31166302, end=31166304, strand="-"
    )
    assert position.to_dna() == [position]


def test_negative_strand_transcript_end(ensembl100):
    position = DnaPosition(_data=ensembl100, contig_id="5", start=1253167, end=1253169, strand="-")
    assert position.to_dna() == [position]


def test_negative_strand_transcript_start(ensembl100):
    position = DnaPosition(_data=ensembl100, contig_id="5", start=1295066, end=1295068, strand="-")
    assert position.to_dna() == [position]


def test_positive_strand(ensembl100):
    position = DnaPosition(
        _data=ensembl100, contig_id="13", start=32371034, end=32371036, strand="+"
    )
    assert position.to_dna() == [position]


def test_positive_strand_across_exon_boundary(ensembl100):
    position = DnaPosition(
        _data=ensembl100, contig_id="4", start=54658081, end=54695513, strand="+"
    )
    assert position.to_dna() == [position]


def test_positive_strand_cdna_protein_end(ensembl100):
    position = DnaPosition(
        _data=ensembl100, contig_id="13", start=32398768, end=32398770, strand="+"
    )
    assert position.to_dna() == [position]


def test_positive_strand_cdna_protein_start(ensembl100):
    position = DnaPosition(
        _data=ensembl100, contig_id="13", start=32316461, end=32316463, strand="+"
    )
    assert position.to_dna() == [position]


def test_positive_strand_overlapping_genes_different_strands(ensembl100):
    position = DnaPosition(
        _data=ensembl100, contig_id="6", start=31166302, end=31166304, strand="+"
    )
    assert position.to_dna() == [position]


def test_positive_strand_transcript_end(ensembl100):
    position = DnaPosition(
        _data=ensembl100, contig_id="13", start=32400264, end=32400266, strand="+"
    )
    assert position.to_dna() == [position]


def test_positive_strand_transcript_start(ensembl100):
    position = DnaPosition(
        _data=ensembl100, contig_id="13", start=32315474, end=32315476, strand="+"
    )
    assert position.to_dna() == [position]
