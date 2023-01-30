from variant_map.core import DnaPosition


def test_is_cdna(ensembl100):
    position = DnaPosition(
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
    position = DnaPosition(
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
    position = DnaPosition(
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
    position = DnaPosition(
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
    position = DnaPosition(
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
    position = DnaPosition(
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
    position = DnaPosition(
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
    position = DnaPosition(
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


def test_str_mutli_position(ensembl100):
    position = DnaPosition(
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
    position = DnaPosition(
        _data=ensembl100,
        contig_id="5",
        start=1282623,
        start_offset=0,
        end=1282623,
        end_offset=0,
        strand="-",
    )
    assert str(position) == "5:g.1282623"
