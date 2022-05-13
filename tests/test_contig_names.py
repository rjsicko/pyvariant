from coordinate_mapper.contig_names import ContigNames

from . import CONTIG_FILE


def test_normalize():
    ContigNames.load(CONTIG_FILE)
    assert ContigNames.normalize("chr3") == "3"


def test_normalize_unknown():
    ContigNames.load("")
    assert ContigNames.normalize("chr3") is None
