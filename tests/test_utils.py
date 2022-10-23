from ensembl_map.utils import reverse_complement, strip_version


def test_reverse_complement():
    assert reverse_complement("AGCT") == "AGCT"


def test_strip_version_with_ver():
    assert strip_version("NM_000546.5") == "NM_000546"


def test_strip_version_no_ver():
    assert strip_version("NM_000546") == "NM_000546"
