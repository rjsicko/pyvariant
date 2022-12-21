import pytest

from ensembl_map.utils import (
    collapse_mutation,
    expand_nt,
    expand_pep,
    format_hgvs_position,
    reverse_complement,
    reverse_translate,
    split_by_codon,
    strip_version,
)


def test_collapse_mutation():
    assert collapse_mutation("GGT", "GTG") == ("GT", "TG", 1, 0)
    assert collapse_mutation("GGT", "GTT") == ("G", "T", 1, 1)
    assert collapse_mutation("GTT", "CTT") == ("G", "C", 0, 2)
    assert collapse_mutation("TTG", "TTA") == ("G", "A", 2, 0)
    assert collapse_mutation("ATG", "A") == ("ATG", "A", 0, 0)
    assert collapse_mutation("ATG", "ATG") == ("ATG", "ATG", 0, 0)


def test_expand_nt():
    assert list(expand_nt("A")) == ["A"]
    assert list(expand_nt("R")) == ["A", "G"]
    assert list(expand_nt("N")) == ["A", "C", "G", "T"]


def test_expand_pep():
    assert list(expand_pep("A")) == ["A"]
    assert list(expand_pep("Z")) == ["E", "Q"]
    assert len(list(expand_pep("X"))) == 20  # all amino acids except '*'


def test_format_hgvs_position():
    assert format_hgvs_position(start=5) == "5"
    assert format_hgvs_position(start=5, start_offset=1) == "5+1"
    assert format_hgvs_position(start=5, start_offset=-1) == "5-1"
    assert format_hgvs_position(start=1, start_offset=1) == "+1"
    assert format_hgvs_position(start=1, start_offset=-1) == "-1"
    assert format_hgvs_position(end=5) == "5"
    assert format_hgvs_position(end=5, end_offset=1) == "5+1"
    assert format_hgvs_position(end=5, end_offset=-1) == "5-1"
    assert format_hgvs_position(end=1, end_offset=1) == "+1"
    assert format_hgvs_position(end=1, end_offset=-1) == "-1"
    assert format_hgvs_position(start=5, end=6) == "(5_6)"
    assert format_hgvs_position(start=5, end=6, start_offset=1) == "(5+1_6)"
    assert format_hgvs_position(start=5, end=6, end_offset=1) == "(5_6+1)"
    assert format_hgvs_position(start=5, end=6, start_offset=1, end_offset=1) == "(5+1_6+1)"


def test_reverse_complement():
    assert reverse_complement("AGCT") == "AGCT"


def test_reverse_translate():
    assert list(reverse_translate("A")) == ["GCA", "GCC", "GCG", "GCT"]
    assert list(reverse_translate("CQ")) == ["TGCCAA", "TGCCAG", "TGTCAA", "TGTCAG"]
    assert list(reverse_translate("B")) == ["GAC", "GAT", "AAC", "AAT"]


def test_split_by_codon():
    assert list(split_by_codon("ABC")) == ["ABC"]
    assert list(split_by_codon("ABCDEF")) == ["ABC", "DEF"]
    with pytest.raises(ValueError):
        # Error raised when iterable is not divisible by 3
        list(split_by_codon("ABCD"))


def test_strip_version():
    assert strip_version("NM_000546.5") == "NM_000546"
    assert strip_version("NM_000546") == "NM_000546"
