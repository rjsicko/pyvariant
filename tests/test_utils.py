import pytest

from ensembl_map.utils import (
    collapse_mutation,
    collapse_mutation_ends,
    expand_nt,
    expand_pep,
    format_hgvs_position,
    is_frameshift,
    is_insertion,
    reverse_complement,
    reverse_translate,
    split_by_codon,
    split_insertion,
    strip_version,
)


def test_collapse_mutation():
    assert collapse_mutation("GGT", "GTT") == ("G", "T", 1, 1)  # substitution
    assert collapse_mutation("GTT", "CTT") == ("G", "C", 0, 2)  # substitution
    assert collapse_mutation("TTG", "TTA") == ("G", "A", 2, 0)  # substitution
    assert collapse_mutation("ATG", "A") == ("TG", "", 1, 0)  # deletion
    assert collapse_mutation("ATTTA", "AA") == ("TTT", "", 1, 1)  # deletion
    assert collapse_mutation("AATTTC", "AAC") == ("TTT", "", 2, 1)  # deletion
    assert collapse_mutation("CTTTAA", "CA") == ("TTTA", "", 1, 1)  # deletion
    assert collapse_mutation("A", "AA") == ("A", "AA", 0, 0)  # duplication
    assert collapse_mutation("ATG", "ATGATG") == ("ATG", "ATGATG", 0, 0)  # duplication
    assert collapse_mutation("TCT", "AGCTCT") == ("TCT", "AGCTCT", 0, 0)  # amino acid duplication
    assert collapse_mutation("AA", "ATTTA") == ("AA", "ATTTA", 0, 0)  # insertion
    assert collapse_mutation("AAC", "AATTTC") == ("AC", "ATTTC", 1, 0)  # insertion
    assert collapse_mutation("CAA", "CTTTAA") == ("CA", "CTTTA", 0, 1)  # insertion
    assert collapse_mutation("GGT", "GTG") == ("GT", "TG", 1, 0)  # delins
    assert collapse_mutation("ATCG", "ATT") == ("CG", "T", 2, 0)  # delins
    assert collapse_mutation("ATCG", "TTG") == ("ATC", "TT", 0, 1)  # delins
    assert collapse_mutation("ATG", "ATG") == ("ATG", "ATG", 0, 0)  # synonymous


def test_collapse_mutation_ends():
    assert collapse_mutation_ends("AATTTC", "AAGC") == ("TTT", "G", "AA", "C")
    assert collapse_mutation_ends("ATTTC", "AC") == ("TTT", "", "A", "C")
    assert collapse_mutation_ends("CGAA", "CTTTAA") == ("G", "TTT", "C", "AA")
    assert collapse_mutation_ends("ATTTC", "AC") == ("TTT", "", "A", "C")
    assert collapse_mutation_ends("ATG", "GTA") == ("ATG", "GTA", "", "")
    assert collapse_mutation_ends("GCTGGT", "GCTACTGGT") == ("", "ACT", "GCT", "GGT")


def test_expand_nt():
    assert list(expand_nt("A")) == ["A"]
    assert list(expand_nt("R")) == ["A", "G"]
    assert list(expand_nt("N")) == ["A", "C", "G", "T"]


def test_expand_pep():
    assert list(expand_pep("A")) == ["A"]
    assert list(expand_pep("Z")) == ["E", "Q"]
    assert len(list(expand_pep("X"))) == 20  # all amino acids except '*'


def test_format_hgvs_position():
    assert format_hgvs_position(5, 0) == "5"
    assert format_hgvs_position(5, 1) == "5+1"
    assert format_hgvs_position(5, -1) == "5-1"
    assert format_hgvs_position(1, 1) == "+1"
    assert format_hgvs_position(1, -1) == "-1"
    assert format_hgvs_position(5, 0, is_3_prime_utr=True) == "*5"
    assert format_hgvs_position(5, 1, is_3_prime_utr=True) == "*5+1"
    assert format_hgvs_position(5, -1, is_3_prime_utr=True) == "*5-1"
    assert format_hgvs_position(1, 1, is_3_prime_utr=True) == "*1+1"
    assert format_hgvs_position(1, -1, is_3_prime_utr=True) == "*1-1"


def test_is_frameshift():
    assert is_frameshift("AAA", "A")
    assert is_frameshift("AAA", "AA")
    assert not is_frameshift("AAA", "AAA")


def test_is_insertion():
    assert is_insertion("GT", "GAT")
    assert is_insertion("GGGTTT", "GGGAAATTT")
    assert not is_insertion("GT", "GAA")
    assert not is_insertion("GT", "G")
    assert not is_insertion("GT", "T")
    assert not is_insertion("", "")


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


def test_split_insertion():
    assert split_insertion("GT", "GAT") == ("G", "A", "T")
    assert split_insertion("GGGTTT", "GGGAAATTT") == ("GGG", "AAA", "TTT")
    assert split_insertion("GT", "GAA") is None
    assert split_insertion("GT", "G") is None
    assert split_insertion("GT", "T") is None
    assert split_insertion("", "") is None


def test_strip_version():
    assert strip_version("NM_000546.5") == "NM_000546"
    assert strip_version("NM_000546") == "NM_000546"
