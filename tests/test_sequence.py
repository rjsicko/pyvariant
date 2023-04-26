import pytest

from variant_map.sequence import get_sequence, mutate_sequence


@pytest.mark.parametrize(
    "transcript, start, end, window, sequence",
    [
        ("ENST00000288135", 7, 9, -1, "GGG"),
        ("ENST00000288135", 7, 9, 3, "GGG"),
        ("ENST00000288135", 7, 9, 6, "TGGGCG"),
        ("ENST00000288135", 7, 9, 7, "TTGGGCG"),
    ],
)
def test_get_sequence(ensembl100, transcript, start, end, window, sequence):
    assert get_sequence(ensembl100.rna_fasta[0], transcript, start, end, window) == sequence


def test_get_sequence_window_error(ensembl100):
    with pytest.raises(ValueError):
        get_sequence(ensembl100.rna_fasta[0], "ENST00000257430", 290, 291, 1)


@pytest.mark.parametrize(
    "transcript, start, end, window, altseq, insertion, sequence",
    [
        ("ENST00000288135", 7, 7, 7, "", False, "CTTGGCG"),  # ENST00000288135:r.7delG
        ("ENST00000288135", 7, 8, 7, "", False, "CTTGCGA"),  # ENST00000288135:r.7_8delGG
        ("ENST00000288135", 7, 9, 7, "", False, "CTTCGAG"),  # ENST00000288135:r.7_9delGGG
        ("ENST00000288135", 7, 8, 7, "T", True, "TTGTGGC"),  # ENST00000288135:r.7_8insT
        ("ENST00000288135", 7, 8, 7, "TT", True, "TGTTGGC"),  # ENST00000288135:r.7_8insTT
        ("ENST00000288135", 7, 8, 7, "TTT", True, "TGTTTGG"),  # ENST00000288135:r.7_8insTTT
        ("ENST00000288135", 7, 7, 7, "C", False, "CTTCGGC"),  # ENST00000288135:r.7G>C
        ("ENST00000288135", 7, 9, 0, "", False, ""),  # ENST00000288135:r.7_9delGGG
        ("ENST00000288135", 7, 9, 1, "", False, "C"),  # ENST00000288135:r.7_9delGGG
    ],
)
def test_mutate_sequence(ensembl100, transcript, start, end, window, altseq, insertion, sequence):
    assert (
        mutate_sequence(
            ensembl100.rna_fasta[0], transcript, start, end, window, altseq, insertion=insertion
        )
        == sequence
    )


def test_mutate_sequence_window_error(ensembl100):
    with pytest.raises(ValueError):
        mutate_sequence(
            ensembl100.rna_fasta[0], "ENST00000257430", 290, 291, 1, "TT", insertion=True
        )
