import pytest
from constants import TEST_ENS100_CDNA_FASTA
from pyfaidx import Fasta

from pyvariant.sequence import PyfaidxFasta, get_sequence, mutate_sequence


def test_pyfaidx_fasta_load():
    result = PyfaidxFasta.load(TEST_ENS100_CDNA_FASTA)
    assert isinstance(result, PyfaidxFasta)
    assert isinstance(result.fasta, Fasta)
    assert result["ENST00000256078"]


def test_pyfaidx_fasta_load_empty():
    result = PyfaidxFasta.load("")
    assert isinstance(result, PyfaidxFasta)
    assert isinstance(result.fasta, Fasta)


@pytest.mark.parametrize(
    "contig, start, end, strand, window, floor, ceiling, sequence",
    [
        # on positive strand
        ("17", 7675238, 7675238, "+", 7, -1, -1, "TACTGTA"),
        # on reverse strand
        ("17", 7675238, 7675238, "-", 7, -1, -1, "TACAGTA"),
    ],
)
def test_get_sequence_dna(ensembl100, contig, start, end, strand, window, floor, ceiling, sequence):
    assert (
        get_sequence(ensembl100.dna_fasta[0], contig, start, end, strand, window, floor, ceiling)
        == sequence
    )


@pytest.mark.parametrize(
    "transcript, start, end, strand, window, floor, ceiling, sequence",
    [
        ("ENST00000288135", 7, 9, "+", -1, -1, -1, "GGG"),
        ("ENST00000288135", 7, 9, "+", 3, -1, -1, "GGG"),
        ("ENST00000288135", 7, 9, "+", 6, -1, -1, "TGGGCG"),
        ("ENST00000288135", 7, 9, "+", 7, -1, -1, "TTGGGCG"),
        ("ENST00000288135", 7, 9, "+", 7, 7, -1, "GGGCGAG"),
        ("ENST00000288135", 7, 9, "+", 7, -1, 9, "ACTTGGG"),
        # at transcript start
        ("ENST00000288135", 1, 1, "+", 7, -1, -1, "GCACTTG"),
        # at transcript end
        ("ENST00000288135", 5147, 5147, "+", 7, -1, -1, "AAATCAA"),
    ],
)
def test_get_sequence_rna(
    ensembl100, transcript, start, end, strand, window, floor, ceiling, sequence
):
    assert (
        get_sequence(
            ensembl100.rna_fasta[0], transcript, start, end, strand, window, floor, ceiling
        )
        == sequence
    )


def test_get_sequence_window_error(ensembl100):
    with pytest.raises(ValueError):
        get_sequence(ensembl100.rna_fasta[0], "ENST00000257430", 290, 291, "+", 1, -1, -1)


@pytest.mark.parametrize(
    "contig, start, end, strand, window, floor, ceiling, altseq, insertion, sequence",
    [
        # on positive strand
        ("17", 7675238, 7675238, "+", 7, -1, -1, "G", False, "TACGGTA"),
        # on reverse strand
        ("17", 7675238, 7675238, "-", 7, -1, -1, "G", False, "TACGGTA"),
    ],
)
def test_mutate_sequence_dna(
    ensembl100, contig, start, end, strand, window, floor, ceiling, altseq, insertion, sequence
):
    assert (
        mutate_sequence(
            ensembl100.dna_fasta[0],
            contig,
            start,
            end,
            strand,
            window,
            floor,
            ceiling,
            altseq,
            insertion=insertion,
        )
        == sequence
    )


@pytest.mark.parametrize(
    "transcript, start, end, strand, window, floor, ceiling, altseq, sequence",
    [
        # ENST00000288135:r.7delG
        ("ENST00000288135", 7, 7, "+", 7, -1, -1, "", "CTTGGCG"),
        # ENST00000288135:r.7_8delGG
        ("ENST00000288135", 7, 8, "+", 7, -1, -1, "", "CTTGCGA"),
        # ENST00000288135:r.7_9delGGG
        ("ENST00000288135", 7, 9, "+", 7, -1, -1, "", "CTTCGAG"),
        # ENST00000288135:r.7_8insT
        ("ENST00000288135", 7, 8, "+", 7, -1, -1, "GTG", "TTGTGGC"),
        # ENST00000288135:r.7_8insTT
        ("ENST00000288135", 7, 8, "+", 7, -1, -1, "GTTG", "TGTTGGC"),
        # ENST00000288135:r.7_8insTTT
        ("ENST00000288135", 7, 8, "+", 7, -1, -1, "GTTTG", "TGTTTGG"),
        # ENST00000288135:r.7G>C
        ("ENST00000288135", 7, 7, "+", 7, -1, -1, "C", "CTTCGGC"),
        # ENST00000288135:r.7_9delGGG
        ("ENST00000288135", 7, 9, "+", 0, -1, -1, "", ""),
        # ENST00000288135:r.7_9delGGG
        ("ENST00000288135", 7, 9, "+", 1, -1, -1, "", "C"),
        ("ENST00000288135", 7, 9, "+", 7, 7, -1, "", "CGAGAGC"),
        ("ENST00000288135", 7, 9, "+", 6, -1, 9, "", "GCACTT"),
        # ENST00000288135:r.7_8insTTT
        ("ENST00000288135", 7, 8, "+", 7, 7, -1, "GTTTG", "GTTTGGC"),
        ("ENST00000288135", 7, 8, "+", 7, -1, 8, "GTTTG", "TTGTTTG"),
        # ENST00000288135:r.7_8delinsTTT
        ("ENST00000288135", 7, 8, "+", 7, -1, -1, "TTT", "TTTTTGC"),
        # ENST00000288135:r.7G>C
        ("ENST00000288135", 7, 7, "+", 7, 7, -1, "C", "CGGCGAG"),
        ("ENST00000288135", 7, 7, "+", 7, -1, 7, "C", "GCACTTC"),
        # near transcript start
        ("ENST00000288135", 2, 2, "+", 7, -1, -1, "A", "GAACTTG"),
        # at transcript end
        ("ENST00000288135", 5147, 5147, "+", 7, -1, -1, "C", "AAATCAC"),
    ],
)
def test_mutate_sequence_rna(
    ensembl100, transcript, start, end, strand, window, floor, ceiling, altseq, sequence
):
    assert (
        mutate_sequence(
            ensembl100.rna_fasta[0], transcript, start, end, strand, window, floor, ceiling, altseq
        )
        == sequence
    )


def test_mutate_sequence_window_error(ensembl100):
    with pytest.raises(ValueError):
        mutate_sequence(
            ensembl100.rna_fasta[0],
            "ENST00000257430",
            290,
            291,
            "+",
            1,
            -1,
            -1,
            "TT",
            insertion=True,
        )
