from ensembl_map.sequences import (
    cdna_sequence,
    contig_sequence,
    gene_sequence,
    protein_sequence,
    transcript_sequence,
)


def test_cdna_sequence():
    assert cdna_sequence("ENST00000288135", 1) == "A"


def test_cdna_sequence_2():
    assert cdna_sequence("ENST00000288135", 1, 3) == "ATG"


def test_contig_sequence():
    assert contig_sequence("12", 25378706) == "A"


def test_contig_sequence_2():
    assert contig_sequence("12", 25378704, 25378706) == "TAA"


def test_gene_sequence():
    assert gene_sequence("ENSG00000133703", 25205346) == "A"


def test_gene_sequence_2():
    assert gene_sequence("ENSG00000133703", 25205346, 25205348) == "TTA"


def test_protein_sequence():
    assert protein_sequence("ENSP00000288135", 1) == "M"


def test_protein_sequence_2():
    assert protein_sequence("ENSP00000288135", 3, 7) == "GARGA"


def test_transcript_sequence():
    assert transcript_sequence("ENST00000288135", 1) == "G"


def test_transcript_sequence_2():
    assert transcript_sequence("ENST00000288135", 166, 168) == "ACC"
