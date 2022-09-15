import pytest

from ensembl_map.constants import CONTIG
from ensembl_map.core import EnsemblRelease
from ensembl_map.map import (
    cdna_to_cdna,
    cdna_to_dna,
    cdna_to_exon,
    cdna_to_protein,
    cdna_to_rna,
    dna_to_cdna,
    dna_to_dna,
    dna_to_exon,
    dna_to_protein,
    dna_to_rna,
    exon_to_cdna,
    exon_to_dna,
    exon_to_exon,
    exon_to_protein,
    exon_to_rna,
    protein_to_cdna,
    protein_to_dna,
    protein_to_exon,
    protein_to_protein,
    protein_to_rna,
    rna_to_cdna,
    rna_to_dna,
    rna_to_exon,
    rna_to_protein,
    rna_to_rna,
)

from . import (
    CACHE_DIR,
    CANONICAL_TRANSCRIPT,
    CONTIG_ALIAS,
    EXON_ALIAS,
    GENE_ALIAS,
    PROTEIN_ALIAS,
    TRANSCRIPT_ALIAS,
)


@pytest.fixture(autouse=True, scope="module")
def ensembl100():
    return EnsemblRelease(
        species="homo_sapiens",
        reference="GRCh38",
        release=100,
        cache_dir=CACHE_DIR,
        canonical_transcript=CANONICAL_TRANSCRIPT,
        contig_alias=CONTIG_ALIAS,
        exon_alias=EXON_ALIAS,
        gene_alias=GENE_ALIAS,
        protein_alias=PROTEIN_ALIAS,
        transcript_alias=TRANSCRIPT_ALIAS,
    )


def test_cdna_to_cdna_pos_strand_by_transcript_id():
    pos = cdna_to_cdna("ENST00000380152", 1)
    assert pos
    assert isinstance(pos, list)
    assert pos[0].transcript_id, "ENST00000380152"
    assert pos[0].start == 1
    assert pos[0].end == 1
    assert pos[0].strand == "+"


def test_cdna_to_cdna_pos_strand_by_transcript_id_2():
    pos = cdna_to_cdna("ENST00000380152", 100)
    assert pos
    assert isinstance(pos, list)
    assert pos[0].transcript_id, "ENST00000380152"
    assert pos[0].start == 100
    assert pos[0].end == 100
    assert pos[0].strand == "+"


def test_cdna_to_cdna_pos_strand_by_transcript_id_3():
    pos = cdna_to_cdna("ENST00000380152", 100, 101)
    assert pos
    assert isinstance(pos, list)
    assert pos[0].transcript_id, "ENST00000380152"
    assert pos[0].start == 100
    assert pos[0].end == 101
    assert pos[0].strand == "+"


def test_cdna_to_cdna_pos_strand_by_transcript_name():
    pos = cdna_to_cdna("BRCA2-201", 1)
    assert pos
    assert isinstance(pos, list)
    assert pos[0].transcript_id, "ENST00000380152"
    assert pos[0].start == 1
    assert pos[0].end == 1
    assert pos[0].strand == "+"


def test_cdna_to_cdna_pos_strand_by_transcript_name_2():
    pos = cdna_to_cdna("BRCA2-201", 360)
    assert pos
    assert isinstance(pos, list)
    assert pos[0].transcript_id, "ENST00000380152"
    assert pos[0].start == 360
    assert pos[0].end == 360
    assert pos[0].strand == "+"


def test_cdna_to_cdna_pos_strand_by_transcript_name_3():
    pos = cdna_to_cdna("BRCA2-201", 360, 363)
    assert pos
    assert isinstance(pos, list)
    assert pos[0].transcript_id, "ENST00000380152"
    assert pos[0].start == 360
    assert pos[0].end == 363
    assert pos[0].strand == "+"


def test_cdna_to_cdna_by_contig_id():
    # NOTE: I don't see this query being useful because it will return coordinates for every transcript on chr12.
    # result = sorted([i.transcript_id for i in cdna_to_cdna("12", 1)])
    pass


def test_cdna_to_cdna_by_exon_id():
    result = sorted([i.transcript_id for i in cdna_to_cdna("ENSE00000936617", 1)])
    assert result == ["ENST00000256078", "ENST00000311936", "ENST00000556131", "ENST00000557334"]


def test_cdna_to_cdna_by_gene_id():
    result = sorted([i.transcript_id for i in cdna_to_cdna("ENSG00000133703", 1)])
    assert result == ["ENST00000256078", "ENST00000311936", "ENST00000556131", "ENST00000557334"]


def test_cdna_to_cdna_by_gene_name():
    result = sorted([i.transcript_id for i in cdna_to_cdna("KRAS", 1)])
    assert result == ["ENST00000256078", "ENST00000311936", "ENST00000556131", "ENST00000557334"]


def test_cdna_to_cdna_by_protein_id():
    result = sorted([i.transcript_id for i in cdna_to_cdna("ENSP00000256078", 1)])
    assert result == ["ENST00000256078"]


def test_cdna_to_dna_pos_strand_by_transcript_id():
    pos = cdna_to_dna("ENST00000380152", 1)
    assert pos
    assert isinstance(pos, list)
    assert pos[0].contig_id == "13"
    assert pos[0].start == 32316461
    assert pos[0].end == 32316461
    assert pos[0].strand == "+"


def test_cdna_to_dna_pos_strand_by_transcript_id_2():
    pos = cdna_to_dna("ENST00000380152", 100)
    assert pos
    assert isinstance(pos, list)
    assert pos[0].contig_id == "13"
    assert pos[0].start == 32319109
    assert pos[0].end == 32319109
    assert pos[0].strand == "+"


def test_cdna_to_dna_pos_strand_by_transcript_id_3():
    pos = cdna_to_dna("ENST00000380152", 100, 103)
    assert pos
    assert isinstance(pos, list)
    assert pos[0].contig_id == "13"
    assert pos[0].start == 32319109
    assert pos[0].end == 32319112
    assert pos[0].strand == "+"


def test_cdna_to_dna_pos_strand_by_transcript_name():
    pos = cdna_to_dna("BRCA2-201", 1)
    assert pos
    assert isinstance(pos, list)
    assert pos[0].contig_id == "13"
    assert pos[0].start == 32316461
    assert pos[0].end == 32316461
    assert pos[0].strand == "+"


def test_cdna_to_dna_pos_strand_by_transcript_name_2():
    pos = cdna_to_dna("BRCA2-201", 100)
    assert pos
    assert isinstance(pos, list)
    assert pos[0].contig_id == "13"
    assert pos[0].start == 32319109
    assert pos[0].end == 32319109
    assert pos[0].strand == "+"


def test_cdna_to_dna_pos_strand_by_transcript_name_3():
    pos = cdna_to_dna("BRCA2-201", 100, 103)
    assert pos
    assert isinstance(pos, list)
    assert pos[0].contig_id == "13"
    assert pos[0].start == 32319109
    assert pos[0].end == 32319112
    assert pos[0].strand == "+"


def test_cdna_to_dna_neg_strand_by_transcript_id():
    pos = cdna_to_dna("ENST00000256078", 1)
    assert pos
    assert isinstance(pos, list)
    assert pos[0].contig_id == "12"
    assert pos[0].start == 25245384
    assert pos[0].end == 25245384
    assert pos[0].strand == "-"


def test_cdna_to_dna_neg_strand_by_transcript_id_2():
    pos = cdna_to_dna("ENST00000256078", 100)
    assert pos
    assert isinstance(pos, list)
    assert pos[0].contig_id == "12"
    assert pos[0].start == 25245285
    assert pos[0].end == 25245285
    assert pos[0].strand == "-"


def test_cdna_to_dna_neg_strand_by_transcript_id_3():
    pos = cdna_to_dna("ENST00000256078", 100, 103)
    assert pos
    assert isinstance(pos, list)
    assert pos[0].contig_id == "12"
    assert pos[0].start == 25245282
    assert pos[0].end == 25245285
    assert pos[0].strand == "-"


def test_cdna_to_dna_neg_strand_by_transcript_name():
    pos = cdna_to_dna("KRAS-201", 1)
    assert pos
    assert isinstance(pos, list)
    assert pos[0].contig_id == "12"
    assert pos[0].start == 25245384
    assert pos[0].end == 25245384
    assert pos[0].strand == "-"


def test_cdna_to_dna_neg_strand_by_transcript_name_2():
    pos = cdna_to_dna("KRAS-201", 100)
    assert pos
    assert isinstance(pos, list)
    assert pos[0].contig_id == "12"
    assert pos[0].start == 25245285
    assert pos[0].end == 25245285
    assert pos[0].strand == "-"


def test_cdna_to_dna_neg_strand_by_transcript_name_3():
    pos = cdna_to_dna("KRAS-201", 100, 103)
    assert pos
    assert isinstance(pos, list)
    assert pos[0].contig_id == "12"
    assert pos[0].start == 25245282
    assert pos[0].end == 25245285
    assert pos[0].strand == "-"


def test_cdna_to_dna_by_contig_id():
    result = sorted([i.contig_id for i in cdna_to_dna("ENSE00000936617", 1)])
    assert result == ["12"]


def test_cdna_to_dna_by_exon_id():
    result = sorted([i.contig_id for i in cdna_to_dna("ENSE00000936617", 1)])
    assert result == ["12"]


def test_cdna_to_dna_by_gene_id():
    result = sorted([i.contig_id for i in cdna_to_dna("ENSG00000133703", 1)])
    assert result == ["12"]


def test_cdna_to_dna_by_gene_name():
    result = sorted([i.contig_id for i in cdna_to_dna("KRAS", 1)])
    assert result == ["12"]


def test_cdna_to_dna_by_protein_id():
    result = sorted([i.contig_id for i in cdna_to_dna("ENSP00000256078", 1)])
    assert result == ["12"]


def test_cdna_to_exon_pos_strand_by_transcript_id():
    pos = cdna_to_exon("ENST00000380152", 1)
    assert pos
    assert isinstance(pos, list)
    assert pos[0].transcript_id, "ENST00000380152"
    assert pos[0].exon_id, "ENSE00001484009"
    assert pos[0].start == 2
    assert pos[0].end == 2
    assert pos[0].strand == "+"


def test_cdna_to_exon_pos_strand_by_transcript_id_2():
    pos = cdna_to_exon("ENST00000380152", 360)
    assert pos
    assert isinstance(pos, list)
    assert pos[0].transcript_id, "ENST00000380152"
    assert pos[0].exon_id, "ENSE00003659301"
    assert pos[0].start == 4
    assert pos[0].end == 4
    assert pos[0].strand == "+"


def test_cdna_to_exon_pos_strand_by_transcript_id_3():
    pos = cdna_to_exon("ENST00000380152", 360, 380)
    assert pos
    assert isinstance(pos, list)
    assert pos[0].transcript_id, "ENST00000380152"
    assert pos[0].exon_id, "ENSE00003659301"
    assert pos[0].start == 4
    assert pos[0].end == 4
    assert pos[0].strand == "+"


def test_cdna_to_exon_pos_strand_by_transcript_name():
    pos = cdna_to_exon("BRCA2-201", 1)
    assert pos
    assert isinstance(pos, list)
    assert pos[0].transcript_id, "ENST00000380152"
    assert pos[0].exon_id, "ENSE00001484009"
    assert pos[0].start == 2
    assert pos[0].end == 2
    assert pos[0].strand == "+"


def test_cdna_to_exon_pos_strand_by_transcript_name_2():
    pos = cdna_to_exon("BRCA2-201", 360)
    assert pos
    assert isinstance(pos, list)
    assert pos[0].transcript_id, "ENST00000380152"
    assert pos[0].exon_id, "ENSE00003659301"
    assert pos[0].start == 4
    assert pos[0].end == 4
    assert pos[0].strand == "+"


def test_cdna_to_exon_pos_strand_by_transcript_name_3():
    pos = cdna_to_exon("BRCA2-201", 360, 380)
    assert pos
    assert isinstance(pos, list)
    assert pos[0].transcript_id, "ENST00000380152"
    assert pos[0].exon_id, "ENSE00003659301"
    assert pos[0].start == 4
    assert pos[0].end == 4
    assert pos[0].strand == "+"


def test_cdna_to_exon_neg_strand_by_transcript_id():
    pos = cdna_to_exon("ENST00000256078", 1)
    assert pos
    assert isinstance(pos, list)
    assert pos[0].transcript_id, "ENST00000256078"
    assert pos[0].exon_id, "ENSE00000936617"
    assert pos[0].start == 2
    assert pos[0].end == 2
    assert pos[0].strand == "-"


def test_cdna_to_exon_neg_strand_by_transcript_id_2():
    pos = cdna_to_exon("ENST00000256078", 144)
    assert pos
    assert isinstance(pos, list)
    assert pos[0].transcript_id, "ENST00000256078"
    assert pos[0].exon_id, "ENSE00001719809"
    assert pos[0].start == 3
    assert pos[0].end == 3
    assert pos[0].strand == "-"


def test_cdna_to_exon_neg_strand_by_transcript_id_3():
    pos = cdna_to_exon("ENST00000256078", 144, 156)
    assert pos
    assert isinstance(pos, list)
    assert pos[0].transcript_id, "ENST00000256078"
    assert pos[0].exon_id, "ENSE00001719809"
    assert pos[0].start == 3
    assert pos[0].end == 3
    assert pos[0].strand == "-"


def test_cdna_to_exon_neg_strand_by_transcript_name():
    pos = cdna_to_exon("KRAS-201", 1)
    assert pos
    assert isinstance(pos, list)
    assert pos[0].transcript_id, "ENST00000256078"
    assert pos[0].exon_id, "ENSE00000936617"
    assert pos[0].start == 2
    assert pos[0].end == 2
    assert pos[0].strand == "-"


def test_cdna_to_exon_neg_strand_by_transcript_name_2():
    pos = cdna_to_exon("KRAS-201", 144)
    assert pos
    assert isinstance(pos, list)
    assert pos[0].transcript_id, "ENST00000256078"
    assert pos[0].exon_id, "ENSE00001719809"
    assert pos[0].start == 3
    assert pos[0].end == 3
    assert pos[0].strand == "-"


def test_cdna_to_exon_neg_strand_by_transcript_name_3():
    pos = cdna_to_exon("KRAS-201", 144, 156)
    assert pos
    assert isinstance(pos, list)
    assert pos[0].transcript_id, "ENST00000256078"
    assert pos[0].exon_id, "ENSE00001719809"
    assert pos[0].start == 3
    assert pos[0].end == 3
    assert pos[0].strand == "-"


def test_cdna_to_exon_different_exons():
    # positions are on different exons
    with pytest.raises(ValueError):
        _ = cdna_to_exon("ENST00000380152", 300, 1000)


def test_cdna_to_exon_by_contig_id():
    # NOTE: I don't see this query being useful because it will return coordinates for every transcript on chr12.
    # result = sorted([i.exon_id for i in cdna_to_exon("12", 1)])
    pass


def test_cdna_to_exon_by_exon_id():
    result = sorted([i.exon_id for i in cdna_to_exon("ENSE00000936617", 1)])
    # the exon id appears 4 times because the exon belongs to 4 different transcripts
    assert result == ["ENSE00000936617", "ENSE00000936617", "ENSE00000936617", "ENSE00000936617"]


def test_cdna_to_exon_by_gene_id():
    result = sorted([i.exon_id for i in cdna_to_exon("ENSG00000133703", 1)])
    # the exon id appears 4 times because the exon belongs to 4 different transcripts
    assert result == ["ENSE00000936617", "ENSE00000936617", "ENSE00000936617", "ENSE00000936617"]


def test_cdna_to_exon_by_gene_name():
    result = sorted([i.exon_id for i in cdna_to_exon("KRAS", 1)])
    # the exon id appears 4 times because the exon belongs to 4 different transcripts
    assert result == ["ENSE00000936617", "ENSE00000936617", "ENSE00000936617", "ENSE00000936617"]


def test_cdna_to_exon_by_protein_id():
    result = sorted([i.exon_id for i in cdna_to_exon("ENSP00000256078", 1)])
    assert result == ["ENSE00000936617"]


def test_cdna_to_protein_pos_strand_by_transcript_id():
    pos = cdna_to_protein("ENST00000288135", 1)
    assert pos
    assert isinstance(pos, list)
    assert pos[0].protein_id, "ENSP00000288135"
    assert pos[0].start == 1
    assert pos[0].end == 1
    assert pos[0].strand == "+"


def test_cdna_to_protein_pos_strand_by_transcript_id_2():
    pos = cdna_to_protein("ENST00000288135", 213)
    assert pos
    assert isinstance(pos, list)
    assert pos[0].protein_id, "ENSP00000288135"
    assert pos[0].start == 71
    assert pos[0].end == 71
    assert pos[0].strand == "+"


def test_cdna_to_protein_pos_strand_by_transcript_id_3():
    pos = cdna_to_protein("ENST00000288135", 213, 221)
    assert pos
    assert isinstance(pos, list)
    assert pos[0].protein_id, "ENSP00000288135"
    assert pos[0].start == 71
    assert pos[0].end == 74
    assert pos[0].strand == "+"


def test_cdna_to_protein_pos_strand_by_transcript_name():
    pos = cdna_to_protein("KIT-201", 1)
    assert pos
    assert isinstance(pos, list)
    assert pos[0].protein_id, "ENSP00000288135"
    assert pos[0].start == 1
    assert pos[0].end == 1
    assert pos[0].strand == "+"


def test_cdna_to_protein_pos_strand_by_transcript_name_2():
    pos = cdna_to_protein("KIT-201", 213)
    assert pos
    assert isinstance(pos, list)
    assert pos[0].protein_id, "ENSP00000288135"
    assert pos[0].start == 71
    assert pos[0].end == 71
    assert pos[0].strand == "+"


def test_cdna_to_protein_pos_strand_by_transcript_name_3():
    pos = cdna_to_protein("KIT-201", 213, 221)
    assert pos
    assert isinstance(pos, list)
    assert pos[0].protein_id, "ENSP00000288135"
    assert pos[0].start == 71
    assert pos[0].end == 74
    assert pos[0].strand == "+"


def test_cdna_to_protein_neg_strand_by_transcript_id():
    pos = cdna_to_protein("ENST00000256078", 1)
    assert pos
    assert isinstance(pos, list)
    assert pos[0].protein_id, "ENSP00000256078"
    assert pos[0].start == 1
    assert pos[0].end == 1
    assert pos[0].strand == "-"


def test_cdna_to_protein_neg_strand_by_transcript_id_2():
    pos = cdna_to_protein("ENST00000256078", 26)
    assert pos
    assert isinstance(pos, list)
    assert pos[0].protein_id, "ENSP00000256078"
    assert pos[0].start == 9
    assert pos[0].end == 9
    assert pos[0].strand == "-"


def test_cdna_to_protein_neg_strand_by_transcript_id_3():
    pos = cdna_to_protein("ENST00000256078", 26, 56)
    assert pos
    assert isinstance(pos, list)
    assert pos[0].protein_id, "ENSP00000256078"
    assert pos[0].start == 9
    assert pos[0].end == 19
    assert pos[0].strand == "-"


def test_cdna_to_protein_neg_strand_by_transcript_name():
    pos = cdna_to_protein("KRAS-201", 1)
    assert pos
    assert isinstance(pos, list)
    assert pos[0].protein_id, "ENSP00000256078"
    assert pos[0].start == 1
    assert pos[0].end == 1
    assert pos[0].strand == "-"


def test_cdna_to_protein_neg_strand_by_transcript_name_2():
    pos = cdna_to_protein("KRAS-201", 26)
    assert pos
    assert isinstance(pos, list)
    assert pos[0].protein_id, "ENSP00000256078"
    assert pos[0].start == 9
    assert pos[0].end == 9
    assert pos[0].strand == "-"


def test_cdna_to_protein_neg_strand_by_transcript_name_3():
    pos = cdna_to_protein("KRAS-201", 26, 56)
    assert pos
    assert isinstance(pos, list)
    assert pos[0].protein_id, "ENSP00000256078"
    assert pos[0].start == 9
    assert pos[0].end == 19
    assert pos[0].strand == "-"


def test_cdna_to_protein_by_contig_id():
    # NOTE: I don't see this query being useful because it will return coordinates for every transcript on chr12.
    # result = sorted([i.protein_id for i in cdna_to_protein("12", 1)])
    pass


def test_cdna_to_protein_by_exon_id():
    result = sorted([i.protein_id for i in cdna_to_protein("ENSE00000936617", 1)])
    assert result == ["ENSP00000256078", "ENSP00000308495", "ENSP00000451856", "ENSP00000452512"]


def test_cdna_to_protein_by_gene_id():
    result = sorted([i.protein_id for i in cdna_to_protein("ENSG00000133703", 1)])
    assert result == ["ENSP00000256078", "ENSP00000308495", "ENSP00000451856", "ENSP00000452512"]


def test_cdna_to_protein_by_gene_name():
    result = sorted([i.protein_id for i in cdna_to_protein("KRAS", 1)])
    assert result == ["ENSP00000256078", "ENSP00000308495", "ENSP00000451856", "ENSP00000452512"]


def test_cdna_to_protein_by_protein_id():
    result = sorted([i.protein_id for i in cdna_to_protein("ENSP00000256078", 1)])
    assert result == ["ENSP00000256078"]


def test_cdna_to_rna_pos_strand_by_transcript_id():
    pos = cdna_to_rna("ENST00000288135", 1)
    assert pos
    assert isinstance(pos, list)
    assert pos[0].transcript_id, "ENST00000288135"
    assert pos[0].start == 59
    assert pos[0].end == 59
    assert pos[0].strand == "+"


def test_cdna_to_rna_pos_strand_by_transcript_id_2():
    pos = cdna_to_rna("ENST00000288135", 213)
    assert pos
    assert isinstance(pos, list)
    assert pos[0].transcript_id, "ENST00000288135"
    assert pos[0].start == 271
    assert pos[0].end == 271
    assert pos[0].strand == "+"


def test_cdna_to_rna_pos_strand_by_transcript_id_3():
    pos = cdna_to_rna("ENST00000288135", 213, 221)
    assert pos
    assert isinstance(pos, list)
    assert pos[0].transcript_id, "ENST00000288135"
    assert pos[0].start == 271
    assert pos[0].end == 279
    assert pos[0].strand == "+"


def test_cdna_to_rna_pos_strand_by_transcript_name():
    pos = cdna_to_rna("KIT-201", 1)
    assert pos
    assert isinstance(pos, list)
    assert pos[0].transcript_id, "ENST00000288135"
    assert pos[0].start == 59
    assert pos[0].end == 59
    assert pos[0].strand == "+"


def test_cdna_to_rna_pos_strand_by_transcript_name_2():
    pos = cdna_to_rna("KIT-201", 213)
    assert pos
    assert isinstance(pos, list)
    assert pos[0].transcript_id, "ENST00000288135"
    assert pos[0].start == 271
    assert pos[0].end == 271
    assert pos[0].strand == "+"


def test_cdna_to_rna_pos_strand_by_transcript_name_3():
    pos = cdna_to_rna("KIT-201", 213, 221)
    assert pos
    assert isinstance(pos, list)
    assert pos[0].transcript_id, "ENST00000288135"
    assert pos[0].start == 271
    assert pos[0].end == 279
    assert pos[0].strand == "+"


def test_cdna_to_rna_neg_strand_by_transcript_id():
    pos = cdna_to_rna("ENST00000256078", 1)
    assert pos
    assert isinstance(pos, list)
    assert pos[0].transcript_id, "ENST00000256078"
    assert pos[0].start == 191
    assert pos[0].end == 191
    assert pos[0].strand == "-"


def test_cdna_to_rna_neg_strand_by_transcript_id_2():
    pos = cdna_to_rna("ENST00000256078", 87)
    assert pos
    assert isinstance(pos, list)
    assert pos[0].transcript_id, "ENST00000256078"
    assert pos[0].start == 277
    assert pos[0].end == 277
    assert pos[0].strand == "-"


def test_cdna_to_rna_neg_strand_by_transcript_id_3():
    pos = cdna_to_rna("ENST00000256078", 87, 92)
    assert pos
    assert isinstance(pos, list)
    assert pos[0].transcript_id, "ENST00000256078"
    assert pos[0].start == 277
    assert pos[0].end == 282
    assert pos[0].strand == "-"


def test_cdna_to_rna_neg_strand_by_transcript_name():
    pos = cdna_to_rna("KRAS-201", 1)
    assert pos
    assert isinstance(pos, list)
    assert pos[0].transcript_id, "ENST00000256078"
    assert pos[0].start == 191
    assert pos[0].end == 191
    assert pos[0].strand == "-"


def test_cdna_to_rna_neg_strand_by_transcript_name_2():
    pos = cdna_to_rna("KRAS-201", 87)
    assert pos
    assert isinstance(pos, list)
    assert pos[0].transcript_id, "ENST00000256078"
    assert pos[0].start == 277
    assert pos[0].end == 277
    assert pos[0].strand == "-"


def test_cdna_to_rna_neg_strand_by_transcript_name_3():
    pos = cdna_to_rna("KRAS-201", 87, 92)
    assert pos
    assert isinstance(pos, list)
    assert pos[0].transcript_id, "ENST00000256078"
    assert pos[0].start == 277
    assert pos[0].end == 282
    assert pos[0].strand == "-"


def test_cdna_to_rna_by_contig_id():
    # NOTE: I don't see this query being useful because it will return coordinates for every transcript on chr12.
    # result = sorted([i.transcript_id for i in cdna_to_rna("12", 1)])
    pass


def test_cdna_to_rna_by_exon_id():
    result = sorted([i.transcript_id for i in cdna_to_rna("ENSE00000936617", 1)])
    assert result == ["ENST00000256078", "ENST00000311936", "ENST00000556131", "ENST00000557334"]


def test_cdna_to_rna_by_gene_id():
    result = sorted([i.transcript_id for i in cdna_to_rna("ENSG00000133703", 1)])
    assert result == ["ENST00000256078", "ENST00000311936", "ENST00000556131", "ENST00000557334"]


def test_cdna_to_rna_by_gene_name():
    result = sorted([i.transcript_id for i in cdna_to_rna("KRAS", 1)])
    assert result == ["ENST00000256078", "ENST00000311936", "ENST00000556131", "ENST00000557334"]


def test_cdna_to_rna_by_protein_id():
    result = sorted([i.transcript_id for i in cdna_to_rna("ENSP00000256078", 1)])
    assert result == ["ENST00000256078"]


def test_dna_to_cdna():
    pos = dna_to_cdna("5", 1294501)
    assert pos
    assert isinstance(pos, list)
    assert pos[0].transcript_id, "ENST00000310581"
    assert pos[0].start == 385
    assert pos[0].end == 385
    assert pos[0].strand == "-"


def test_dna_to_cdna_2():
    pos = dna_to_cdna("5", 1294497, 1294501)
    assert pos
    assert isinstance(pos, list)
    assert pos[0].transcript_id, "ENST00000310581"
    assert pos[0].start == 385
    assert pos[0].end == 389
    assert pos[0].strand == "-"


def test_dna_to_cdna_by_exon_id():
    result = sorted([i.transcript_id for i in dna_to_cdna("ENSE00000936617", 25245275)])
    assert result == ["ENST00000256078", "ENST00000311936", "ENST00000556131", "ENST00000557334"]


def test_dna_to_cdna_by_gene_id():
    result = sorted([i.transcript_id for i in dna_to_cdna("ENSG00000133703", 25245275)])
    assert result == ["ENST00000256078", "ENST00000311936", "ENST00000556131", "ENST00000557334"]


def test_dna_to_cdna_by_gene_name():
    result = sorted([i.transcript_id for i in dna_to_cdna("KRAS", 25245275)])
    assert result == ["ENST00000256078", "ENST00000311936", "ENST00000556131", "ENST00000557334"]


def test_dna_to_cdna_by_protein_id():
    result = sorted([i.transcript_id for i in dna_to_cdna("ENSP00000256078", 25245275)])
    assert result == ["ENST00000256078", "ENST00000311936", "ENST00000556131", "ENST00000557334"]


def test_dna_to_cdna_by_transcript_id():
    result = sorted([i.transcript_id for i in dna_to_cdna("ENST00000256078", 25245275)])
    assert result == ["ENST00000256078", "ENST00000311936", "ENST00000556131", "ENST00000557334"]


def test_dna_to_cdna_by_transcript_name():
    result = sorted([i.transcript_id for i in dna_to_cdna("KRAS-201", 25245275)])
    assert result == ["ENST00000256078", "ENST00000311936", "ENST00000556131", "ENST00000557334"]


def test_dna_to_dna():
    pos = dna_to_dna("5", 1253147, 1295068)
    assert pos
    assert isinstance(pos, list)
    assert pos[0].contig_id == "5"
    assert pos[0].start == 1253147
    assert pos[0].end == 1295068
    assert pos[0].strand == "+"


def test_dna_to_dna_2():
    pos = dna_to_dna("5", 1254000)
    assert pos
    assert isinstance(pos, list)
    assert pos[0].contig_id == "5"
    assert pos[0].start == 1254000
    assert pos[0].end == 1254000
    assert pos[0].strand == "+"


def test_dna_to_dna_by_exon_id():
    result = sorted([i.contig_id for i in dna_to_dna("ENSE00000936617", 25245275)])
    assert result == ["12"]


def test_dna_to_dna_by_gene_id():
    result = sorted([i.contig_id for i in dna_to_dna("ENSG00000133703", 25245275)])
    assert result == ["12"]


def test_dna_to_dna_by_gene_name():
    result = sorted([i.contig_id for i in dna_to_dna("KRAS", 25245275)])
    assert result == ["12"]


def test_dna_to_dna_by_protein_id():
    result = sorted([i.contig_id for i in dna_to_dna("ENSP00000256078", 25245275)])
    assert result == ["12"]


def test_dna_to_dna_by_transcript_id():
    result = sorted([i.contig_id for i in dna_to_dna("ENST00000256078", 25245275)])
    assert result == ["12"]


def test_dna_to_dna_by_transcript_name():
    result = sorted([i.contig_id for i in dna_to_dna("KRAS-201", 25245275)])
    assert result == ["12"]


def test_dna_to_dna_with_chr():
    pos = dna_to_dna("chr5", 1253147, 1295068)
    assert pos
    assert isinstance(pos, list)
    assert pos[0].contig_id == "5"
    assert pos[0].start == 1253147
    assert pos[0].end == 1295068
    assert pos[0].strand == "+"


def test_dna_to_protein():
    pos = dna_to_protein("5", 1294501)
    assert pos
    assert isinstance(pos, list)
    assert pos[0].protein_id, "ENSP00000309572"
    assert pos[0].start == 129
    assert pos[0].end == 129
    assert pos[0].strand == "-"


def test_dna_to_protein_2():
    pos = dna_to_protein("5", 1294497, 1294501)
    assert pos
    assert isinstance(pos, list)
    assert pos[0].protein_id, "ENSP00000309572"
    assert pos[0].start == 129
    assert pos[0].end == 130
    assert pos[0].strand == "-"


def test_dna_to_protein_by_exon_id():
    result = sorted([i.protein_id for i in dna_to_protein("ENSE00000936617", 25245275)])
    assert result == ["ENSP00000256078", "ENSP00000308495", "ENSP00000451856", "ENSP00000452512"]


def test_dna_to_protein_by_gene_id():
    result = sorted([i.protein_id for i in dna_to_protein("ENSG00000133703", 25245275)])
    assert result == ["ENSP00000256078", "ENSP00000308495", "ENSP00000451856", "ENSP00000452512"]


def test_dna_to_protein_by_gene_name():
    result = sorted([i.protein_id for i in dna_to_protein("KRAS", 25245275)])
    assert result == ["ENSP00000256078", "ENSP00000308495", "ENSP00000451856", "ENSP00000452512"]


def test_dna_to_protein_by_protein_id():
    result = sorted([i.protein_id for i in dna_to_protein("ENSP00000256078", 25245275)])
    assert result == ["ENSP00000256078", "ENSP00000308495", "ENSP00000451856", "ENSP00000452512"]


def test_dna_to_protein_by_transcript_id():
    result = sorted([i.protein_id for i in dna_to_protein("ENST00000256078", 25245275)])
    assert result == ["ENSP00000256078", "ENSP00000308495", "ENSP00000451856", "ENSP00000452512"]


def test_dna_to_protein_by_transcript_name():
    result = sorted([i.protein_id for i in dna_to_protein("KRAS-201", 25245275)])
    assert result == ["ENSP00000256078", "ENSP00000308495", "ENSP00000451856", "ENSP00000452512"]


def test_dna_to_rna():
    pos = dna_to_rna("5", 1294501)
    assert pos
    assert isinstance(pos, list)
    assert pos[0].transcript_id, "ENST00000310581"
    assert pos[0].start == 464
    assert pos[0].end == 464
    assert pos[0].strand == "-"


def test_dna_to_rna_2():
    pos = dna_to_rna("5", 1294497, 1294501)
    assert pos
    assert isinstance(pos, list)
    assert pos[0].transcript_id, "ENST00000310581"
    assert pos[0].start == 464
    assert pos[0].end == 468
    assert pos[0].strand == "-"


def test_dna_to_rna_by_exon_id():
    result = sorted([i.transcript_id for i in dna_to_rna("ENSE00000936617", 25245275)])

    assert result == ["ENST00000256078", "ENST00000311936", "ENST00000556131", "ENST00000557334"]


def test_dna_to_rna_by_gene_id():
    result = sorted([i.transcript_id for i in dna_to_rna("ENSG00000133703", 25245275)])

    assert result == ["ENST00000256078", "ENST00000311936", "ENST00000556131", "ENST00000557334"]


def test_dna_to_rna_by_gene_name():
    result = sorted([i.transcript_id for i in dna_to_rna("KRAS", 25245275)])
    assert result == ["ENST00000256078", "ENST00000311936", "ENST00000556131", "ENST00000557334"]


def test_dna_to_rna_by_protein_id():
    result = sorted([i.transcript_id for i in dna_to_rna("ENSP00000256078", 25245275)])

    assert result == ["ENST00000256078", "ENST00000311936", "ENST00000556131", "ENST00000557334"]


def test_dna_to_rna_by_transcript_id():
    result = sorted([i.transcript_id for i in dna_to_rna("ENST00000256078", 25245275)])

    assert result == ["ENST00000256078", "ENST00000311936", "ENST00000556131", "ENST00000557334"]


def test_dna_to_rna_by_transcript_name():
    result = sorted([i.transcript_id for i in dna_to_rna("KRAS-201", 25245275)])
    assert result == ["ENST00000256078", "ENST00000311936", "ENST00000556131", "ENST00000557334"]


def test_exon_to_cdna_pos_strand_by_exon_id():
    pos = exon_to_cdna("ENSE00003659301")
    assert pos
    assert isinstance(pos, list)
    assert pos[0].transcript_id, "ENST00000380152"
    assert pos[0].start == 317
    assert pos[0].end == 425
    assert pos[0].strand == "+"


def test_exon_to_cdna_neg_strand_by_exon_id():
    pos = exon_to_cdna("ENSE00001719809")
    assert pos
    assert isinstance(pos, list)
    assert pos[0].transcript_id, "ENST00000256078"
    assert pos[0].start == 112
    assert pos[0].end == 290
    assert pos[0].strand == "-"


def test_exon_to_cdna_neg_strand_by_exon_id_2():
    pos = exon_to_cdna("ENSE00000936617")
    assert pos
    assert isinstance(pos, list)
    assert pos[0].transcript_id, "ENST00000256078"
    assert pos[0].start == 1
    assert pos[0].end == 111
    assert pos[0].strand == "-"


def test_exon_to_cdna_neg_strand_by_exon_id_3():
    # 'ENSE00001189804' is not part of the CDNA
    with pytest.raises(ValueError):
        _ = exon_to_cdna("ENSE00001189804")


def test_exon_to_cdna_by_contig_id():
    # NOTE: I don't see this query being useful because it will return coordinates for every exon on chr12.
    # result = sorted([i.transcript_id for i in exon_to_cdna("12")])
    pass


def test_exon_to_cdna_by_gene_id():
    result = sorted([i.transcript_id for i in exon_to_cdna("ENSG00000133703")])
    assert result == [
        "ENST00000256078",
        "ENST00000256078",
        "ENST00000256078",
        "ENST00000256078",
        "ENST00000311936",
        "ENST00000311936",
        "ENST00000311936",
        "ENST00000311936",
        "ENST00000556131",
        "ENST00000556131",
        "ENST00000557334",
        "ENST00000557334",
    ]


def test_exon_to_cdna_by_gene_name():
    result = sorted([i.transcript_id for i in exon_to_cdna("KRAS")])
    assert result == [
        "ENST00000256078",
        "ENST00000256078",
        "ENST00000256078",
        "ENST00000256078",
        "ENST00000311936",
        "ENST00000311936",
        "ENST00000311936",
        "ENST00000311936",
        "ENST00000556131",
        "ENST00000556131",
        "ENST00000557334",
        "ENST00000557334",
    ]


def test_exon_to_cdna_by_protein_id():
    result = sorted([i.transcript_id for i in exon_to_cdna("ENSP00000256078")])
    assert result == [
        "ENST00000256078",
        "ENST00000256078",
        "ENST00000256078",
        "ENST00000256078",
        "ENST00000311936",
        "ENST00000311936",
        "ENST00000311936",
        "ENST00000311936",
        "ENST00000556131",
        "ENST00000557334",
        "ENST00000557334",
    ]


def test_exon_to_cdna_by_transcript_id():
    result = sorted([i.transcript_id for i in exon_to_cdna("ENST00000256078")])
    assert result == [
        "ENST00000256078",
        "ENST00000256078",
        "ENST00000256078",
        "ENST00000256078",
        "ENST00000311936",
        "ENST00000311936",
        "ENST00000311936",
        "ENST00000311936",
        "ENST00000556131",
        "ENST00000557334",
        "ENST00000557334",
    ]


def test_exon_to_cdna_by_transcript_name():
    result = sorted([i.transcript_id for i in exon_to_cdna("KRAS-201")])
    assert result == [
        "ENST00000256078",
        "ENST00000256078",
        "ENST00000256078",
        "ENST00000256078",
        "ENST00000311936",
        "ENST00000311936",
        "ENST00000311936",
        "ENST00000311936",
        "ENST00000556131",
        "ENST00000557334",
        "ENST00000557334",
    ]


def test_exon_to_dna_pos_strand_by_exon_id():
    pos = exon_to_dna("ENSE00003659301")
    assert pos
    assert isinstance(pos, list)
    assert pos[0].contig_id == "13"
    assert pos[0].start == 32325076
    assert pos[0].end == 32325184
    assert pos[0].strand == "+"


def test_exon_to_dna_neg_strand_by_exon_id():
    pos = exon_to_dna("ENSE00001189807")
    assert pos
    assert isinstance(pos, list)
    assert pos[0].contig_id == "12"
    assert pos[0].start == 25215437
    assert pos[0].end == 25215560
    assert pos[0].strand == "+"


def test_exon_to_dna_by_contig_id():
    # NOTE: I don't see this query being useful because it will return coordinates for every exon on chr12.
    # result = sorted([i.contig_id for i in exon_to_dna("12")])
    pass


def test_exon_to_dna_by_gene_id():
    result = sorted([i.contig_id for i in exon_to_dna("ENSG00000133703")])
    # multiple results because they are for different positions on chr12
    assert result == ["12", "12", "12", "12", "12", "12", "12", "12", "12", "12"]


def test_exon_to_dna_by_gene_name():
    result = sorted([i.contig_id for i in exon_to_dna("KRAS")])
    # multiple results because they are for different positions on chr12
    assert result == ["12", "12", "12", "12", "12", "12", "12", "12", "12", "12"]


def test_exon_to_dna_by_protein_id():
    result = sorted([i.contig_id for i in exon_to_dna("ENSP00000256078")])
    # multiple results because they are for different positions on chr12
    assert result == ["12", "12", "12", "12", "12", "12"]


def test_exon_to_dna_by_transcript_id():
    result = sorted([i.contig_id for i in exon_to_dna("ENST00000256078")])
    # multiple results because they are for different positions on chr12
    assert result == ["12", "12", "12", "12", "12", "12"]


def test_exon_to_dna_by_transcript_name():
    result = sorted([i.contig_id for i in exon_to_dna("KRAS-201")])
    # multiple results because they are for different positions on chr12
    assert result == ["12", "12", "12", "12", "12", "12"]


def test_exon_to_exon_pos_strand_by_exon_id():
    pos = exon_to_exon("ENSE00003659301")
    assert pos
    assert isinstance(pos, list)
    assert pos[0].transcript_id, "ENST00000380152"
    assert pos[0].exon_id, "ENSE00003659301"
    assert pos[0].start == 32325076
    assert pos[0].end == 32325184
    assert pos[0].strand == "+"


def test_exon_to_exon_neg_strand_by_exon_id():
    pos = exon_to_exon("ENSE00001189807")
    assert pos
    assert isinstance(pos, list)
    assert pos[0].transcript_id, "ENST00000256078"
    assert pos[0].exon_id, "ENSE00001189807"
    assert pos[0].start == 25215437
    assert pos[0].end == 25215560
    assert pos[0].strand == "-"


def test_exon_to_exon_by_contig_id():
    # NOTE: I don't see this query being useful because it will return coordinates for every exon on chr12.
    # result = sorted([i.exon_id for i in exon_to_exon("12")])
    pass


def test_exon_to_exon_by_gene_id():
    result = sorted([i.exon_id for i in exon_to_exon("ENSG00000133703")])
    assert result == [
        "ENSE00000000028",
        "ENSE00000936617",
        "ENSE00000936617",
        "ENSE00000936617",
        "ENSE00000936617",
        "ENSE00001189804",
        "ENSE00001189807",
        "ENSE00001644818",
        "ENSE00001644818",
        "ENSE00001719809",
        "ENSE00001719809",
        "ENSE00002446502",
        "ENSE00002456976",
        "ENSE00002464674",
        "ENSE00002477035",
        "ENSE00002478081",
        "ENSE00002530521",
    ]


def test_exon_to_exon_by_gene_name():
    result = sorted([i.exon_id for i in exon_to_exon("KRAS")])
    assert result == [
        "ENSE00000000028",
        "ENSE00000936617",
        "ENSE00000936617",
        "ENSE00000936617",
        "ENSE00000936617",
        "ENSE00001189804",
        "ENSE00001189807",
        "ENSE00001644818",
        "ENSE00001644818",
        "ENSE00001719809",
        "ENSE00001719809",
        "ENSE00002446502",
        "ENSE00002456976",
        "ENSE00002464674",
        "ENSE00002477035",
        "ENSE00002478081",
        "ENSE00002530521",
    ]


def test_exon_to_exon_by_protein_id():
    result = sorted([i.exon_id for i in exon_to_exon("ENSP00000256078")])
    assert result == [
        "ENSE00000000028",
        "ENSE00000936617",
        "ENSE00000936617",
        "ENSE00000936617",
        "ENSE00000936617",
        "ENSE00001189807",
        "ENSE00001644818",
        "ENSE00001644818",
        "ENSE00001719809",
        "ENSE00001719809",
        "ENSE00002477035",
    ]


def test_exon_to_exon_by_transcript_id():
    result = sorted([i.exon_id for i in exon_to_exon("ENST00000256078")])
    assert result == [
        "ENSE00000000028",
        "ENSE00000936617",
        "ENSE00000936617",
        "ENSE00000936617",
        "ENSE00000936617",
        "ENSE00001189807",
        "ENSE00001644818",
        "ENSE00001644818",
        "ENSE00001719809",
        "ENSE00001719809",
        "ENSE00002477035",
    ]


def test_exon_to_exon_by_transcript_name():
    result = sorted([i.exon_id for i in exon_to_exon("KRAS-201")])
    assert result == [
        "ENSE00000000028",
        "ENSE00000936617",
        "ENSE00000936617",
        "ENSE00000936617",
        "ENSE00000936617",
        "ENSE00001189807",
        "ENSE00001644818",
        "ENSE00001644818",
        "ENSE00001719809",
        "ENSE00001719809",
        "ENSE00002477035",
    ]


def test_exon_to_protein_pos_strand_by_exon_id():
    pos = exon_to_protein("ENSE00003659301")
    assert pos
    assert isinstance(pos, list)
    assert pos[0].protein_id, "ENSP00000369497"
    assert pos[0].start == 106
    assert pos[0].end == 142
    assert pos[0].strand == "+"


def test_exon_to_protein_neg_strand_by_exon_id():
    pos = exon_to_protein("ENSE00001719809")
    assert pos
    assert isinstance(pos, list)
    assert pos[0].protein_id, "ENSP00000256078"
    assert pos[0].start == 38
    assert pos[0].end == 97
    assert pos[0].strand == "-"


def test_exon_to_protein_by_contig_id():
    # NOTE: I don't see this query being useful because it will return coordinates for every exon on chr12.
    # result = sorted([i.protein_id for i in exon_to_protein("12")])
    pass


def test_exon_to_protein_by_gene_id():
    result = sorted([i.protein_id for i in exon_to_protein("ENSG00000133703")])
    assert result == [
        "ENSP00000256078",
        "ENSP00000256078",
        "ENSP00000256078",
        "ENSP00000256078",
        "ENSP00000308495",
        "ENSP00000308495",
        "ENSP00000308495",
        "ENSP00000308495",
        "ENSP00000451856",
        "ENSP00000451856",
        "ENSP00000452512",
        "ENSP00000452512",
    ]


def test_exon_to_protein_by_gene_name():
    result = sorted([i.protein_id for i in exon_to_protein("KRAS")])
    assert result == [
        "ENSP00000256078",
        "ENSP00000256078",
        "ENSP00000256078",
        "ENSP00000256078",
        "ENSP00000308495",
        "ENSP00000308495",
        "ENSP00000308495",
        "ENSP00000308495",
        "ENSP00000451856",
        "ENSP00000451856",
        "ENSP00000452512",
        "ENSP00000452512",
    ]


def test_exon_to_protein_by_protein_id():
    result = sorted([i.protein_id for i in exon_to_protein("ENSP00000256078")])
    assert result == [
        "ENSP00000256078",
        "ENSP00000256078",
        "ENSP00000256078",
        "ENSP00000256078",
        "ENSP00000308495",
        "ENSP00000308495",
        "ENSP00000308495",
        "ENSP00000308495",
        "ENSP00000451856",
        "ENSP00000452512",
        "ENSP00000452512",
    ]


def test_exon_to_protein_by_transcript_id():
    result = sorted([i.protein_id for i in exon_to_protein("ENST00000256078")])
    assert result == [
        "ENSP00000256078",
        "ENSP00000256078",
        "ENSP00000256078",
        "ENSP00000256078",
        "ENSP00000308495",
        "ENSP00000308495",
        "ENSP00000308495",
        "ENSP00000308495",
        "ENSP00000451856",
        "ENSP00000452512",
        "ENSP00000452512",
    ]


def test_exon_to_protein_by_transcript_name():
    result = sorted([i.protein_id for i in exon_to_protein("KRAS-201")])
    assert result == [
        "ENSP00000256078",
        "ENSP00000256078",
        "ENSP00000256078",
        "ENSP00000256078",
        "ENSP00000308495",
        "ENSP00000308495",
        "ENSP00000308495",
        "ENSP00000308495",
        "ENSP00000451856",
        "ENSP00000452512",
        "ENSP00000452512",
    ]


def test_exon_to_rna_pos_strand_by_exon_id():
    pos = exon_to_rna("ENSE00003659301")
    assert pos
    assert isinstance(pos, list)
    assert pos[0].transcript_id, "ENST00000380152"
    assert pos[0].start == 550
    assert pos[0].end == 658
    assert pos[0].strand == "+"


def test_exon_to_rna_neg_strand_by_exon_id():
    pos = exon_to_rna("ENSE00001719809")
    assert pos
    assert isinstance(pos, list)
    assert pos[0].transcript_id, "ENST00000256078"
    assert pos[0].start == 302
    assert pos[0].end == 480
    assert pos[0].strand == "-"


def test_exon_to_rna_by_contig_id():
    # NOTE: I don't see this query being useful because it will return coordinates for every exon on chr12.
    # result = sorted([i.transcript_id for i in exon_to_rna("12")])
    pass


def test_exon_to_rna_by_gene_id():
    result = sorted([i.transcript_id for i in exon_to_rna("ENSG00000133703")])
    assert result == [
        "ENST00000256078",
        "ENST00000256078",
        "ENST00000256078",
        "ENST00000256078",
        "ENST00000256078",
        "ENST00000256078",
        "ENST00000311936",
        "ENST00000311936",
        "ENST00000311936",
        "ENST00000311936",
        "ENST00000311936",
        "ENST00000556131",
        "ENST00000556131",
        "ENST00000556131",
        "ENST00000557334",
        "ENST00000557334",
        "ENST00000557334",
    ]


def test_exon_to_rna_by_gene_name():
    result = sorted([i.transcript_id for i in exon_to_rna("KRAS")])
    assert result == [
        "ENST00000256078",
        "ENST00000256078",
        "ENST00000256078",
        "ENST00000256078",
        "ENST00000256078",
        "ENST00000256078",
        "ENST00000311936",
        "ENST00000311936",
        "ENST00000311936",
        "ENST00000311936",
        "ENST00000311936",
        "ENST00000556131",
        "ENST00000556131",
        "ENST00000556131",
        "ENST00000557334",
        "ENST00000557334",
        "ENST00000557334",
    ]


def test_exon_to_rna_by_protein_id():
    result = sorted([i.transcript_id for i in exon_to_rna("ENSP00000256078")])
    assert result == [
        "ENST00000256078",
        "ENST00000256078",
        "ENST00000256078",
        "ENST00000256078",
        "ENST00000256078",
        "ENST00000256078",
        "ENST00000311936",
        "ENST00000311936",
        "ENST00000311936",
        "ENST00000556131",
        "ENST00000557334",
    ]


def test_exon_to_rna_by_transcript_id():
    result = sorted([i.transcript_id for i in exon_to_rna("ENST00000256078")])
    assert result == [
        "ENST00000256078",
        "ENST00000256078",
        "ENST00000256078",
        "ENST00000256078",
        "ENST00000256078",
        "ENST00000256078",
        "ENST00000311936",
        "ENST00000311936",
        "ENST00000311936",
        "ENST00000556131",
        "ENST00000557334",
    ]


def test_exon_to_rna_by_transcript_name():
    result = sorted([i.transcript_id for i in exon_to_rna("KRAS-201")])
    assert result == [
        "ENST00000256078",
        "ENST00000256078",
        "ENST00000256078",
        "ENST00000256078",
        "ENST00000256078",
        "ENST00000256078",
        "ENST00000311936",
        "ENST00000311936",
        "ENST00000311936",
        "ENST00000556131",
        "ENST00000557334",
    ]


def test_protein_to_cdna_pos_strand_by_protein_id():
    pos = protein_to_cdna("ENSP00000288135", 1)
    assert pos
    assert isinstance(pos, list)
    assert pos[0].transcript_id, "ENST00000288135"
    assert pos[0].start == 1
    assert pos[0].end == 3
    assert pos[0].strand == "+"


def test_protein_to_cdna_pos_strand_by_protein_id_2():
    pos = protein_to_cdna("ENSP00000288135", 12)
    assert pos
    assert isinstance(pos, list)
    assert pos[0].transcript_id, "ENST00000288135"
    assert pos[0].start == 34
    assert pos[0].end == 36
    assert pos[0].strand == "+"


def test_protein_to_cdna_pos_strand_by_protein_id_3():
    pos = protein_to_cdna("ENSP00000288135", 12, 15)
    assert pos
    assert isinstance(pos, list)
    assert pos[0].transcript_id, "ENST00000288135"
    assert pos[0].start == 34
    assert pos[0].end == 45
    assert pos[0].strand == "+"


def test_protein_to_cdna_neg_strand_by_protein_id():
    pos = protein_to_cdna("ENSP00000308495", 1)
    assert pos
    assert isinstance(pos, list)
    assert pos[0].transcript_id, "ENST00000311936"
    assert pos[0].start == 1
    assert pos[0].end == 3
    assert pos[0].strand == "-"


def test_protein_to_cdna_neg_strand_by_protein_id_2():
    pos = protein_to_cdna("ENSP00000308495", 123)
    assert pos
    assert isinstance(pos, list)
    assert pos[0].transcript_id, "ENST00000311936"
    assert pos[0].start == 367
    assert pos[0].end == 369
    assert pos[0].strand == "-"


def test_protein_to_cdna_neg_strand_by_protein_id_3():
    pos = protein_to_cdna("ENSP00000308495", 123, 124)
    assert pos
    assert isinstance(pos, list)
    assert pos[0].transcript_id, "ENST00000311936"
    assert pos[0].start == 367
    assert pos[0].end == 372
    assert pos[0].strand == "-"


def test_protein_to_cdna_by_contig_id():
    # NOTE: I don't see this query being useful because it will return coordinates for every protein on chr12.
    # result = sorted([i.transcript_id for i in protein_to_cdna("12", 1)])
    pass


def test_protein_to_cdna_by_exon_id():
    result = sorted([i.transcript_id for i in protein_to_cdna("ENSE00000936617", 1)])
    assert result == ["ENST00000256078", "ENST00000311936", "ENST00000556131", "ENST00000557334"]


def test_protein_to_cdna_by_gene_id():
    result = sorted([i.transcript_id for i in protein_to_cdna("ENSG00000133703", 1)])
    assert result == ["ENST00000256078", "ENST00000311936", "ENST00000556131", "ENST00000557334"]


def test_protein_to_cdna_by_gene_name():
    result = sorted([i.transcript_id for i in protein_to_cdna("KRAS", 1)])
    assert result == ["ENST00000256078", "ENST00000311936", "ENST00000556131", "ENST00000557334"]


def test_protein_to_cdna_by_transcript_id():
    result = sorted([i.transcript_id for i in protein_to_cdna("ENST00000256078", 1)])
    assert result == ["ENST00000256078"]


def test_protein_to_cdna_by_transcript_name():
    result = sorted([i.transcript_id for i in protein_to_cdna("KRAS-201", 1)])
    assert result == ["ENST00000256078"]


def test_protein_to_dna_pos_strand_by_protein_id():
    pos = protein_to_dna("ENSP00000288135", 1)
    assert pos
    assert isinstance(pos, list)
    assert pos[0].contig_id == "4"
    assert pos[0].start == 54658015
    assert pos[0].end == 54658017
    assert pos[0].strand == "+"


def test_protein_to_dna_pos_strand_by_protein_id_2():
    pos = protein_to_dna("ENSP00000288135", 71)
    assert pos
    assert isinstance(pos, list)
    assert pos[0].contig_id == "4"
    assert pos[0].start == 54695655
    assert pos[0].end == 54695657
    assert pos[0].strand == "+"


def test_protein_to_dna_pos_strand_by_protein_id_3():
    pos = protein_to_dna("ENSP00000288135", 71, 74)
    assert pos
    assert isinstance(pos, list)
    assert pos[0].contig_id == "4"
    assert pos[0].start == 54695655
    assert pos[0].end == 54695666
    assert pos[0].strand == "+"


def test_protein_to_dna_neg_strand_by_protein_id():
    pos = protein_to_dna("ENSP00000256078", 1)
    assert pos
    assert isinstance(pos, list)
    assert pos[0].contig_id == "12"
    assert pos[0].start == 25245382
    assert pos[0].end == 25245384
    assert pos[0].strand == "-"


def test_protein_to_dna_neg_strand_by_protein_id_2():
    pos = protein_to_dna("ENSP00000256078", 6)
    assert pos
    assert isinstance(pos, list)
    assert pos[0].contig_id == "12"
    assert pos[0].start == 25245367
    assert pos[0].end == 25245369
    assert pos[0].strand == "-"


def test_protein_to_dna_neg_strand_by_protein_id_3():
    pos = protein_to_dna("ENSP00000256078", 6, 9)
    assert pos
    assert isinstance(pos, list)
    assert pos[0].contig_id == "12"
    assert pos[0].start == 25245358
    assert pos[0].end == 25245369
    assert pos[0].strand == "-"


def test_protein_to_dna_by_contig_id():
    # NOTE: I don't see this query being useful because it will return coordinates for every protein on chr12.
    # result = sorted([i.contig_id for i in protein_to_dna("12", 1)])
    pass


def test_protein_to_dna_by_exon_id():
    result = sorted([i.contig_id for i in protein_to_dna("ENSE00000936617", 1)])
    assert result == ["12"]


def test_protein_to_dna_by_gene_id():
    result = sorted([i.contig_id for i in protein_to_dna("ENSG00000133703", 1)])
    assert result == ["12"]


def test_protein_to_dna_by_gene_name():
    result = sorted([i.contig_id for i in protein_to_dna("KRAS", 1)])
    assert result == ["12"]


def test_protein_to_dna_by_transcript_id():
    result = sorted([i.contig_id for i in protein_to_dna("ENST00000256078", 1)])
    assert result == ["12"]


def test_protein_to_dna_by_transcript_name():
    result = sorted([i.contig_id for i in protein_to_dna("KRAS-201", 1)])
    assert result == ["12"]


def test_protein_to_exon_pos_strand_by_protein_id():
    pos = protein_to_exon("ENSP00000288135", 1)
    assert pos
    assert isinstance(pos, list)
    assert pos[0].transcript_id, "ENST00000288135"
    assert pos[0].exon_id, "ENSE00000000233"
    assert pos[0].start == 1
    assert pos[0].end == 1
    assert pos[0].strand == "+"


def test_protein_to_exon_pos_strand_by_protein_id_2():
    pos = protein_to_exon("ENSP00000288135", 4)
    assert pos
    assert isinstance(pos, list)
    assert pos[0].transcript_id, "ENST00000288135"
    assert pos[0].exon_id, "ENSE00000000233"
    assert pos[0].start == 1
    assert pos[0].end == 1
    assert pos[0].strand == "+"


def test_protein_to_exon_pos_strand_by_protein_id_3():
    pos = protein_to_exon("ENSP00000288135", 4, 6)
    assert pos
    assert isinstance(pos, list)
    assert pos[0].transcript_id, "ENST00000288135"
    assert pos[0].exon_id, "ENSE00000000233"
    assert pos[0].start == 1
    assert pos[0].end == 1
    assert pos[0].strand == "+"


def test_protein_to_exon_neg_strand_by_protein_id():
    pos = protein_to_exon("ENSP00000308495", 1)
    assert pos
    assert isinstance(pos, list)
    assert pos[0].transcript_id, "ENST00000311936"
    assert pos[0].exon_id, "ENSE00000936617"
    assert pos[0].start == 2
    assert pos[0].end == 2
    assert pos[0].strand == "-"


def test_protein_to_exon_neg_strand_by_protein_id_2():
    pos = protein_to_exon("ENSP00000308495", 4)
    assert pos
    assert isinstance(pos, list)
    assert pos[0].transcript_id, "ENST00000311936"
    assert pos[0].exon_id, "ENSE00000936617"
    assert pos[0].start == 2
    assert pos[0].end == 2
    assert pos[0].strand == "-"


def test_protein_to_exon_neg_strand_by_protein_id_3():
    pos = protein_to_exon("ENSP00000308495", 4, 6)
    assert pos
    assert isinstance(pos, list)
    assert pos[0].transcript_id, "ENST00000311936"
    assert pos[0].exon_id, "ENSE00000936617"
    assert pos[0].start == 2
    assert pos[0].end == 2
    assert pos[0].strand == "-"


def test_protein_to_exon_by_contig_id():
    # NOTE: I don't see this query being useful because it will return coordinates for every protein on chr12.
    # result = sorted([i.exon_id for i in protein_to_exon("12", 1)])
    pass


def test_protein_to_exon_by_exon_id():
    result = sorted([i.exon_id for i in protein_to_exon("ENSE00000936617", 1)])
    assert result == ["ENSE00000936617", "ENSE00000936617", "ENSE00000936617", "ENSE00000936617"]


def test_protein_to_exon_by_gene_id():
    result = sorted([i.exon_id for i in protein_to_exon("ENSG00000133703", 1)])
    assert result == ["ENSE00000936617", "ENSE00000936617", "ENSE00000936617", "ENSE00000936617"]


def test_protein_to_exon_by_gene_name():
    result = sorted([i.exon_id for i in protein_to_exon("KRAS", 1)])
    assert result == ["ENSE00000936617", "ENSE00000936617", "ENSE00000936617", "ENSE00000936617"]


def test_protein_to_exon_by_transcript_id():
    result = sorted([i.exon_id for i in protein_to_exon("ENST00000256078", 1)])
    assert result == ["ENSE00000936617"]


def test_protein_to_exon_by_transcript_name():
    result = sorted([i.exon_id for i in protein_to_exon("KRAS-201", 1)])
    assert result == ["ENSE00000936617"]


def test_protein_to_protein_pos_strand_by_protein_id():
    pos = protein_to_protein("ENSP00000288135", 1)
    assert pos
    assert isinstance(pos, list)
    assert pos[0].protein_id, "ENSP00000288135"
    assert pos[0].start == 1
    assert pos[0].end == 1
    assert pos[0].strand == "+"


def test_protein_to_protein_pos_strand_by_protein_id_2():
    pos = protein_to_protein("ENSP00000288135", 12)
    assert pos
    assert isinstance(pos, list)
    assert pos[0].protein_id, "ENSP00000288135"
    assert pos[0].start == 12
    assert pos[0].end == 12
    assert pos[0].strand == "+"


def test_protein_to_protein_pos_strand_by_protein_id_3():
    pos = protein_to_protein("ENSP00000288135", 12, 21)
    assert pos
    assert isinstance(pos, list)
    assert pos[0].protein_id, "ENSP00000288135"
    assert pos[0].start == 12
    assert pos[0].end == 21
    assert pos[0].strand == "+"


def test_protein_to_protein_neg_strand_by_protein_id():
    pos = protein_to_protein("ENSP00000260947", 1)
    assert pos
    assert isinstance(pos, list)
    assert pos[0].protein_id, "ENSP00000260947"
    assert pos[0].start == 1
    assert pos[0].end == 1
    assert pos[0].strand == "-"


def test_protein_to_protein_neg_strand_by_protein_id_2():
    pos = protein_to_protein("ENSP00000260947", 71)
    assert pos
    assert isinstance(pos, list)
    assert pos[0].protein_id, "ENSP00000260947"
    assert pos[0].start == 71
    assert pos[0].end == 71
    assert pos[0].strand == "-"


def test_protein_to_protein_neg_strand_by_protein_id_3():
    pos = protein_to_protein("ENSP00000260947", 71, 74)
    assert pos
    assert isinstance(pos, list)
    assert pos[0].protein_id, "ENSP00000260947"
    assert pos[0].start == 71
    assert pos[0].end == 74
    assert pos[0].strand == "-"


def test_protein_to_protein_by_contig_id():
    # NOTE: I don't see this query being useful because it will return coordinates for every protein on chr12.
    # result = sorted([i.protein_id for i in protein_to_protein("12", 1)])
    pass


def test_protein_to_protein_by_exon_id():
    result = sorted([i.protein_id for i in protein_to_protein("ENSE00000936617", 1)])
    assert result == ["ENSP00000256078", "ENSP00000308495", "ENSP00000451856", "ENSP00000452512"]


def test_protein_to_protein_by_gene_id():
    result = sorted([i.protein_id for i in protein_to_protein("ENSG00000133703", 1)])
    assert result == ["ENSP00000256078", "ENSP00000308495", "ENSP00000451856", "ENSP00000452512"]


def test_protein_to_protein_by_gene_name():
    result = sorted([i.protein_id for i in protein_to_protein("KRAS", 1)])
    assert result == ["ENSP00000256078", "ENSP00000308495", "ENSP00000451856", "ENSP00000452512"]


def test_protein_to_protein_by_transcript_id():
    result = sorted([i.protein_id for i in protein_to_protein("ENST00000256078", 1)])
    assert result == ["ENSP00000256078"]


def test_protein_to_protein_by_transcript_name():
    result = sorted([i.protein_id for i in protein_to_protein("KRAS-201", 1)])
    assert result == ["ENSP00000256078"]


def test_protein_to_rna_pos_strand_by_protein_id():
    pos = protein_to_rna("ENSP00000288135", 1)
    assert pos
    assert isinstance(pos, list)
    assert pos[0].transcript_id, "ENST00000288135"
    assert pos[0].start == 59
    assert pos[0].end == 61
    assert pos[0].strand == "+"


def test_protein_to_rna_pos_strand_by_protein_id_2():
    pos = protein_to_rna("ENSP00000288135", 12)
    assert pos
    assert isinstance(pos, list)
    assert pos[0].transcript_id, "ENST00000288135"
    assert pos[0].start == 92
    assert pos[0].end == 94
    assert pos[0].strand == "+"


def test_protein_to_rna_pos_strand_by_protein_id_3():
    pos = protein_to_rna("ENSP00000288135", 12, 21)
    assert pos
    assert isinstance(pos, list)
    assert pos[0].transcript_id, "ENST00000288135"
    assert pos[0].start == 92
    assert pos[0].end == 121
    assert pos[0].strand == "+"


def test_protein_to_rna_neg_strand_by_protein_id():
    pos = protein_to_rna("ENSP00000260947", 1)
    assert pos
    assert isinstance(pos, list)
    assert pos[0].transcript_id, "ENST00000260947"
    assert pos[0].start == 115
    assert pos[0].end == 117
    assert pos[0].strand == "-"


def test_protein_to_rna_neg_strand_by_protein_id_2():
    pos = protein_to_rna("ENSP00000260947", 71)
    assert pos
    assert isinstance(pos, list)
    assert pos[0].transcript_id, "ENST00000260947"
    assert pos[0].start == 325
    assert pos[0].end == 327
    assert pos[0].strand == "-"


def test_protein_to_rna_neg_strand_by_protein_id_3():
    pos = protein_to_rna("ENSP00000260947", 71, 74)
    assert pos
    assert isinstance(pos, list)
    assert pos[0].transcript_id, "ENST00000260947"
    assert pos[0].start == 325
    assert pos[0].end == 336
    assert pos[0].strand == "-"


def test_protein_to_rna_by_contig_id():
    # NOTE: I don't see this query being useful because it will return coordinates for every protein on chr12.
    # result = sorted([i.transcript_id for i in protein_to_rna("12", 1)])
    pass


def test_protein_to_rna_by_exon_id():
    result = sorted([i.transcript_id for i in protein_to_rna("ENSE00000936617", 1)])
    assert result == ["ENST00000256078", "ENST00000311936", "ENST00000556131", "ENST00000557334"]


def test_protein_to_rna_by_gene_id():
    result = sorted([i.transcript_id for i in protein_to_rna("ENSG00000133703", 1)])
    assert result == ["ENST00000256078", "ENST00000311936", "ENST00000556131", "ENST00000557334"]


def test_protein_to_rna_by_gene_name():
    result = sorted([i.transcript_id for i in protein_to_rna("KRAS", 1)])
    assert result == ["ENST00000256078", "ENST00000311936", "ENST00000556131", "ENST00000557334"]


def test_protein_to_rna_by_transcript_id():
    result = sorted([i.transcript_id for i in protein_to_rna("ENST00000256078", 1)])
    assert result == ["ENST00000256078"]


def test_protein_to_rna_by_transcript_name():
    result = sorted([i.transcript_id for i in protein_to_rna("KRAS-201", 1)])
    assert result == ["ENST00000256078"]


def test_rna_to_cdna_pos_strand_by_transcript_id():
    pos = rna_to_cdna("ENST00000288135", 98)
    assert pos
    assert isinstance(pos, list)
    assert pos[0].transcript_id, "ENST00000288135"
    assert pos[0].start == 40
    assert pos[0].end == 40
    assert pos[0].strand == "+"


def test_rna_to_cdna_pos_strand_by_transcript_id_2():
    pos = rna_to_cdna("ENST00000288135", 100)
    assert pos
    assert isinstance(pos, list)
    assert pos[0].transcript_id, "ENST00000288135"
    assert pos[0].start == 42
    assert pos[0].end == 42
    assert pos[0].strand == "+"


def test_rna_to_cdna_pos_strand_by_transcript_id_3():
    pos = rna_to_cdna("ENST00000288135", 98, 100)
    assert pos
    assert isinstance(pos, list)
    assert pos[0].transcript_id, "ENST00000288135"
    assert pos[0].start == 40
    assert pos[0].end == 42
    assert pos[0].strand == "+"


def test_rna_to_cdna_pos_strand_by_transcript_name():
    pos = rna_to_cdna("KIT-201", 98)
    assert pos
    assert isinstance(pos, list)
    assert pos[0].transcript_id, "ENST00000288135"
    assert pos[0].start == 40
    assert pos[0].end == 40
    assert pos[0].strand == "+"


def test_rna_to_cdna_pos_strand_by_transcript_name_2():
    pos = rna_to_cdna("KIT-201", 100)
    assert pos
    assert isinstance(pos, list)
    assert pos[0].transcript_id, "ENST00000288135"
    assert pos[0].start == 42
    assert pos[0].end == 42
    assert pos[0].strand == "+"


def test_rna_to_cdna_pos_strand_by_transcript_name_3():
    pos = rna_to_cdna("KIT-201", 98, 100)
    assert pos
    assert isinstance(pos, list)
    assert pos[0].transcript_id, "ENST00000288135"
    assert pos[0].start == 40
    assert pos[0].end == 42
    assert pos[0].strand == "+"


def test_rna_to_cdna_neg_strand_by_transcript_id():
    pos = rna_to_cdna("ENST00000256078", 191)
    assert pos
    assert isinstance(pos, list)
    assert pos[0].transcript_id, "ENST00000256078"
    assert pos[0].start == 1
    assert pos[0].end == 1
    assert pos[0].strand == "-"


def test_rna_to_cdna_neg_strand_by_transcript_id_2():
    pos = rna_to_cdna("ENST00000256078", 301)
    assert pos
    assert isinstance(pos, list)
    assert pos[0].transcript_id, "ENST00000256078"
    assert pos[0].start == 111
    assert pos[0].end == 111
    assert pos[0].strand == "-"


def test_rna_to_cdna_neg_strand_by_transcript_id_3():
    pos = rna_to_cdna("ENST00000256078", 301, 323)
    assert pos
    assert isinstance(pos, list)
    assert pos[0].transcript_id, "ENST00000256078"
    assert pos[0].start == 111
    assert pos[0].end == 133
    assert pos[0].strand == "-"


def test_rna_to_cdna_neg_strand_by_transcript_name():
    pos = rna_to_cdna("KRAS-201", 191)
    assert pos
    assert isinstance(pos, list)
    assert pos[0].transcript_id, "ENST00000256078"
    assert pos[0].start == 1
    assert pos[0].end == 1
    assert pos[0].strand == "-"


def test_rna_to_cdna_neg_strand_by_transcript_name_2():
    pos = rna_to_cdna("KRAS-201", 301)
    assert pos
    assert isinstance(pos, list)
    assert pos[0].transcript_id, "ENST00000256078"
    assert pos[0].start == 111
    assert pos[0].end == 111
    assert pos[0].strand == "-"


def test_rna_to_cdna_neg_strand_by_transcript_name_3():
    pos = rna_to_cdna("KRAS-201", 301, 323)
    assert pos
    assert isinstance(pos, list)
    assert pos[0].transcript_id, "ENST00000256078"
    assert pos[0].start == 111
    assert pos[0].end == 133
    assert pos[0].strand == "-"


def test_rna_to_cdna_by_contig_id():
    # NOTE: I don't see this query being useful because it will return coordinates for every transcript on chr12.
    # result = sorted([i.transcript_id for i in rna_to_cdna("12", 301)])
    pass


def test_rna_to_cdna_by_exon_id():
    result = sorted([i.transcript_id for i in rna_to_cdna("ENSE00000936617", 301)])
    assert result == ["ENST00000256078", "ENST00000311936", "ENST00000556131", "ENST00000557334"]


def test_rna_to_cdna_by_gene_id():
    result = sorted([i.transcript_id for i in rna_to_cdna("ENSG00000133703", 301)])
    assert result == ["ENST00000256078", "ENST00000311936", "ENST00000556131", "ENST00000557334"]


def test_rna_to_cdna_by_gene_name():
    result = sorted([i.transcript_id for i in rna_to_cdna("KRAS", 301)])
    assert result == ["ENST00000256078", "ENST00000311936", "ENST00000556131", "ENST00000557334"]


def test_rna_to_cdna_by_protein_id():
    result = sorted([i.transcript_id for i in rna_to_cdna("ENSP00000256078", 301)])
    assert result == ["ENST00000256078"]


def test_rna_to_dna_pos_strand_by_transcript_id():
    pos = rna_to_dna("ENST00000341165", 1)
    assert pos
    assert isinstance(pos, list)
    assert pos[0].contig_id == "17"
    assert pos[0].start == 43170481
    assert pos[0].end == 43170481
    assert pos[0].strand == "+"


def test_rna_to_dna_pos_strand_by_transcript_id_2():
    pos = rna_to_dna("ENST00000341165", 731)
    assert pos
    assert isinstance(pos, list)
    assert pos[0].contig_id == "17"
    assert pos[0].start == 43189698
    assert pos[0].end == 43189698
    assert pos[0].strand == "+"


def test_rna_to_dna_pos_strand_by_transcript_id_3():
    pos = rna_to_dna("ENST00000341165", 731, 817)
    assert pos
    assert isinstance(pos, list)
    assert pos[0].contig_id == "17"
    assert pos[0].start == 43189698
    assert pos[0].end == 43189784
    assert pos[0].strand == "+"


def test_rna_to_dna_pos_strand_by_transcript_name():
    pos = rna_to_dna("NBR1-201", 1)
    assert pos
    assert isinstance(pos, list)
    assert pos[0].contig_id == "17"
    assert pos[0].start == 43170481
    assert pos[0].end == 43170481
    assert pos[0].strand == "+"


def test_rna_to_dna_pos_strand_by_transcript_name_2():
    pos = rna_to_dna("NBR1-201", 731)
    assert pos
    assert isinstance(pos, list)
    assert pos[0].contig_id == "17"
    assert pos[0].start == 43189698
    assert pos[0].end == 43189698
    assert pos[0].strand == "+"


def test_rna_to_dna_pos_strand_by_transcript_name_3():
    pos = rna_to_dna("NBR1-201", 731, 817)
    assert pos
    assert isinstance(pos, list)
    assert pos[0].contig_id == "17"
    assert pos[0].start == 43189698
    assert pos[0].end == 43189784
    assert pos[0].strand == "+"


def test_rna_to_dna_neg_strand_by_transcript_id():
    pos = rna_to_dna("ENST00000256078", 1)
    assert pos
    assert isinstance(pos, list)
    assert pos[0].contig_id == "12"
    assert pos[0].start == 25250929
    assert pos[0].end == 25250929
    assert pos[0].strand == "-"


def test_rna_to_dna_neg_strand_by_transcript_id_2():
    pos = rna_to_dna("ENST00000256078", 466)
    assert pos
    assert isinstance(pos, list)
    assert pos[0].contig_id == "12"
    assert pos[0].start == 25227248
    assert pos[0].end == 25227248
    assert pos[0].strand == "-"


def test_rna_to_dna_neg_strand_by_transcript_id_3():
    pos = rna_to_dna("ENST00000256078", 466, 501)
    assert pos
    assert isinstance(pos, list)
    assert pos[0].contig_id == "12"
    assert pos[0].start == 25225753
    assert pos[0].end == 25227248
    assert pos[0].strand == "-"


def test_rna_to_dna_neg_strand_by_transcript_name():
    pos = rna_to_dna("KRAS-201", 1)
    assert pos
    assert isinstance(pos, list)
    assert pos[0].contig_id == "12"
    assert pos[0].start == 25250929
    assert pos[0].end == 25250929
    assert pos[0].strand == "-"


def test_rna_to_dna_neg_strand_by_transcript_name_2():
    pos = rna_to_dna("KRAS-201", 466)
    assert pos
    assert isinstance(pos, list)
    assert pos[0].contig_id == "12"
    assert pos[0].start == 25227248
    assert pos[0].end == 25227248
    assert pos[0].strand == "-"


def test_rna_to_dna_neg_strand_by_transcript_name_3():
    pos = rna_to_dna("KRAS-201", 466, 501)
    assert pos
    assert isinstance(pos, list)
    assert pos[0].contig_id == "12"
    assert pos[0].start == 25225753
    assert pos[0].end == 25227248
    assert pos[0].strand == "-"


def test_rna_to_dna_by_refseq_id():
    pos = rna_to_dna("NM_001004491", 1)
    assert pos
    assert isinstance(pos, list)
    assert pos[0].contig_id == "1"
    assert pos[0].start == 247965233
    assert pos[0].end == 247965233
    assert pos[0].strand == "+"


def test_rna_to_dna_by_contig_id():
    # NOTE: I don't see this query being useful because it will return coordinates for every transcript on chr12.
    # result = sorted([i.contig_id for i in rna_to_dna("12", 301)])
    pass


def test_rna_to_dna_by_exon_id():
    result = sorted([i.contig_id for i in rna_to_dna("ENSE00000936617", 301)])
    assert result == ["12", "12", "12"]


def test_rna_to_dna_by_gene_id():
    result = sorted([i.contig_id for i in rna_to_dna("ENSG00000133703", 301)])
    assert result == ["12", "12", "12"]


def test_rna_to_dna_by_gene_name():
    result = sorted([i.contig_id for i in rna_to_dna("KRAS", 301)])
    print(rna_to_dna("KRAS", 301))
    assert result == ["12", "12", "12"]


def test_rna_to_dna_by_protein_id():
    result = sorted([i.contig_id for i in rna_to_dna("ENSP00000256078", 301)])
    assert result == ["12"]


def test_rna_to_exon_pos_strand_by_transcript_id():
    pos = rna_to_exon("ENST00000380152", 1)
    assert pos
    assert isinstance(pos, list)
    assert pos[0].transcript_id, "ENST00000380152"
    assert pos[0].exon_id, "ENSE00001184784"
    assert pos[0].start == 1
    assert pos[0].end == 1
    assert pos[0].strand == "+"


def test_rna_to_exon_pos_strand_by_transcript_id_2():
    pos = rna_to_exon("ENST00000380152", 500)
    assert pos
    assert isinstance(pos, list)
    assert pos[0].transcript_id, "ENST00000380152"
    assert pos[0].exon_id, "ENSE00003666217"
    assert pos[0].start == 3
    assert pos[0].end == 3
    assert pos[0].strand == "+"


def test_rna_to_exon_pos_strand_by_transcript_id_3():
    pos = rna_to_exon("ENST00000380152", 500, 510)
    assert pos
    assert isinstance(pos, list)
    assert pos[0].transcript_id, "ENST00000380152"
    assert pos[0].exon_id, "ENSE00003666217"
    assert pos[0].start == 3
    assert pos[0].end == 3
    assert pos[0].strand == "+"


def test_rna_to_exon_pos_strand_by_transcript_name():
    pos = rna_to_exon("BRCA2-201", 1)
    assert pos
    assert isinstance(pos, list)
    assert pos[0].transcript_id, "ENST00000380152"
    assert pos[0].exon_id, "ENSE00001184784"
    assert pos[0].start == 1
    assert pos[0].end == 1
    assert pos[0].strand == "+"


def test_rna_to_exon_pos_strand_by_transcript_name_2():
    pos = rna_to_exon("BRCA2-201", 500)
    assert pos
    assert isinstance(pos, list)
    assert pos[0].transcript_id, "ENST00000380152"
    assert pos[0].exon_id, "ENSE00003666217"
    assert pos[0].start == 3
    assert pos[0].end == 3
    assert pos[0].strand == "+"


def test_rna_to_exon_pos_strand_by_transcript_name_3():
    pos = rna_to_exon("BRCA2-201", 500, 510)
    assert pos
    assert isinstance(pos, list)
    assert pos[0].transcript_id, "ENST00000380152"
    assert pos[0].exon_id, "ENSE00003666217"
    assert pos[0].start == 3
    assert pos[0].end == 3
    assert pos[0].strand == "+"


def test_rna_to_exon_neg_strand_by_transcript_id():
    pos = rna_to_exon("ENST00000256078", 1)
    assert pos
    assert isinstance(pos, list)
    assert pos[0].transcript_id, "ENST00000256078"
    assert pos[0].exon_id, "ENSE00000000028"
    assert pos[0].start == 1
    assert pos[0].end == 1
    assert pos[0].strand == "-"


def test_rna_to_exon_neg_strand_by_transcript_id_2():
    pos = rna_to_exon("ENST00000256078", 400)
    assert pos
    assert isinstance(pos, list)
    assert pos[0].transcript_id, "ENST00000256078"
    assert pos[0].exon_id, "ENSE00001719809"
    assert pos[0].start == 3
    assert pos[0].end == 3
    assert pos[0].strand == "-"


def test_rna_to_exon_neg_strand_by_transcript_id_3():
    pos = rna_to_exon("ENST00000256078", 400, 420)
    assert pos
    assert isinstance(pos, list)
    assert pos[0].transcript_id, "ENST00000256078"
    assert pos[0].exon_id, "ENSE00001719809"
    assert pos[0].start == 3
    assert pos[0].end == 3
    assert pos[0].strand == "-"


def test_rna_to_exon_neg_strand_by_transcript_name():
    pos = rna_to_exon("KRAS-201", 1)
    assert pos
    assert isinstance(pos, list)
    assert pos[0].transcript_id, "ENST00000256078"
    assert pos[0].exon_id, "ENSE00000000028"
    assert pos[0].start == 1
    assert pos[0].end == 1
    assert pos[0].strand == "-"


def test_rna_to_exon_neg_strand_by_transcript_name_2():
    pos = rna_to_exon("KRAS-201", 400)
    assert pos
    assert isinstance(pos, list)
    assert pos[0].transcript_id, "ENST00000256078"
    assert pos[0].exon_id, "ENSE00001719809"
    assert pos[0].start == 3
    assert pos[0].end == 3
    assert pos[0].strand == "-"


def test_rna_to_exon_neg_strand_by_transcript_name_3():
    pos = rna_to_exon("KRAS-201", 400, 420)
    assert pos
    assert isinstance(pos, list)
    assert pos[0].transcript_id, "ENST00000256078"
    assert pos[0].exon_id, "ENSE00001719809"
    assert pos[0].start == 3
    assert pos[0].end == 3
    assert pos[0].strand == "-"


def test_rna_to_exon_by_contig_id():
    # NOTE: I don't see this query being useful because it will return coordinates for every transcript on chr12.
    # result = sorted([i.exon_id for i in rna_to_exon("12", 1)])
    pass


def test_rna_to_exon_by_exon_id():
    result = sorted([i.exon_id for i in rna_to_exon("ENSE00000936617", 1)])
    assert result == ["ENSE00000000028", "ENSE00001189804", "ENSE00002446502", "ENSE00002530521"]


def test_rna_to_exon_by_gene_id():
    result = sorted([i.exon_id for i in rna_to_exon("ENSG00000133703", 1)])
    assert result == ["ENSE00000000028", "ENSE00001189804", "ENSE00002446502", "ENSE00002530521"]


def test_rna_to_exon_by_gene_name():
    result = sorted([i.exon_id for i in rna_to_exon("KRAS", 1)])
    assert result == ["ENSE00000000028", "ENSE00001189804", "ENSE00002446502", "ENSE00002530521"]


def test_rna_to_exon_by_protein_id():
    result = sorted([i.exon_id for i in rna_to_exon("ENSP00000256078", 1)])
    assert result == ["ENSE00000000028"]


def test_rna_to_protein_pos_strand_by_transcript_id():
    pos = rna_to_protein("ENST00000288135", 98)
    assert pos
    assert isinstance(pos, list)
    assert pos[0].protein_id, "ENSP00000288135"
    assert pos[0].start == 14
    assert pos[0].end == 14
    assert pos[0].strand == "+"


def test_rna_to_protein_pos_strand_by_transcript_id_2():
    pos = rna_to_protein("ENST00000288135", 128)
    assert pos
    assert isinstance(pos, list)
    assert pos[0].protein_id, "ENSP00000288135"
    assert pos[0].start == 24
    assert pos[0].end == 24
    assert pos[0].strand == "+"


def test_rna_to_protein_pos_strand_by_transcript_id_3():
    pos = rna_to_protein("ENST00000288135", 128, 145)
    assert pos
    assert isinstance(pos, list)
    assert pos[0].protein_id, "ENSP00000288135"
    assert pos[0].start == 24
    assert pos[0].end == 29
    assert pos[0].strand == "+"


def test_rna_to_protein_pos_strand_by_transcript_name():
    pos = rna_to_protein("KIT-201", 98)
    assert pos
    assert isinstance(pos, list)
    assert pos[0].protein_id, "ENSP00000288135"
    assert pos[0].start == 14
    assert pos[0].end == 14
    assert pos[0].strand == "+"


def test_rna_to_protein_pos_strand_by_transcript_name_2():
    pos = rna_to_protein("KIT-201", 128)
    assert pos
    assert isinstance(pos, list)
    assert pos[0].protein_id, "ENSP00000288135"
    assert pos[0].start == 24
    assert pos[0].end == 24
    assert pos[0].strand == "+"


def test_rna_to_protein_pos_strand_by_transcript_name_3():
    pos = rna_to_protein("KIT-201", 128, 145)
    assert pos
    assert isinstance(pos, list)
    assert pos[0].protein_id, "ENSP00000288135"
    assert pos[0].start == 24
    assert pos[0].end == 29
    assert pos[0].strand == "+"


def test_rna_to_protein_neg_strand_by_transcript_id():
    pos = rna_to_protein("ENST00000260947", 115)
    assert pos
    assert isinstance(pos, list)
    assert pos[0].protein_id, "ENSP00000260947"
    assert pos[0].start == 1
    assert pos[0].end == 1
    assert pos[0].strand == "-"


def test_rna_to_protein_neg_strand_by_transcript_id_2():
    pos = rna_to_protein("ENST00000260947", 346)
    assert pos
    assert isinstance(pos, list)
    assert pos[0].protein_id, "ENSP00000260947"
    assert pos[0].start == 78
    assert pos[0].end == 78
    assert pos[0].strand == "-"


def test_rna_to_protein_neg_strand_by_transcript_id_3():
    pos = rna_to_protein("ENST00000260947", 345, 356)
    assert pos
    assert isinstance(pos, list)
    assert pos[0].protein_id, "ENSP00000260947"
    assert pos[0].start == 77
    assert pos[0].end == 81
    assert pos[0].strand == "-"


def test_rna_to_protein_neg_strand_by_transcript_name():
    pos = rna_to_protein("BARD1-201", 115)
    assert pos
    assert isinstance(pos, list)
    assert pos[0].protein_id, "ENSP00000260947"
    assert pos[0].start == 1
    assert pos[0].end == 1
    assert pos[0].strand == "-"


def test_rna_to_protein_neg_strand_by_transcript_name_2():
    pos = rna_to_protein("BARD1-201", 346)
    assert pos
    assert isinstance(pos, list)
    assert pos[0].protein_id, "ENSP00000260947"
    assert pos[0].start == 78
    assert pos[0].end == 78
    assert pos[0].strand == "-"


def test_rna_to_protein_neg_strand_by_transcript_name_3():
    pos = rna_to_protein("BARD1-201", 345, 356)
    assert pos
    assert isinstance(pos, list)
    assert pos[0].protein_id, "ENSP00000260947"
    assert pos[0].start == 77
    assert pos[0].end == 81
    assert pos[0].strand == "-"


def test_rna_to_protein_by_contig_id():
    # NOTE: I don't see this query being useful because it will return coordinates for every transcript on chr12.
    # result = sorted([i.protein_id for i in rna_to_protein("12", 1)])
    pass


def test_rna_to_protein_by_exon_id():
    result = sorted([i.protein_id for i in rna_to_protein("ENSE00000936617", 301)])
    assert result == ["ENSP00000256078", "ENSP00000308495", "ENSP00000451856", "ENSP00000452512"]


def test_rna_to_protein_by_gene_id():
    result = sorted([i.protein_id for i in rna_to_protein("ENSG00000133703", 301)])
    assert result == ["ENSP00000256078", "ENSP00000308495", "ENSP00000451856", "ENSP00000452512"]


def test_rna_to_protein_by_gene_name():
    result = sorted([i.protein_id for i in rna_to_protein("KRAS", 301)])
    assert result == ["ENSP00000256078", "ENSP00000308495", "ENSP00000451856", "ENSP00000452512"]


def test_rna_to_protein_by_protein_id():
    result = sorted([i.protein_id for i in rna_to_protein("ENSP00000256078", 301)])
    assert result == ["ENSP00000256078"]


def test_rna_to_rna_pos_strand_by_transcript_id():
    pos = rna_to_rna("ENST00000288135", 1)
    assert pos
    assert isinstance(pos, list)
    assert pos[0].transcript_id, "ENST00000288135"
    assert pos[0].start == 1
    assert pos[0].end == 1
    assert pos[0].strand == "+"


def test_rna_to_rna_pos_strand_by_transcript_id_2():
    pos = rna_to_rna("ENST00000288135", 100)
    assert pos
    assert isinstance(pos, list)
    assert pos[0].transcript_id, "ENST00000288135"
    assert pos[0].start == 100
    assert pos[0].end == 100
    assert pos[0].strand == "+"


def test_rna_to_rna_pos_strand_by_transcript_id_3():
    pos = rna_to_rna("ENST00000288135", 98, 100)
    assert pos
    assert isinstance(pos, list)
    assert pos[0].transcript_id, "ENST00000288135"
    assert pos[0].start == 98
    assert pos[0].end == 100
    assert pos[0].strand == "+"


def test_rna_to_rna_pos_strand_by_transcript_name():
    pos = rna_to_rna("KIT-201", 1)
    assert pos
    assert isinstance(pos, list)
    assert pos[0].transcript_id, "ENST00000288135"
    assert pos[0].start == 1
    assert pos[0].end == 1
    assert pos[0].strand == "+"


def test_rna_to_rna_pos_strand_by_transcript_name_2():
    pos = rna_to_rna("KIT-201", 100)
    assert pos
    assert isinstance(pos, list)
    assert pos[0].transcript_id, "ENST00000288135"
    assert pos[0].start == 100
    assert pos[0].end == 100
    assert pos[0].strand == "+"


def test_rna_to_rna_pos_strand_by_transcript_name_3():
    pos = rna_to_rna("KIT-201", 98, 100)
    assert pos
    assert isinstance(pos, list)
    assert pos[0].transcript_id, "ENST00000288135"
    assert pos[0].start == 98
    assert pos[0].end == 100
    assert pos[0].strand == "+"


def test_rna_to_rna_neg_strand_by_transcript_id():
    pos = rna_to_rna("ENST00000256078", 1)
    assert pos
    assert isinstance(pos, list)
    assert pos[0].transcript_id, "ENST00000256078"
    assert pos[0].start == 1
    assert pos[0].end == 1
    assert pos[0].strand == "-"


def test_rna_to_rna_neg_strand_by_transcript_id_2():
    pos = rna_to_rna("ENST00000256078", 111)
    assert pos
    assert isinstance(pos, list)
    assert pos[0].transcript_id, "ENST00000256078"
    assert pos[0].start == 111
    assert pos[0].end == 111
    assert pos[0].strand == "-"


def test_rna_to_rna_neg_strand_by_transcript_id_3():
    pos = rna_to_rna("ENST00000256078", 111, 133)
    assert pos
    assert isinstance(pos, list)
    assert pos[0].transcript_id, "ENST00000256078"
    assert pos[0].start == 111
    assert pos[0].end == 133
    assert pos[0].strand == "-"


def test_rna_to_rna_neg_strand_by_transcript_name():
    pos = rna_to_rna("KRAS-201", 1)
    assert pos
    assert isinstance(pos, list)
    assert pos[0].transcript_id, "ENST00000256078"
    assert pos[0].start == 1
    assert pos[0].end == 1
    assert pos[0].strand == "-"


def test_rna_to_rna_neg_strand_by_transcript_name_2():
    pos = rna_to_rna("KRAS-201", 111)
    assert pos
    assert isinstance(pos, list)
    assert pos[0].transcript_id, "ENST00000256078"
    assert pos[0].start == 111
    assert pos[0].end == 111
    assert pos[0].strand == "-"


def test_rna_to_rna_neg_strand_by_transcript_name_3():
    pos = rna_to_rna("KRAS-201", 111, 133)
    assert pos
    assert isinstance(pos, list)
    assert pos[0].transcript_id, "ENST00000256078"
    assert pos[0].start == 111
    assert pos[0].end == 133
    assert pos[0].strand == "-"


def test_rna_to_rna_by_contig_id():
    # NOTE: I don't see this query being useful because it will return coordinates for essentially every transcript on chr12.
    # result = sorted([i.transcript_id for i in rna_to_rna("12", 301)])
    pass


def test_rna_to_rna_by_exon_id():
    result = sorted([i.transcript_id for i in rna_to_rna("ENSE00000936617", 301)])
    assert result == ["ENST00000256078", "ENST00000311936", "ENST00000556131", "ENST00000557334"]


def test_rna_to_rna_by_gene_id():
    result = sorted([i.transcript_id for i in rna_to_rna("ENSG00000133703", 301)])
    assert result == ["ENST00000256078", "ENST00000311936", "ENST00000556131", "ENST00000557334"]


def test_rna_to_rna_by_gene_name():
    result = sorted([i.transcript_id for i in rna_to_rna("KRAS", 301)])
    assert result == ["ENST00000256078", "ENST00000311936", "ENST00000556131", "ENST00000557334"]


def test_rna_to_rna_by_protein_id():
    result = sorted([i.transcript_id for i in rna_to_rna("ENSP00000256078", 301)])
    assert result == ["ENST00000256078"]
