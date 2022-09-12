import pytest

from ensembl_map.constants import CONTIG, EXON, GENE, PROTEIN, TRANSCRIPT
from ensembl_map.core import EnsemblRelease
from . import (
    CACHE_DIR,
    CANONICAL_TRANSCRIPT,
    CONTIG_ALIAS,
    EXON_ALIAS,
    GENE_ALIAS,
    PROTEIN_ALIAS,
    TRANSCRIPT_ALIAS,
)

EXON_ID = "ENSE00003826864"
CONTIG_ID = "16"
GENE_ID = "ENSG00000149925"
GENE_NAME = "ALDOA"
PROTEIN_ID = "ENSP00000494188"  # version 2
TRANSCRIPT_ID = "ENST00000643777"  # version 3
TRANSCRIPT_NAME = "ALDOA-219"


def test_unique_instances():
    a = EnsemblRelease(species="homo_sapiens", reference="GRCh38", release=100, cache_dir=CACHE_DIR)
    b = EnsemblRelease(species="homo_sapiens", reference="GRCh38", release=100, cache_dir=CACHE_DIR)
    c = EnsemblRelease(species="homo_sapiens", reference="GRCh37", release=69, cache_dir=CACHE_DIR)
    assert a is b
    assert a is not c


@pytest.fixture
def ensembl69():
    return EnsemblRelease(
        species="homo_sapiens", reference="GRCh37", release=69, cache_dir=CACHE_DIR
    )


@pytest.fixture
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


def test_all_contig_ids(ensembl100):
    result = ensembl100.all_contig_ids()
    assert isinstance(result, list)
    assert CONTIG_ID in result


def test_all_exon_ids(ensembl100):
    result = ensembl100.all_exon_ids()
    assert isinstance(result, list)
    assert EXON_ID in result


def test_all_gene_ids(ensembl100):
    result = ensembl100.all_gene_ids()
    assert isinstance(result, list)
    assert GENE_ID in result


def test_all_gene_names(ensembl100):
    result = ensembl100.all_gene_names()
    assert isinstance(result, list)
    assert GENE_NAME in result


def test_all_protein_ids(ensembl100):
    result = ensembl100.all_protein_ids()
    assert isinstance(result, list)
    assert PROTEIN_ID in result


def test_all_transcript_ids(ensembl100):
    result = ensembl100.all_transcript_ids()
    assert isinstance(result, list)
    assert TRANSCRIPT_ID in result


def test_all_transcript_names(ensembl100):
    result = ensembl100.all_transcript_names()
    assert isinstance(result, list)
    assert TRANSCRIPT_NAME in result


def test_contig_ids_from_contig(ensembl100):
    ret_by_id = ensembl100.contig_ids(CONTIG_ID)
    assert CONTIG_ID in ret_by_id


def test_contig_ids_from_contig_id(ensembl100):
    result = ensembl100.contig_ids(CONTIG_ID)
    assert isinstance(result, list)
    assert CONTIG_ID in result


def test_contig_ids_from_contig_with_chr(ensembl100):
    ret_by_id = ensembl100.contig_ids("chr12")
    assert "12" in ret_by_id


def test_contig_ids_from_exon(ensembl100):
    ret_by_id = ensembl100.contig_ids(EXON_ID, EXON)
    assert CONTIG_ID in ret_by_id


def test_contig_ids_from_exon_id(ensembl100):
    result = ensembl100.contig_ids(EXON_ID, EXON)
    assert isinstance(result, list)
    assert CONTIG_ID in result


def test_contig_ids_from_gene(ensembl100):
    ret_by_id = ensembl100.contig_ids(GENE_ID)
    ret_by_name = ensembl100.contig_ids(GENE_NAME)
    assert CONTIG_ID in ret_by_id
    assert CONTIG_ID in ret_by_name
    assert ret_by_id == ret_by_name


def test_contig_ids_from_gene_id(ensembl100):
    result = ensembl100.contig_ids(GENE_ID)
    assert isinstance(result, list)
    assert CONTIG_ID in result


def test_contig_ids_from_gene_name(ensembl100):
    result = ensembl100.contig_ids(GENE_NAME)
    assert isinstance(result, list)
    assert CONTIG_ID in result


def test_contig_ids_from_protein(ensembl100):
    ret_by_id = ensembl100.contig_ids(PROTEIN_ID)
    assert CONTIG_ID in ret_by_id


def test_contig_ids_from_protein_id(ensembl100):
    result = ensembl100.contig_ids(PROTEIN_ID)
    assert isinstance(result, list)
    assert CONTIG_ID in result


def test_contig_ids_from_transcript(ensembl100):
    ret_by_id = ensembl100.contig_ids(TRANSCRIPT_ID)
    ret_by_name = ensembl100.contig_ids(TRANSCRIPT_NAME)
    assert CONTIG_ID in ret_by_id
    assert CONTIG_ID in ret_by_name
    assert ret_by_id == ret_by_name


def test_contig_ids_from_transcript_id(ensembl100):
    result = ensembl100.contig_ids(TRANSCRIPT_ID)
    assert isinstance(result, list)
    assert CONTIG_ID in result


def test_contig_ids_from_transcript_name(ensembl100):
    result = ensembl100.contig_ids(TRANSCRIPT_NAME)
    assert isinstance(result, list)
    assert CONTIG_ID in result


def test_exon_ids_from_contig_id(ensembl100):
    result = ensembl100.exon_ids(CONTIG_ID)
    assert isinstance(result, list)
    assert EXON_ID in result


def test_exon_ids_from_exon(ensembl100):
    ret_by_id = ensembl100.exon_ids(EXON_ID, EXON)
    assert EXON_ID in ret_by_id


def test_exon_ids_from_exon_id(ensembl100):
    result = ensembl100.exon_ids(EXON_ID)
    assert isinstance(result, list)
    assert EXON_ID in result


def test_exon_ids_from_gene(ensembl100):
    ret_by_id = ensembl100.exon_ids(GENE_ID)
    ret_by_name = ensembl100.exon_ids(GENE_NAME)
    assert EXON_ID in ret_by_id
    assert EXON_ID in ret_by_name
    assert ret_by_id == ret_by_name


def test_exon_ids_from_gene_id(ensembl100):
    result = ensembl100.exon_ids(GENE_ID)
    assert isinstance(result, list)
    assert EXON_ID in result


def test_exon_ids_from_gene_name(ensembl100):
    result = ensembl100.exon_ids(GENE_NAME)
    assert isinstance(result, list)
    assert EXON_ID in result


def test_exon_ids_from_protein(ensembl100):
    ret_by_id = ensembl100.exon_ids(PROTEIN_ID)
    assert EXON_ID in ret_by_id


def test_exon_ids_from_protein_id(ensembl100):
    result = ensembl100.exon_ids(PROTEIN_ID)
    assert isinstance(result, list)
    assert EXON_ID in result


def test_exon_ids_from_transcript(ensembl100):
    ret_by_id = ensembl100.exon_ids(TRANSCRIPT_ID)
    ret_by_name = ensembl100.exon_ids(TRANSCRIPT_NAME)
    assert EXON_ID in ret_by_id
    assert EXON_ID in ret_by_name
    assert ret_by_id == ret_by_name


def test_exon_ids_from_transcript_id(ensembl100):
    result = ensembl100.exon_ids(TRANSCRIPT_ID)
    assert isinstance(result, list)
    assert EXON_ID in result


def test_exon_ids_from_transcript_name(ensembl100):
    result = ensembl100.exon_ids(TRANSCRIPT_NAME)
    assert isinstance(result, list)
    assert EXON_ID in result


def test_gene_ids_different_gene_ids(ensembl69, ensembl100):
    # Gene ID changed bewteen 69 and 100
    # ENSG00000166748 AGBL1
    # ENSG00000273540 AGBL1
    assert ensembl69.gene_ids("AGBL1") == ["ENSG00000166748"]
    assert ensembl100.gene_ids("AGBL1") == ["ENSG00000273540"]


def test_gene_ids_different_gene_names(ensembl69, ensembl100):
    # Gene name changed bewteen 69 and 100
    # ENSG00000118058 MLL
    # ENSG00000118058 KMT2A
    assert ensembl69.gene_ids("MLL") == ["ENSG00000118058"]
    assert ensembl100.gene_ids("KMT2A") == ["ENSG00000118058"]


def test_gene_ids_from_contig_id(ensembl100):
    result = ensembl100.gene_ids(CONTIG_ID)
    assert isinstance(result, list)
    assert GENE_ID in result


def test_gene_ids_from_exon(ensembl100):
    ret_by_id = ensembl100.gene_ids(EXON_ID)
    assert GENE_ID in ret_by_id


def test_gene_ids_from_exon_id(ensembl100):
    result = ensembl100.gene_ids(EXON_ID)
    assert isinstance(result, list)
    assert GENE_ID in result


def test_gene_ids_from_gene(ensembl100):
    ret_by_id = ensembl100.gene_ids(GENE_ID)
    ret_by_name = ensembl100.gene_ids(GENE_NAME)
    assert GENE_ID in ret_by_id
    assert GENE_ID in ret_by_name
    assert ret_by_id == ret_by_name


def test_gene_ids_from_gene_id(ensembl100):
    result = ensembl100.gene_ids(GENE_ID)
    assert isinstance(result, list)
    assert GENE_ID in result


def test_gene_ids_from_gene_name(ensembl100):
    result = ensembl100.gene_ids(GENE_NAME)
    assert isinstance(result, list)
    assert GENE_ID in result


def test_gene_ids_from_protein(ensembl100):
    ret_by_id = ensembl100.gene_ids(PROTEIN_ID)
    assert GENE_ID in ret_by_id


def test_gene_ids_from_protein_id(ensembl100):
    result = ensembl100.gene_ids(PROTEIN_ID)
    assert isinstance(result, list)
    assert GENE_ID in result


def test_gene_ids_from_transcript(ensembl100):
    ret_by_id = ensembl100.gene_ids(TRANSCRIPT_ID)
    ret_by_name = ensembl100.gene_ids(TRANSCRIPT_NAME)
    assert GENE_ID in ret_by_id
    assert GENE_ID in ret_by_name
    assert ret_by_id == ret_by_name


def test_gene_ids_from_transcript_id(ensembl100):
    result = ensembl100.gene_ids(TRANSCRIPT_ID)
    assert isinstance(result, list)
    assert GENE_ID in result


def test_gene_ids_from_transcript_name(ensembl100):
    result = ensembl100.gene_ids(TRANSCRIPT_NAME)
    assert isinstance(result, list)
    assert GENE_ID in result


def test_gene_names_different_gene_ids(ensembl69, ensembl100):
    # Gene ID changed bewteen 69 and 100
    # ENSG00000204645 SSX4
    # ENSG00000268009 SSX4
    assert ensembl69.gene_names("ENSG00000204645") == ["SSX4"]
    assert ensembl100.gene_names("ENSG00000268009") == ["SSX4"]


def test_gene_names_different_gene_names(ensembl69, ensembl100):
    # Gene name changed bewteen 69 and 100
    # ENSG00000118058 MLL
    # ENSG00000118058 KMT2A
    assert ensembl69.gene_names("ENSG00000118058") == ["MLL"]
    assert ensembl100.gene_names("ENSG00000118058") == ["KMT2A"]


def test_gene_names_from_contig_id(ensembl100):
    result = ensembl100.gene_names(CONTIG_ID)
    assert isinstance(result, list)
    assert GENE_NAME in result


def test_gene_names_from_exon(ensembl100):
    ret_by_id = ensembl100.gene_names(EXON_ID)
    assert GENE_NAME in ret_by_id


def test_gene_names_from_exon_id(ensembl100):
    result = ensembl100.gene_names(EXON_ID)
    assert isinstance(result, list)
    assert GENE_NAME in result


def test_gene_names_from_gene(ensembl100):
    ret_by_id = ensembl100.gene_names(GENE_ID)
    ret_by_name = ensembl100.gene_names(GENE_NAME)
    assert GENE_NAME in ret_by_id
    assert GENE_NAME in ret_by_name
    assert ret_by_id == ret_by_name


def test_gene_names_from_gene_id(ensembl100):
    result = ensembl100.gene_names(GENE_ID)
    assert isinstance(result, list)
    assert GENE_NAME in result


def test_gene_names_from_gene_name(ensembl100):
    result = ensembl100.gene_names(GENE_NAME)
    assert isinstance(result, list)
    assert GENE_NAME in result


def test_gene_names_from_protein(ensembl100):
    ret_by_id = ensembl100.gene_names(PROTEIN_ID)
    assert GENE_NAME in ret_by_id


def test_gene_names_from_protein_id(ensembl100):
    result = ensembl100.gene_names(PROTEIN_ID)
    assert isinstance(result, list)
    assert GENE_NAME in result


def test_gene_names_from_refseq_transcript(ensembl100):
    ret_by_id = ensembl100.gene_names("NM_000244")
    assert "MEN1" in ret_by_id
    ret_by_id = ensembl100.gene_names("NM_007194")
    assert "CHEK2" in ret_by_id
    ret_by_id = ensembl100.gene_names("NM_001128849")
    assert "SMARCA4" in ret_by_id
    ret_by_id = ensembl100.gene_names("NM_000314")
    assert "PTEN" in ret_by_id


def test_gene_names_from_transcript(ensembl100):
    ret_by_id = ensembl100.gene_names(TRANSCRIPT_ID)
    ret_by_name = ensembl100.gene_names(TRANSCRIPT_NAME)
    assert GENE_NAME in ret_by_id
    assert GENE_NAME in ret_by_name
    assert ret_by_id == ret_by_name


def test_gene_names_from_transcript_id(ensembl100):
    result = ensembl100.gene_names(TRANSCRIPT_ID)
    assert isinstance(result, list)
    assert GENE_NAME in result


def test_gene_names_from_transcript_name(ensembl100):
    result = ensembl100.gene_names(TRANSCRIPT_NAME)
    assert isinstance(result, list)
    assert GENE_NAME in result


def test_is_canonical_transcript(ensembl100):
    assert ensembl100.is_canonical_transcript("ENST00000000233") is True
    assert ensembl100.is_canonical_transcript("spam") is False


def test_is_contig_false(ensembl100):
    assert ensembl100.is_contig("spam") is False


def test_is_contig_true(ensembl100):
    assert ensembl100.is_contig("12") is True


def test_is_contig_true_alias(ensembl100):
    assert ensembl100.is_contig("chr12") is True


def test_is_exon_false(ensembl100):
    assert ensembl100.is_exon("spam") is False


def test_is_exon_true(ensembl100):
    assert ensembl100.is_exon(EXON_ID) is True


def test_is_gene_false(ensembl100):
    assert ensembl100.is_gene("spam") is False


def test_is_gene_true_from_id(ensembl100):
    assert ensembl100.is_gene(GENE_ID) is True


def test_is_gene_true_from_name(ensembl100):
    assert ensembl100.is_gene(GENE_NAME) is True


def test_is_protein_false(ensembl100):
    assert ensembl100.is_protein("spam") is False


def test_is_protein_true(ensembl100):
    assert ensembl100.is_protein(PROTEIN_ID) is True


def test_is_transcript_false(ensembl100):
    assert ensembl100.is_transcript("spam") is False


def test_is_transcript_true_from_id(ensembl100):
    assert ensembl100.is_transcript(TRANSCRIPT_ID) is True


def test_is_transcript_true_from_name(ensembl100):
    assert ensembl100.is_transcript(TRANSCRIPT_NAME) is True


def test_normalize_feature_contig_id(ensembl100):
    assert ensembl100.normalize_feature("17") == [("17", CONTIG)]


def test_normalize_feature_exon_id(ensembl100):
    assert ensembl100.normalize_feature(EXON_ID) == [(EXON_ID, EXON)]


def test_normalize_feature_gene_id(ensembl100):
    assert ensembl100.normalize_feature(GENE_ID) == [(GENE_ID, GENE)]


def test_normalize_feature_gene_name(ensembl100):
    assert ensembl100.normalize_feature(GENE_NAME) == [(GENE_NAME, GENE)]


def test_normalize_feature_protein_id(ensembl100):
    assert ensembl100.normalize_feature(PROTEIN_ID) == [(PROTEIN_ID, PROTEIN)]


def test_normalize_feature_refseq_transcript_id(ensembl100):
    assert ensembl100.normalize_feature("NM_000314.4") == [("ENST00000371953", TRANSCRIPT)]


def test_normalize_feature_transcript_id(ensembl100):
    assert ensembl100.normalize_feature(TRANSCRIPT_ID) == [(TRANSCRIPT_ID, TRANSCRIPT)]


def test_normalize_feature_transcript_name(ensembl100):
    assert ensembl100.normalize_feature(TRANSCRIPT_NAME) == [(TRANSCRIPT_NAME, TRANSCRIPT)]


def test_protein_ids_from_contig_id(ensembl100):
    result = ensembl100.protein_ids(CONTIG_ID)
    assert isinstance(result, list)
    assert PROTEIN_ID in result


def test_protein_ids_from_exon(ensembl100):
    ret_by_id = ensembl100.protein_ids(EXON_ID, EXON)
    assert PROTEIN_ID in ret_by_id


def test_protein_ids_from_exon_id(ensembl100):
    result = ensembl100.protein_ids(EXON_ID)
    assert isinstance(result, list)
    assert PROTEIN_ID in result


def test_protein_ids_from_gene(ensembl100):
    ret_by_id = ensembl100.protein_ids(GENE_ID)
    ret_by_name = ensembl100.protein_ids(GENE_NAME)
    assert PROTEIN_ID in ret_by_id
    assert PROTEIN_ID in ret_by_name
    assert ret_by_id == ret_by_name


def test_protein_ids_from_gene_id(ensembl100):
    result = ensembl100.protein_ids(GENE_ID)
    assert isinstance(result, list)
    assert PROTEIN_ID in result


def test_protein_ids_from_gene_name(ensembl100):
    result = ensembl100.protein_ids(GENE_NAME)
    assert isinstance(result, list)
    assert PROTEIN_ID in result


def test_protein_ids_from_protein(ensembl100):
    ret_by_id = ensembl100.protein_ids(PROTEIN_ID)
    assert PROTEIN_ID in ret_by_id


def test_protein_ids_from_protein_id(ensembl100):
    result = ensembl100.protein_ids(PROTEIN_ID)
    assert isinstance(result, list)
    assert PROTEIN_ID in result


def test_protein_ids_from_transcript(ensembl100):
    ret_by_id = ensembl100.protein_ids(TRANSCRIPT_ID)
    ret_by_name = ensembl100.protein_ids(TRANSCRIPT_NAME)
    assert PROTEIN_ID in ret_by_id
    assert PROTEIN_ID in ret_by_name
    assert ret_by_id == ret_by_name


def test_protein_ids_from_transcript_id(ensembl100):
    result = ensembl100.protein_ids(TRANSCRIPT_ID)
    assert isinstance(result, list)
    assert PROTEIN_ID in result


def test_protein_ids_from_transcript_name(ensembl100):
    result = ensembl100.protein_ids(TRANSCRIPT_NAME)
    assert isinstance(result, list)
    assert PROTEIN_ID in result


def test_transcript_ids_from_contig_id(ensembl100):
    result = ensembl100.transcript_ids(CONTIG_ID)
    assert isinstance(result, list)
    assert TRANSCRIPT_ID in result


def test_transcript_ids_from_exon(ensembl100):
    ret_by_id = ensembl100.transcript_ids(EXON_ID, EXON)
    assert TRANSCRIPT_ID in ret_by_id


def test_transcript_ids_from_exon_id(ensembl100):
    result = ensembl100.transcript_ids(EXON_ID)
    assert isinstance(result, list)
    assert TRANSCRIPT_ID in result


def test_transcript_ids_from_gene(ensembl100):
    ret_by_id = ensembl100.transcript_ids(GENE_ID)
    ret_by_name = ensembl100.transcript_ids(GENE_NAME)
    assert TRANSCRIPT_ID in ret_by_id
    assert TRANSCRIPT_ID in ret_by_name
    assert ret_by_id == ret_by_name


@pytest.mark.xfail(reason="unsure if this will be implemented")
def test_transcript_ids_from_gene_best_only(ensembl100):
    # assumes ENST00000288135 is the best transcript for ENSG00000157404 (KIT)
    ret_all = ensembl100.transcript_ids(GENE_ID, best_only=False)
    ret_best = ensembl100.transcript_ids(GENE_ID, best_only=True)
    assert "ENST00000412167" in ret_all
    assert "ENST00000412167" not in ret_best


def test_transcript_ids_from_gene_id(ensembl100):
    result = ensembl100.transcript_ids(GENE_ID)
    assert isinstance(result, list)
    assert TRANSCRIPT_ID in result


def test_transcript_ids_from_gene_name(ensembl100):
    result = ensembl100.transcript_ids(GENE_NAME)
    assert isinstance(result, list)
    assert TRANSCRIPT_ID in result


def test_transcript_ids_from_protein(ensembl100):
    ret_by_id = ensembl100.transcript_ids(PROTEIN_ID)
    assert TRANSCRIPT_ID in ret_by_id


def test_transcript_ids_from_protein_id(ensembl100):
    result = ensembl100.transcript_ids(PROTEIN_ID)
    assert isinstance(result, list)
    assert TRANSCRIPT_ID in result


def test_transcript_ids_from_transcript(ensembl100):
    ret_by_id = ensembl100.transcript_ids(TRANSCRIPT_ID)
    ret_by_name = ensembl100.transcript_ids(TRANSCRIPT_NAME)
    assert TRANSCRIPT_ID in ret_by_id
    assert TRANSCRIPT_ID in ret_by_name
    assert ret_by_id == ret_by_name


def test_transcript_ids_from_transcript_id(ensembl100):
    result = ensembl100.transcript_ids(TRANSCRIPT_ID)
    assert isinstance(result, list)
    assert TRANSCRIPT_ID in result


def test_transcript_ids_from_transcript_name(ensembl100):
    result = ensembl100.transcript_ids(TRANSCRIPT_NAME)
    assert isinstance(result, list)
    assert TRANSCRIPT_ID in result


def test_transcript_names_from_contig_id(ensembl100):
    result = ensembl100.transcript_names(CONTIG_ID)
    assert isinstance(result, list)
    assert TRANSCRIPT_NAME in result


def test_transcript_names_from_exon(ensembl100):
    ret_by_id = ensembl100.transcript_names(EXON_ID, EXON)
    assert TRANSCRIPT_NAME in ret_by_id


def test_transcript_names_from_exon_id(ensembl100):
    result = ensembl100.transcript_names(EXON_ID)
    assert isinstance(result, list)
    assert TRANSCRIPT_NAME in result


def test_transcript_names_from_gene(ensembl100):
    ret_by_id = ensembl100.transcript_names(GENE_ID)
    ret_by_name = ensembl100.transcript_names(GENE_NAME)
    assert TRANSCRIPT_NAME in ret_by_id
    assert TRANSCRIPT_NAME in ret_by_name
    assert ret_by_id == ret_by_name


def test_transcript_names_from_gene_id(ensembl100):
    result = ensembl100.transcript_names(GENE_ID)
    assert isinstance(result, list)
    assert TRANSCRIPT_NAME in result


def test_transcript_names_from_gene_name(ensembl100):
    result = ensembl100.transcript_names(GENE_NAME)
    assert isinstance(result, list)
    assert TRANSCRIPT_NAME in result


def test_transcript_names_from_protein(ensembl100):
    ret_by_id = ensembl100.transcript_names(PROTEIN_ID)
    assert TRANSCRIPT_NAME in ret_by_id


def test_transcript_names_from_protein_id(ensembl100):
    result = ensembl100.transcript_names(PROTEIN_ID)
    assert isinstance(result, list)
    assert TRANSCRIPT_NAME in result


def test_transcript_names_from_transcript(ensembl100):
    ret_by_id = ensembl100.transcript_names(TRANSCRIPT_ID)
    ret_by_name = ensembl100.transcript_names(TRANSCRIPT_NAME)
    assert TRANSCRIPT_NAME in ret_by_id
    assert TRANSCRIPT_NAME in ret_by_name
    assert ret_by_id == ret_by_name


def test_transcript_names_from_transcript_id(ensembl100):
    result = ensembl100.transcript_names(TRANSCRIPT_ID)
    assert isinstance(result, list)
    assert TRANSCRIPT_NAME in result
