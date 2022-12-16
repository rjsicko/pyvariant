import pandas as pd
import pytest
from constants import (
    CACHE_DIR,
    CANONICAL_TRANSCRIPT,
    CONTIG_ALIAS,
    EXON_ALIAS,
    GENE_ALIAS,
    PROTEIN_ALIAS,
    TEST_CDNA_FASTA,
    TEST_DNA_FASTA,
    TEST_GTF,
    TEST_NCRNA_FASTA,
    TEST_PEP_FASTA,
    TRANSCRIPT_ALIAS,
)
from pyfaidx import Fasta

from ensembl_map.constants import (
    CONTIG_ID,
    EXON_ID,
    GENE_ID,
    GENE_NAME,
    PROTEIN_ID,
    TRANSCRIPT_ID,
    TRANSCRIPT_NAME,
)
from ensembl_map.core import (
    CdnaPosition,
    Core,
    DnaPosition,
    EnsemblRelease,
    ExonPosition,
    ProteinPosition,
    RnaPosition,
)

# features are from Ensembl 100
TEST_EXON_ID = "ENSE00003826864"
TEST_CONTIG_ID = "16"
TEST_GENE_ID = "ENSG00000149925"
TEST_GENE_NAME = "ALDOA"
TEST_PROTEIN_ID = "ENSP00000494188"  # version 2
TEST_TRANSCRIPT_ID = "ENST00000643777"  # version 3
TEST_TRANSCRIPT_NAME = "ALDOA-219"


# -------------------------------------------------------------------------------------------------
# custom core
# -------------------------------------------------------------------------------------------------
@pytest.mark.skip(reason="takes too long")
def test_custom_core():
    obj = Core(
        gtf=TEST_GTF,
        cdna=TEST_CDNA_FASTA,
        dna=TEST_DNA_FASTA,
        pep=TEST_PEP_FASTA,
        ncrna=TEST_NCRNA_FASTA,
        canonical_transcript=CANONICAL_TRANSCRIPT,
        contig_alias=CONTIG_ALIAS,
        exon_alias=EXON_ALIAS,
        gene_alias=GENE_ALIAS,
        protein_alias=PROTEIN_ALIAS,
        transcript_alias=TRANSCRIPT_ALIAS,
    )
    assert isinstance(obj.df, pd.DataFrame)
    assert isinstance(obj.cdna, Fasta)
    assert isinstance(obj.dna, Fasta)
    assert isinstance(obj.pep, Fasta)
    assert isinstance(obj.ncrna, Fasta)
    assert isinstance(obj.canonical_transcript, list)
    assert isinstance(obj.contig_alias, dict)
    assert isinstance(obj.exon_alias, dict)
    assert isinstance(obj.gene_alias, dict)
    assert isinstance(obj.protein_alias, dict)
    assert isinstance(obj.transcript_alias, dict)


# -------------------------------------------------------------------------------------------------
# cache logic
# -------------------------------------------------------------------------------------------------
def test_unique_instances():
    a = EnsemblRelease(species="homo_sapiens", release=100, cache_dir=CACHE_DIR)
    b = EnsemblRelease(species="homo_sapiens", release=100, cache_dir=CACHE_DIR)
    c = EnsemblRelease(species="homo_sapiens", release=69, cache_dir=CACHE_DIR)
    assert a is b
    assert a is not c


# -------------------------------------------------------------------------------------------------
# <feature_symbol>()
# -------------------------------------------------------------------------------------------------
def test_all_contig_ids(ensembl100):
    result = ensembl100.all_contig_ids()
    assert isinstance(result, list)
    assert TEST_CONTIG_ID in result


def test_all_exon_ids(ensembl100):
    result = ensembl100.all_exon_ids()
    assert isinstance(result, list)
    assert TEST_EXON_ID in result


def test_all_gene_ids(ensembl100):
    result = ensembl100.all_gene_ids()
    assert isinstance(result, list)
    assert TEST_GENE_ID in result


def test_all_gene_names(ensembl100):
    result = ensembl100.all_gene_names()
    assert isinstance(result, list)
    assert TEST_GENE_NAME in result


def test_all_protein_ids(ensembl100):
    result = ensembl100.all_protein_ids()
    assert isinstance(result, list)
    assert TEST_PROTEIN_ID in result


def test_all_transcript_ids(ensembl100):
    result = ensembl100.all_transcript_ids()
    assert isinstance(result, list)
    assert TEST_TRANSCRIPT_ID in result


def test_all_transcript_names(ensembl100):
    result = ensembl100.all_transcript_names()
    assert isinstance(result, list)
    assert TEST_TRANSCRIPT_NAME in result


# -------------------------------------------------------------------------------------------------
# get_<feature_symbol>(feature)
# -------------------------------------------------------------------------------------------------
def test_contig_ids_from_contig_id(ensembl100):
    result = ensembl100.contig_ids(TEST_CONTIG_ID)
    assert isinstance(result, list)
    assert TEST_CONTIG_ID in result


def test_contig_ids_from_contig_with_chr(ensembl100):
    ret_by_id = ensembl100.contig_ids("chr4")
    assert "4" in ret_by_id


def test_contig_ids_from_exon_id(ensembl100):
    result = ensembl100.contig_ids(TEST_EXON_ID)
    assert isinstance(result, list)
    assert TEST_CONTIG_ID in result


def test_contig_ids_from_gene(ensembl100):
    ret_by_id = ensembl100.contig_ids(TEST_GENE_ID)
    ret_by_name = ensembl100.contig_ids(TEST_GENE_NAME)
    assert TEST_CONTIG_ID in ret_by_id
    assert TEST_CONTIG_ID in ret_by_name
    assert ret_by_id == ret_by_name


def test_contig_ids_from_protein_id(ensembl100):
    result = ensembl100.contig_ids(TEST_PROTEIN_ID)
    assert isinstance(result, list)
    assert TEST_CONTIG_ID in result


def test_contig_ids_from_transcript(ensembl100):
    ret_by_id = ensembl100.contig_ids(TEST_TRANSCRIPT_ID)
    ret_by_name = ensembl100.contig_ids(TEST_TRANSCRIPT_NAME)
    assert TEST_CONTIG_ID in ret_by_id
    assert TEST_CONTIG_ID in ret_by_name
    assert ret_by_id == ret_by_name


def test_exon_ids_from_contig_id(ensembl100):
    result = ensembl100.exon_ids(TEST_CONTIG_ID)
    assert isinstance(result, list)
    assert TEST_EXON_ID in result


def test_exon_ids_from_exon_id(ensembl100):
    result = ensembl100.exon_ids(TEST_EXON_ID)
    assert isinstance(result, list)
    assert TEST_EXON_ID in result


def test_exon_ids_from_gene(ensembl100):
    ret_by_id = ensembl100.exon_ids(TEST_GENE_ID)
    ret_by_name = ensembl100.exon_ids(TEST_GENE_NAME)
    assert TEST_EXON_ID in ret_by_id
    assert TEST_EXON_ID in ret_by_name
    assert ret_by_id == ret_by_name


def test_exon_ids_from_protein_id(ensembl100):
    result = ensembl100.exon_ids(TEST_PROTEIN_ID)
    assert isinstance(result, list)
    assert TEST_EXON_ID in result


def test_exon_ids_from_transcript(ensembl100):
    ret_by_id = ensembl100.exon_ids(TEST_TRANSCRIPT_ID)
    ret_by_name = ensembl100.exon_ids(TEST_TRANSCRIPT_NAME)
    assert TEST_EXON_ID in ret_by_id
    assert TEST_EXON_ID in ret_by_name
    assert ret_by_id == ret_by_name


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
    result = ensembl100.gene_ids(TEST_CONTIG_ID)
    assert isinstance(result, list)
    assert TEST_GENE_ID in result


def test_gene_ids_from_exon_id(ensembl100):
    result = ensembl100.gene_ids(TEST_EXON_ID)
    assert isinstance(result, list)
    assert TEST_GENE_ID in result


def test_gene_ids_from_gene(ensembl100):
    ret_by_id = ensembl100.gene_ids(TEST_GENE_ID)
    ret_by_name = ensembl100.gene_ids(TEST_GENE_NAME)
    assert TEST_GENE_ID in ret_by_id
    assert TEST_GENE_ID in ret_by_name
    assert ret_by_id == ret_by_name


def test_gene_ids_from_protein_id(ensembl100):
    result = ensembl100.gene_ids(TEST_PROTEIN_ID)
    assert isinstance(result, list)
    assert TEST_GENE_ID in result


def test_gene_ids_from_transcript(ensembl100):
    ret_by_id = ensembl100.gene_ids(TEST_TRANSCRIPT_ID)
    ret_by_name = ensembl100.gene_ids(TEST_TRANSCRIPT_NAME)
    assert TEST_GENE_ID in ret_by_id
    assert TEST_GENE_ID in ret_by_name
    assert ret_by_id == ret_by_name


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
    result = ensembl100.gene_names(TEST_CONTIG_ID)
    assert isinstance(result, list)
    assert TEST_GENE_NAME in result


def test_gene_names_from_exon_id(ensembl100):
    result = ensembl100.gene_names(TEST_EXON_ID)
    assert isinstance(result, list)
    assert TEST_GENE_NAME in result


def test_gene_names_from_gene(ensembl100):
    ret_by_id = ensembl100.gene_names(TEST_GENE_ID)
    ret_by_name = ensembl100.gene_names(TEST_GENE_NAME)
    assert TEST_GENE_NAME in ret_by_id
    assert TEST_GENE_NAME in ret_by_name
    assert ret_by_id == ret_by_name


def test_gene_names_from_protein_id(ensembl100):
    result = ensembl100.gene_names(TEST_PROTEIN_ID)
    assert isinstance(result, list)
    assert TEST_GENE_NAME in result


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
    ret_by_id = ensembl100.gene_names(TEST_TRANSCRIPT_ID)
    ret_by_name = ensembl100.gene_names(TEST_TRANSCRIPT_NAME)
    assert TEST_GENE_NAME in ret_by_id
    assert TEST_GENE_NAME in ret_by_name
    assert ret_by_id == ret_by_name


def test_protein_ids_from_contig_id(ensembl100):
    result = ensembl100.protein_ids(TEST_CONTIG_ID)
    assert isinstance(result, list)
    assert TEST_PROTEIN_ID in result


def test_protein_ids_from_exon_id(ensembl100):
    result = ensembl100.protein_ids(TEST_EXON_ID)
    assert isinstance(result, list)
    assert TEST_PROTEIN_ID in result


def test_protein_ids_from_gene(ensembl100):
    ret_by_id = ensembl100.protein_ids(TEST_GENE_ID)
    ret_by_name = ensembl100.protein_ids(TEST_GENE_NAME)
    assert TEST_PROTEIN_ID in ret_by_id
    assert TEST_PROTEIN_ID in ret_by_name
    assert ret_by_id == ret_by_name


def test_protein_ids_from_protein_id(ensembl100):
    result = ensembl100.protein_ids(TEST_PROTEIN_ID)
    assert isinstance(result, list)
    assert TEST_PROTEIN_ID in result


def test_protein_ids_from_transcript(ensembl100):
    ret_by_id = ensembl100.protein_ids(TEST_TRANSCRIPT_ID)
    ret_by_name = ensembl100.protein_ids(TEST_TRANSCRIPT_NAME)
    assert TEST_PROTEIN_ID in ret_by_id
    assert TEST_PROTEIN_ID in ret_by_name
    assert ret_by_id == ret_by_name


def test_transcript_ids_from_contig_id(ensembl100):
    result = ensembl100.transcript_ids(TEST_CONTIG_ID)
    assert isinstance(result, list)
    assert TEST_TRANSCRIPT_ID in result


def test_transcript_ids_from_exon_id(ensembl100):
    result = ensembl100.transcript_ids(TEST_EXON_ID)
    assert isinstance(result, list)
    assert TEST_TRANSCRIPT_ID in result


def test_transcript_ids_from_gene(ensembl100):
    ret_by_id = ensembl100.transcript_ids(TEST_GENE_ID)
    ret_by_name = ensembl100.transcript_ids(TEST_GENE_NAME)
    assert TEST_TRANSCRIPT_ID in ret_by_id
    assert TEST_TRANSCRIPT_ID in ret_by_name
    assert ret_by_id == ret_by_name


def test_transcript_ids_from_protein_id(ensembl100):
    result = ensembl100.transcript_ids(TEST_PROTEIN_ID)
    assert isinstance(result, list)
    assert TEST_TRANSCRIPT_ID in result


def test_transcript_ids_from_transcript(ensembl100):
    ret_by_id = ensembl100.transcript_ids(TEST_TRANSCRIPT_ID)
    ret_by_name = ensembl100.transcript_ids(TEST_TRANSCRIPT_NAME)
    assert TEST_TRANSCRIPT_ID in ret_by_id
    assert TEST_TRANSCRIPT_ID in ret_by_name
    assert ret_by_id == ret_by_name


def test_transcript_names_from_contig_id(ensembl100):
    result = ensembl100.transcript_names(TEST_CONTIG_ID)
    assert isinstance(result, list)
    assert TEST_TRANSCRIPT_NAME in result


def test_transcript_names_from_exon_id(ensembl100):
    result = ensembl100.transcript_names(TEST_EXON_ID)
    assert isinstance(result, list)
    assert TEST_TRANSCRIPT_NAME in result


def test_transcript_names_from_gene(ensembl100):
    ret_by_id = ensembl100.transcript_names(TEST_GENE_ID)
    ret_by_name = ensembl100.transcript_names(TEST_GENE_NAME)
    assert TEST_TRANSCRIPT_NAME in ret_by_id
    assert TEST_TRANSCRIPT_NAME in ret_by_name
    assert ret_by_id == ret_by_name


def test_transcript_names_from_protein_id(ensembl100):
    result = ensembl100.transcript_names(TEST_PROTEIN_ID)
    assert isinstance(result, list)
    assert TEST_TRANSCRIPT_NAME in result


def test_transcript_names_from_transcript(ensembl100):
    ret_by_id = ensembl100.transcript_names(TEST_TRANSCRIPT_ID)
    ret_by_name = ensembl100.transcript_names(TEST_TRANSCRIPT_NAME)
    assert TEST_TRANSCRIPT_NAME in ret_by_id
    assert TEST_TRANSCRIPT_NAME in ret_by_name
    assert ret_by_id == ret_by_name


# -------------------------------------------------------------------------------------------------
# canonical transcript
# -------------------------------------------------------------------------------------------------
def test_is_canonical_transcript(ensembl100):
    assert ensembl100.is_canonical_transcript("ENST00000000233") is True
    assert ensembl100.is_canonical_transcript("spam") is False


# -------------------------------------------------------------------------------------------------
# is_<feature>(feature)
# -------------------------------------------------------------------------------------------------
def test_is_contig_false(ensembl100):
    assert ensembl100.is_contig("spam") is False


def test_is_contig_true(ensembl100):
    assert ensembl100.is_contig("4") is True


def test_is_contig_true_alias(ensembl100):
    assert ensembl100.is_contig("chr4") is True


def test_is_exon_false(ensembl100):
    assert ensembl100.is_exon("spam") is False


def test_is_exon_true(ensembl100):
    assert ensembl100.is_exon(TEST_EXON_ID) is True


def test_is_gene_false(ensembl100):
    assert ensembl100.is_gene("spam") is False


def test_is_gene_true_from_id(ensembl100):
    assert ensembl100.is_gene(TEST_GENE_ID) is True


def test_is_gene_true_from_name(ensembl100):
    assert ensembl100.is_gene(TEST_GENE_NAME) is True


def test_is_protein_false(ensembl100):
    assert ensembl100.is_protein("spam") is False


def test_is_protein_true(ensembl100):
    assert ensembl100.is_protein(TEST_PROTEIN_ID) is True


def test_is_transcript_false(ensembl100):
    assert ensembl100.is_transcript("spam") is False


def test_is_transcript_true_from_id(ensembl100):
    assert ensembl100.is_transcript(TEST_TRANSCRIPT_ID) is True


def test_is_transcript_true_from_name(ensembl100):
    assert ensembl100.is_transcript(TEST_TRANSCRIPT_NAME) is True


# -------------------------------------------------------------------------------------------------
# normalize_feature(feature, feature_type)
# -------------------------------------------------------------------------------------------------
def test_normalize_feature_contig_id(ensembl100):
    assert ensembl100.normalize_feature("4") == [("4", CONTIG_ID)]


def test_normalize_feature_exon_id(ensembl100):
    assert ensembl100.normalize_feature(TEST_EXON_ID) == [(TEST_EXON_ID, EXON_ID)]


def test_normalize_feature_gene_id(ensembl100):
    assert ensembl100.normalize_feature(TEST_GENE_ID) == [(TEST_GENE_ID, GENE_ID)]


def test_normalize_feature_gene_name(ensembl100):
    assert ensembl100.normalize_feature(TEST_GENE_NAME) == [(TEST_GENE_NAME, GENE_NAME)]


def test_normalize_feature_protein_id(ensembl100):
    assert ensembl100.normalize_feature(TEST_PROTEIN_ID) == [(TEST_PROTEIN_ID, PROTEIN_ID)]


def test_normalize_feature_refseq_transcript_id(ensembl100):
    assert ensembl100.normalize_feature("NM_000314.4") == [("ENST00000371953", TRANSCRIPT_ID)]


def test_normalize_feature_transcript_id(ensembl100):
    assert ensembl100.normalize_feature(TEST_TRANSCRIPT_ID) == [(TEST_TRANSCRIPT_ID, TRANSCRIPT_ID)]


def test_normalize_feature_transcript_name(ensembl100):
    assert ensembl100.normalize_feature(TEST_TRANSCRIPT_NAME) == [
        (TEST_TRANSCRIPT_NAME, TRANSCRIPT_NAME)
    ]


# -------------------------------------------------------------------------------------------------
# str()
# -------------------------------------------------------------------------------------------------
def test_cdna_to_str(ensembl100):
    position1 = CdnaPosition(
        _data=ensembl100,
        contig_id="5",
        start=3399,
        start_offset=0,
        end=3399,
        end_offset=0,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
        protein_id="ENSP00000309572",
    )
    position2 = CdnaPosition(
        _data=ensembl100,
        contig_id="5",
        start=1,
        start_offset=0,
        end=3399,
        end_offset=0,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
        protein_id="ENSP00000309572",
    )
    assert str(position1) == "ENST00000310581:c.3399"
    assert str(position2) == "ENST00000310581:c.1_3399"


def test_dna_to_str(ensembl100):
    position1 = DnaPosition(
        _data=ensembl100,
        contig_id="5",
        start=1234,
        start_offset=0,
        end=1234,
        end_offset=0,
        strand="-",
    )
    position2 = DnaPosition(
        _data=ensembl100, contig_id="5", start=1, start_offset=0, end=1234, end_offset=0, strand="-"
    )
    assert str(position1) == "5:g.1234"
    assert str(position2) == "5:g.1_1234"


def test_exon_to_str(ensembl100):
    position1 = ExonPosition(
        _data=ensembl100,
        contig_id="5",
        start=1,
        start_offset=0,
        end=1,
        end_offset=0,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
        exon_id="ENSE00003896691",
    )
    position2 = ExonPosition(
        _data=ensembl100,
        contig_id="5",
        start=1,
        start_offset=0,
        end=2,
        end_offset=0,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
        exon_id="ENSE00003896691",
    )
    assert str(position1) == "ENST00000310581:e.1"
    assert str(position2) == "ENST00000310581:e.1_2"


def test_protein_to_str(ensembl100):
    position1 = ProteinPosition(
        _data=ensembl100,
        contig_id="5",
        start=524,
        start_offset=0,
        end=524,
        end_offset=0,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
        protein_id="ENSP00000309572",
    )
    position2 = ProteinPosition(
        _data=ensembl100,
        contig_id="5",
        start=74,
        start_offset=0,
        end=524,
        end_offset=0,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
        protein_id="ENSP00000309572",
    )
    assert str(position1) == "ENSP00000309572:p.524"
    assert str(position2) == "ENSP00000309572:p.74_524"


def test_transcript_to_str(ensembl100):
    position1 = RnaPosition(
        _data=ensembl100,
        contig_id="5",
        start=4039,
        start_offset=0,
        end=4039,
        end_offset=0,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
    )
    position2 = RnaPosition(
        _data=ensembl100,
        contig_id="5",
        start=1,
        start_offset=0,
        end=4039,
        end_offset=0,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
    )
    assert str(position1) == "ENST00000310581:r.4039"
    assert str(position2) == "ENST00000310581:r.1_4039"


# -------------------------------------------------------------------------------------------------
# get_<feature>(feature)
# -------------------------------------------------------------------------------------------------
@pytest.fixture
def test_cdna(ensembl100):
    return CdnaPosition(
        _data=ensembl100,
        contig_id="5",
        start=1,
        start_offset=0,
        end=3399,
        end_offset=0,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
        protein_id="ENSP00000309572",
    )


@pytest.mark.parametrize(
    "feature,num_results",
    [
        ("5", 4),
        ("ENSE00003896691", 1),
        ("ENSG00000164362", 4),
        ("ENSP00000309572", 1),
        ("ENST00000310581", 1),
        ("TERT", 4),
        ("TERT-201", 1),
    ],
)
def test_get_cdna(ensembl100, test_cdna, feature, num_results):
    results = ensembl100.get_cdna(feature)
    assert len(results) == num_results
    assert test_cdna in results


@pytest.fixture
def test_dna(ensembl100):
    return DnaPosition(
        _data=ensembl100,
        contig_id="5",
        start=1,
        start_offset=0,
        end=181538259,
        end_offset=0,
        strand="-",
    )


@pytest.mark.parametrize(
    "feature,num_results",
    [
        ("5", 1),
        ("ENSE00003896691", 1),
        ("ENSG00000164362", 1),
        ("ENSP00000309572", 1),
        ("ENST00000310581", 1),
        ("TERT", 1),
        ("TERT-201", 1),
    ],
)
def test_get_dna(ensembl100, test_dna, feature, num_results):
    results = ensembl100.get_dna(feature)
    assert len(results) == num_results
    assert test_dna in results


@pytest.fixture
def test_exon(ensembl100):
    return ExonPosition(
        _data=ensembl100,
        contig_id="5",
        start=1,
        start_offset=0,
        end=1,
        end_offset=0,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
        exon_id="ENSE00003896691",
    )


@pytest.mark.parametrize(
    "feature,num_results",
    [
        ("5", 79),
        ("ENSE00003896691", 1),
        ("ENSG00000164362", 79),
        ("ENSP00000309572", 33),
        ("ENST00000310581", 33),
        ("TERT", 79),
        ("TERT-201", 33),
    ],
)
def test_get_exon(ensembl100, test_exon, feature, num_results):
    results = ensembl100.get_exons(feature)
    assert len(results) == num_results
    assert test_exon in results


@pytest.fixture
def test_gene(ensembl100):
    return DnaPosition(
        _data=ensembl100,
        contig_id="5",
        start=1253147,
        start_offset=0,
        end=1295068,
        end_offset=0,
        strand="-",
    )


@pytest.mark.parametrize(
    "feature,num_results",
    [
        ("5", 1),
        ("ENSE00003896691", 1),
        ("ENSG00000164362", 1),
        ("ENSP00000309572", 1),
        ("ENST00000310581", 1),
        ("TERT", 1),
        ("TERT-201", 1),
    ],
)
def test_get_gene(ensembl100, test_gene, feature, num_results):
    results = ensembl100.get_genes(feature)
    assert len(results) == num_results
    assert test_gene in results


@pytest.fixture
def test_transcript(ensembl100):
    return RnaPosition(
        _data=ensembl100,
        contig_id="5",
        start=1,
        start_offset=0,
        end=4039,
        end_offset=0,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
    )


@pytest.mark.parametrize(
    "feature,num_results",
    [
        ("5", 7),
        ("ENSE00003896691", 1),
        ("ENSG00000164362", 7),
        ("ENSP00000309572", 1),
        ("ENST00000310581", 1),
        ("TERT", 7),
        ("TERT-201", 1),
    ],
)
def test_get_transcript(ensembl100, test_transcript, feature, num_results):
    results = ensembl100.get_transcripts(feature)
    assert len(results) == num_results
    assert test_transcript in results
