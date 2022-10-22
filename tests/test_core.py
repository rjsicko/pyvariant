import pytest

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
    DnaPosition,
    EnsemblRelease,
    ExonPosition,
    ProteinPosition,
    RnaPosition,
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

# features are from Ensembl 100
TEST_EXON_ID = "ENSE00003826864"
TEST_CONTIG_ID = "16"
TEST_GENE_ID = "ENSG00000149925"
TEST_GENE_NAME = "ALDOA"
TEST_PROTEIN_ID = "ENSP00000494188"  # version 2
TEST_TRANSCRIPT_ID = "ENST00000643777"  # version 3
TEST_TRANSCRIPT_NAME = "ALDOA-219"


@pytest.fixture
def ensembl69():
    return EnsemblRelease(
        species="homo_sapiens",
        release=69,
        cache_dir=CACHE_DIR,
        canonical_transcript="",
        contig_alias="",
        exon_alias="",
        gene_alias="",
        protein_alias="",
        transcript_alias="",
    )


@pytest.fixture
def ensembl100():
    return EnsemblRelease(
        species="homo_sapiens",
        release=100,
        cache_dir=CACHE_DIR,
        canonical_transcript=CANONICAL_TRANSCRIPT,
        contig_alias=CONTIG_ALIAS,
        exon_alias=EXON_ALIAS,
        gene_alias=GENE_ALIAS,
        protein_alias=PROTEIN_ALIAS,
        transcript_alias=TRANSCRIPT_ALIAS,
    )


# -------------------------------------------------------------------------------------------------
# cache logic
# -------------------------------------------------------------------------------------------------
def test_unique_instances():
    a = EnsemblRelease(species="homo_sapiens", release=100, cache_dir=CACHE_DIR)
    b = EnsemblRelease(species="homo_sapiens", release=100, cache_dir=CACHE_DIR)
    c = EnsemblRelease(species="homo_sapiens", release=69, cache_dir=CACHE_DIR)
    assert a is b
    assert a is not c


# ---------------------------------------------------------------------------------------------
# <feature_symbol>()
# ---------------------------------------------------------------------------------------------
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


# ---------------------------------------------------------------------------------------------
# get_<feature_symbol>(feature)
# ---------------------------------------------------------------------------------------------
def test_contig_ids_from_contig_id(ensembl100):
    result = ensembl100.contig_ids(TEST_CONTIG_ID)
    assert isinstance(result, list)
    assert TEST_CONTIG_ID in result


def test_contig_ids_from_contig_with_chr(ensembl100):
    ret_by_id = ensembl100.contig_ids("chr12")
    assert "12" in ret_by_id


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


# ---------------------------------------------------------------------------------------------
# canonical transcript
# ---------------------------------------------------------------------------------------------
def test_is_canonical_transcript(ensembl100):
    assert ensembl100.is_canonical_transcript("ENST00000000233") is True
    assert ensembl100.is_canonical_transcript("spam") is False


# ---------------------------------------------------------------------------------------------
# is_<feature>(feature)
# ---------------------------------------------------------------------------------------------
def test_is_contig_false(ensembl100):
    assert ensembl100.is_contig("spam") is False


def test_is_contig_true(ensembl100):
    assert ensembl100.is_contig("12") is True


def test_is_contig_true_alias(ensembl100):
    assert ensembl100.is_contig("chr12") is True


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


# ---------------------------------------------------------------------------------------------
# normalize_feature(feature, feature_type)
# ---------------------------------------------------------------------------------------------
def test_normalize_feature_contig_id(ensembl100):
    assert ensembl100.normalize_feature("17") == [("17", CONTIG_ID)]


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


# ---------------------------------------------------------------------------------------------
# get_<feature>(feature)
# ---------------------------------------------------------------------------------------------
@pytest.fixture
def test_cdna(ensembl100):
    return CdnaPosition(
        _data=ensembl100,
        contig_id="5",
        start=1,
        end=3399,
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
    return DnaPosition(_data=ensembl100, contig_id="5", start=1, end=181538259, strand="-")


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
        end=1,
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
    return DnaPosition(_data=ensembl100, contig_id="5", start=1253147, end=1295068, strand="-")


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
        end=4039,
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


# -------------------------------------------------------------------------------------------------
# <feature>_to_<feature>
# -------------------------------------------------------------------------------------------------
@pytest.fixture
def cdna(ensembl100):
    return CdnaPosition(
        _data=ensembl100,
        contig_id="5",
        start=220,
        end=1572,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
        protein_id="ENSP00000309572",
    )


@pytest.fixture
def cdna_exon_end(ensembl100):
    return CdnaPosition(
        _data=ensembl100,
        contig_id="5",
        start=220,
        end=1573,  # codon breaks across exon boundary
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
        protein_id="ENSP00000309572",
    )


@pytest.fixture
def dna(ensembl100):
    return DnaPosition(_data=ensembl100, contig_id="5", start=1293314, end=1294666, strand="-")


@pytest.fixture
def dna_exon_end(ensembl100):
    return DnaPosition(_data=ensembl100, contig_id="5", start=1293313, end=1294666, strand="-")


@pytest.fixture
def exon(ensembl100):
    return ExonPosition(
        _data=ensembl100,
        contig_id="5",
        start=2,
        end=2,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
        exon_id="ENSE00001197112",
    )


@pytest.fixture
def protein(ensembl100):
    return ProteinPosition(
        _data=ensembl100,
        contig_id="5",
        start=74,
        end=524,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
        protein_id="ENSP00000309572",
    )


@pytest.fixture
def protein_exon_end(ensembl100):
    return ProteinPosition(
        _data=ensembl100,
        contig_id="5",
        start=74,
        end=525,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
        protein_id="ENSP00000309572",
    )


@pytest.fixture
def rna(ensembl100):
    return RnaPosition(
        _data=ensembl100,
        contig_id="5",
        start=299,
        end=1651,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
    )


@pytest.fixture
def rna_exon_end(ensembl100):
    return RnaPosition(
        _data=ensembl100,
        contig_id="5",
        start=299,
        end=1652,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
    )


feature_list = [
    "5",
    "ENSE00001197112",
    "ENSG00000164362",
    "ENSP00000309572",
    "ENST00000310581",
    "TERT",
    "TERT-201",
]


@pytest.mark.parametrize("feature", feature_list)
def test_cdna_to_cdna(ensembl100, cdna, feature):
    results = ensembl100.cdna_to_cdna(feature, cdna.start, cdna.end)
    assert cdna in results


@pytest.mark.parametrize("feature", feature_list)
def test_cdna_to_dna(ensembl100, cdna, dna, feature):
    results = ensembl100.cdna_to_dna(feature, cdna.start, cdna.end)
    assert dna in results


@pytest.mark.parametrize("feature", feature_list)
def test_cdna_to_exon(ensembl100, cdna_exon_end, exon, feature):
    results = ensembl100.cdna_to_exon(feature, cdna_exon_end.start, cdna_exon_end.end)
    assert exon in results


@pytest.mark.parametrize("feature", feature_list)
def test_cdna_to_protein(ensembl100, cdna, protein, feature):
    results = ensembl100.cdna_to_protein(feature, cdna.start, cdna.end)
    assert protein in results


@pytest.mark.parametrize("feature", feature_list)
def test_cdna_to_rna(ensembl100, cdna, rna, feature):
    results = ensembl100.cdna_to_rna(feature, cdna.start, cdna.end)
    assert rna in results


@pytest.mark.parametrize("feature", feature_list)
def test_dna_to_cdna(ensembl100, dna, cdna, feature):
    results = ensembl100.dna_to_cdna(feature, dna.start, dna.end)
    assert cdna in results


@pytest.mark.parametrize("feature", feature_list)
def test_dna_to_dna(ensembl100, dna, feature):
    results = ensembl100.dna_to_dna(feature, dna.start, dna.end)
    assert dna in results


@pytest.mark.parametrize("feature", feature_list)
def test_dna_to_exon(ensembl100, dna, exon, feature):
    results = ensembl100.dna_to_exon(feature, dna.start, dna.end)
    assert exon in results


@pytest.mark.parametrize("feature", feature_list)
def test_dna_to_protein(ensembl100, dna, protein, feature):
    results = ensembl100.dna_to_protein(feature, dna.start, dna.end)
    assert protein in results


@pytest.mark.parametrize("feature", feature_list)
def test_dna_to_rna(ensembl100, dna, rna, feature):
    results = ensembl100.dna_to_rna(feature, dna.start, dna.end)
    assert rna in results


@pytest.mark.parametrize("feature", feature_list)
def test_exon_to_cdna(ensembl100, exon, cdna_exon_end, feature):
    results = ensembl100.exon_to_cdna(feature, exon.start, exon.end)
    assert cdna_exon_end in results


@pytest.mark.parametrize("feature", feature_list)
def test_exon_to_dna(ensembl100, exon, dna_exon_end, feature):
    results = ensembl100.exon_to_dna(feature, exon.start, exon.end)
    assert dna_exon_end in results


@pytest.mark.parametrize("feature", feature_list)
def test_exon_to_exon(ensembl100, exon, feature):
    results = ensembl100.exon_to_exon(feature, exon.start, exon.end)
    assert exon in results


@pytest.mark.parametrize("feature", feature_list)
def test_exon_to_protein(ensembl100, exon, protein_exon_end, feature):
    results = ensembl100.exon_to_protein(feature, exon.start, exon.end)
    assert protein_exon_end in results


@pytest.mark.parametrize("feature", feature_list)
def test_exon_to_rna(ensembl100, exon, rna_exon_end, feature):
    results = ensembl100.exon_to_rna(feature, exon.start, exon.end)
    assert rna_exon_end in results


@pytest.mark.parametrize("feature", feature_list)
def test_protein_to_cdna(ensembl100, protein, cdna, feature):
    results = ensembl100.protein_to_cdna(feature, protein.start, protein.end)
    assert cdna in results


@pytest.mark.parametrize("feature", feature_list)
def test_protein_to_dna(ensembl100, protein, dna, feature):
    results = ensembl100.protein_to_dna(feature, protein.start, protein.end)
    assert dna in results


@pytest.mark.parametrize("feature", feature_list)
def test_protein_to_exon(ensembl100, protein, exon, feature):
    results = ensembl100.protein_to_exon(feature, protein.start, protein.end)
    assert exon in results


@pytest.mark.parametrize("feature", feature_list)
def test_protein_to_protein(ensembl100, protein, feature):
    results = ensembl100.protein_to_protein(feature, protein.start, protein.end)
    assert protein in results


@pytest.mark.parametrize("feature", feature_list)
def test_protein_to_rna(ensembl100, protein, rna, feature):
    results = ensembl100.protein_to_rna(feature, protein.start, protein.end)
    assert rna in results


@pytest.mark.parametrize("feature", feature_list)
def test_rna_to_cdna(ensembl100, rna, cdna, feature):
    results = ensembl100.rna_to_cdna(feature, rna.start, rna.end)
    assert cdna in results


@pytest.mark.parametrize("feature", feature_list)
def test_rna_to_dna(ensembl100, rna, dna, feature):
    results = ensembl100.rna_to_dna(feature, rna.start, rna.end)
    assert dna in results


@pytest.mark.parametrize("feature", feature_list)
def test_rna_to_exon(ensembl100, rna_exon_end, exon, feature):
    results = ensembl100.rna_to_exon(feature, rna_exon_end.start, rna_exon_end.end)
    assert exon in results


@pytest.mark.parametrize("feature", feature_list)
def test_rna_to_protein(ensembl100, rna, protein, feature):
    results = ensembl100.rna_to_protein(feature, rna.start, rna.end)
    assert protein in results


@pytest.mark.parametrize("feature", feature_list)
def test_rna_to_rna(ensembl100, rna, feature):
    results = ensembl100.rna_to_rna(feature, rna.start, rna.end)
    assert rna in results


# -------------------------------------------------------------------------------------------------
# position classes
# -------------------------------------------------------------------------------------------------
def test_cdna_to_cdna_negative_strand(ensembl100):
    position = CdnaPosition(
        _data=ensembl100,
        contig_id="5",
        start=124,
        end=126,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
        protein_id="ENSP00000309572",
    )
    assert position.to_cdna() == [position]


def test_cdna_to_cdna_negative_strand_across_exon_boundary(ensembl100):
    position = CdnaPosition(
        _data=ensembl100,
        contig_id="5",
        start=1573,
        end=1575,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
        protein_id="ENSP00000309572",
    )
    assert position.to_cdna() == [position]


def test_cdna_to_cdna_negative_strand_cDNA_protein_end(ensembl100):
    position = CdnaPosition(
        _data=ensembl100,
        contig_id="5",
        start=3397,
        end=3399,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
        protein_id="ENSP00000309572",
    )
    assert position.to_cdna() == [position]


def test_cdna_to_cdna_negative_strand_cDNA_protein_start(ensembl100):
    position = CdnaPosition(
        _data=ensembl100,
        contig_id="5",
        start=1,
        end=3,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
        protein_id="ENSP00000309572",
    )
    assert position.to_cdna() == [position]


def test_cdna_to_cdna_positive_strand(ensembl100):
    position = CdnaPosition(
        _data=ensembl100,
        contig_id="13",
        start=8566,
        end=8568,
        strand="+",
        gene_id="ENSG00000139618",
        gene_name="BRCA2",
        transcript_id="ENST00000380152",
        transcript_name="BRCA2-201",
        protein_id="ENSP00000369497",
    )
    assert position.to_cdna() == [position]


def test_cdna_to_cdna_positive_strand_across_exon_boundary(ensembl100):
    position = CdnaPosition(
        _data=ensembl100,
        contig_id="4",
        start=67,
        end=69,
        strand="+",
        gene_id="ENSG00000157404",
        gene_name="KIT",
        transcript_id="ENST00000288135",
        transcript_name="KIT-201",
        protein_id="ENSP00000288135",
    )
    assert position.to_cdna() == [position]


def test_cdna_to_cdna_positive_strand_cDNA_protein_end(ensembl100):
    position = CdnaPosition(
        _data=ensembl100,
        contig_id="13",
        start=10255,
        end=10257,
        strand="+",
        gene_id="ENSG00000139618",
        gene_name="BRCA2",
        transcript_id="ENST00000380152",
        transcript_name="BRCA2-201",
        protein_id="ENSP00000369497",
    )
    assert position.to_cdna() == [position]


def test_cdna_to_cdna_positive_strand_cDNA_protein_start(ensembl100):
    position = CdnaPosition(
        _data=ensembl100,
        contig_id="13",
        start=1,
        end=3,
        strand="+",
        gene_id="ENSG00000139618",
        gene_name="BRCA2",
        transcript_id="ENST00000380152",
        transcript_name="BRCA2-201",
        protein_id="ENSP00000369497",
    )
    assert position.to_cdna() == [position]


def test_cdna_to_cdna_positive_strand_overlapping_genes_different_strands(ensembl100):
    position = CdnaPosition(
        _data=ensembl100,
        contig_id="6",
        start=580,
        end=582,
        strand="+",
        gene_id="ENSG00000137310",
        gene_name="TCF19",
        transcript_id="ENST00000542218",
        transcript_name="TCF19-204",
        protein_id="ENSP00000439397",
    )
    assert position.to_cdna() == [position]


def test_cdna_to_dna_negative_strand(ensembl100):
    position = CdnaPosition(
        _data=ensembl100,
        contig_id="5",
        start=124,
        end=126,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
        protein_id="ENSP00000309572",
    )
    expected = DnaPosition(_data=ensembl100, contig_id="5", start=1294864, end=1294866, strand="-")
    assert position.to_dna() == [expected]


def test_cdna_to_dna_negative_strand_across_exon_boundary(ensembl100):
    position = CdnaPosition(
        _data=ensembl100,
        contig_id="5",
        start=1573,
        end=1575,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
        protein_id="ENSP00000309572",
    )
    expected = DnaPosition(_data=ensembl100, contig_id="5", start=1282623, end=1293313, strand="-")
    assert position.to_dna() == [expected]


def test_cdna_to_dna_negative_strand_cDNA_protein_end(ensembl100):
    position = CdnaPosition(
        _data=ensembl100,
        contig_id="5",
        start=3397,
        end=3399,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
        protein_id="ENSP00000309572",
    )
    expected = DnaPosition(_data=ensembl100, contig_id="5", start=1253728, end=1253730, strand="-")
    assert position.to_dna() == [expected]


def test_cdna_to_dna_negative_strand_cDNA_protein_start(ensembl100):
    position = CdnaPosition(
        _data=ensembl100,
        contig_id="5",
        start=1,
        end=3,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
        protein_id="ENSP00000309572",
    )
    expected = DnaPosition(_data=ensembl100, contig_id="5", start=1294987, end=1294989, strand="-")
    assert position.to_dna() == [expected]


def test_cdna_to_dna_positive_strand(ensembl100):
    position = CdnaPosition(
        _data=ensembl100,
        contig_id="13",
        start=8566,
        end=8568,
        strand="+",
        gene_id="ENSG00000139618",
        gene_name="BRCA2",
        transcript_id="ENST00000380152",
        transcript_name="BRCA2-201",
        protein_id="ENSP00000369497",
    )
    expected = DnaPosition(
        _data=ensembl100, contig_id="13", start=32371034, end=32371036, strand="+"
    )
    assert position.to_dna() == [expected]


def test_cdna_to_dna_positive_strand_across_exon_boundary(ensembl100):
    position = CdnaPosition(
        _data=ensembl100,
        contig_id="4",
        start=67,
        end=69,
        strand="+",
        gene_id="ENSG00000157404",
        gene_name="KIT",
        transcript_id="ENST00000288135",
        transcript_name="KIT-201",
        protein_id="ENSP00000288135",
    )
    expected = DnaPosition(
        _data=ensembl100, contig_id="4", start=54658081, end=54695513, strand="+"
    )
    assert position.to_dna() == [expected]


def test_cdna_to_dna_positive_strand_cDNA_protein_end(ensembl100):
    position = CdnaPosition(
        _data=ensembl100,
        contig_id="13",
        start=10255,
        end=10257,
        strand="+",
        gene_id="ENSG00000139618",
        gene_name="BRCA2",
        transcript_id="ENST00000380152",
        transcript_name="BRCA2-201",
        protein_id="ENSP00000369497",
    )
    expected = DnaPosition(
        _data=ensembl100, contig_id="13", start=32398768, end=32398770, strand="+"
    )
    assert position.to_dna() == [expected]


def test_cdna_to_dna_positive_strand_cDNA_protein_start(ensembl100):
    position = CdnaPosition(
        _data=ensembl100,
        contig_id="13",
        start=1,
        end=3,
        strand="+",
        gene_id="ENSG00000139618",
        gene_name="BRCA2",
        transcript_id="ENST00000380152",
        transcript_name="BRCA2-201",
        protein_id="ENSP00000369497",
    )
    expected = DnaPosition(
        _data=ensembl100, contig_id="13", start=32316461, end=32316463, strand="+"
    )
    assert position.to_dna() == [expected]


def test_cdna_to_dna_positive_strand_overlapping_genes_different_strands(ensembl100):
    position = CdnaPosition(
        _data=ensembl100,
        contig_id="6",
        start=580,
        end=582,
        strand="+",
        gene_id="ENSG00000137310",
        gene_name="TCF19",
        transcript_id="ENST00000542218",
        transcript_name="TCF19-204",
        protein_id="ENSP00000439397",
    )
    expected = DnaPosition(
        _data=ensembl100, contig_id="6", start=31166302, end=31166304, strand="+"
    )
    assert position.to_dna() == [expected]


def test_cdna_to_exon_negative_strand(ensembl100):
    position = CdnaPosition(
        _data=ensembl100,
        contig_id="5",
        start=124,
        end=126,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
        protein_id="ENSP00000309572",
    )
    expected = ExonPosition(
        _data=ensembl100,
        contig_id="5",
        start=1,
        end=1,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
        exon_id="ENSE00003896691",
    )
    assert position.to_exon() == [expected]


def test_cdna_to_exon_negative_strand_cDNA_protein_end(ensembl100):
    position = CdnaPosition(
        _data=ensembl100,
        contig_id="5",
        start=3397,
        end=3399,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
        protein_id="ENSP00000309572",
    )
    expected = ExonPosition(
        _data=ensembl100,
        contig_id="5",
        start=16,
        end=16,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
        exon_id="ENSE00001863787",
    )
    assert position.to_exon() == [expected]


def test_cdna_to_exon_negative_strand_cDNA_protein_start(ensembl100):
    position = CdnaPosition(
        _data=ensembl100,
        contig_id="5",
        start=1,
        end=3,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
        protein_id="ENSP00000309572",
    )
    expected = ExonPosition(
        _data=ensembl100,
        contig_id="5",
        start=1,
        end=1,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
        exon_id="ENSE00003896691",
    )
    assert position.to_exon() == [expected]


def test_cdna_to_exon_positive_strand(ensembl100):
    position = CdnaPosition(
        _data=ensembl100,
        contig_id="13",
        start=8566,
        end=8568,
        strand="+",
        gene_id="ENSG00000139618",
        gene_name="BRCA2",
        transcript_id="ENST00000380152",
        transcript_name="BRCA2-201",
        protein_id="ENSP00000369497",
    )
    expected = ExonPosition(
        _data=ensembl100,
        contig_id="13",
        start=20,
        end=20,
        strand="+",
        gene_id="ENSG00000139618",
        gene_name="BRCA2",
        transcript_id="ENST00000380152",
        transcript_name="BRCA2-201",
        exon_id="ENSE00000939180",
    )
    assert position.to_exon() == [expected]


def test_cdna_to_exon_positive_strand_cDNA_protein_end(ensembl100):
    position = CdnaPosition(
        _data=ensembl100,
        contig_id="13",
        start=10255,
        end=10257,
        strand="+",
        gene_id="ENSG00000139618",
        gene_name="BRCA2",
        transcript_id="ENST00000380152",
        transcript_name="BRCA2-201",
        protein_id="ENSP00000369497",
    )
    expected = ExonPosition(
        _data=ensembl100,
        contig_id="13",
        start=27,
        end=27,
        strand="+",
        gene_id="ENSG00000139618",
        gene_name="BRCA2",
        transcript_id="ENST00000380152",
        transcript_name="BRCA2-201",
        exon_id="ENSE00003717596",
    )
    assert position.to_exon() == [expected]


def test_cdna_to_exon_positive_strand_cDNA_protein_start(ensembl100):
    position = CdnaPosition(
        _data=ensembl100,
        contig_id="13",
        start=1,
        end=3,
        strand="+",
        gene_id="ENSG00000139618",
        gene_name="BRCA2",
        transcript_id="ENST00000380152",
        transcript_name="BRCA2-201",
        protein_id="ENSP00000369497",
    )
    expected = ExonPosition(
        _data=ensembl100,
        contig_id="13",
        start=2,
        end=2,
        strand="+",
        gene_id="ENSG00000139618",
        gene_name="BRCA2",
        transcript_id="ENST00000380152",
        transcript_name="BRCA2-201",
        exon_id="ENSE00001484009",
    )
    assert position.to_exon() == [expected]


def test_cdna_to_exon_positive_strand_overlapping_genes_different_strands(ensembl100):
    position = CdnaPosition(
        _data=ensembl100,
        contig_id="6",
        start=580,
        end=582,
        strand="+",
        gene_id="ENSG00000137310",
        gene_name="TCF19",
        transcript_id="ENST00000542218",
        transcript_name="TCF19-204",
        protein_id="ENSP00000439397",
    )
    expected = ExonPosition(
        _data=ensembl100,
        contig_id="6",
        start=2,
        end=2,
        strand="+",
        gene_id="ENSG00000137310",
        gene_name="TCF19",
        transcript_id="ENST00000542218",
        transcript_name="TCF19-204",
        exon_id="ENSE00002283659",
    )
    assert position.to_exon() == [expected]


def test_cdna_to_protein_negative_strand(ensembl100):
    position = CdnaPosition(
        _data=ensembl100,
        contig_id="5",
        start=124,
        end=126,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
        protein_id="ENSP00000309572",
    )
    expected = ProteinPosition(
        _data=ensembl100,
        contig_id="5",
        start=42,
        end=42,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
        protein_id="ENSP00000309572",
    )
    assert position.to_protein() == [expected]


def test_cdna_to_protein_negative_strand_across_exon_boundary(ensembl100):
    position = CdnaPosition(
        _data=ensembl100,
        contig_id="5",
        start=1573,
        end=1575,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
        protein_id="ENSP00000309572",
    )
    expected = ProteinPosition(
        _data=ensembl100,
        contig_id="5",
        start=525,
        end=525,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
        protein_id="ENSP00000309572",
    )
    assert position.to_protein() == [expected]


def test_cdna_to_protein_negative_strand_cDNA_protein_end(ensembl100):
    position = CdnaPosition(
        _data=ensembl100,
        contig_id="5",
        start=3394,
        end=3396,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
        protein_id="ENSP00000309572",
    )
    expected = ProteinPosition(
        _data=ensembl100,
        contig_id="5",
        start=1132,
        end=1132,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
        protein_id="ENSP00000309572",
    )
    assert position.to_protein() == [expected]


def test_cdna_to_protein_negative_strand_stop_codon(ensembl100):
    position = CdnaPosition(
        _data=ensembl100,
        contig_id="5",
        start=3397,
        end=3399,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
        protein_id="ENSP00000309572",
    )
    assert position.to_protein() == []


def test_cdna_to_protein_negative_strand_cDNA_protein_start(ensembl100):
    position = CdnaPosition(
        _data=ensembl100,
        contig_id="5",
        start=1,
        end=3,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
        protein_id="ENSP00000309572",
    )
    expected = ProteinPosition(
        _data=ensembl100,
        contig_id="5",
        start=1,
        end=1,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
        protein_id="ENSP00000309572",
    )
    assert position.to_protein() == [expected]


def test_cdna_to_protein_positive_strand(ensembl100):
    position = CdnaPosition(
        _data=ensembl100,
        contig_id="13",
        start=8566,
        end=8568,
        strand="+",
        gene_id="ENSG00000139618",
        gene_name="BRCA2",
        transcript_id="ENST00000380152",
        transcript_name="BRCA2-201",
        protein_id="ENSP00000369497",
    )
    expected = ProteinPosition(
        _data=ensembl100,
        contig_id="13",
        start=2856,
        end=2856,
        strand="+",
        gene_id="ENSG00000139618",
        gene_name="BRCA2",
        transcript_id="ENST00000380152",
        transcript_name="BRCA2-201",
        protein_id="ENSP00000369497",
    )
    assert position.to_protein() == [expected]


def test_cdna_to_protein_positive_strand_across_exon_boundary(ensembl100):
    position = CdnaPosition(
        _data=ensembl100,
        contig_id="4",
        start=67,
        end=69,
        strand="+",
        gene_id="ENSG00000157404",
        gene_name="KIT",
        transcript_id="ENST00000288135",
        transcript_name="KIT-201",
        protein_id="ENSP00000288135",
    )
    expected = ProteinPosition(
        _data=ensembl100,
        contig_id="4",
        start=23,
        end=23,
        strand="+",
        gene_id="ENSG00000157404",
        gene_name="KIT",
        transcript_id="ENST00000288135",
        transcript_name="KIT-201",
        protein_id="ENSP00000288135",
    )
    assert position.to_protein() == [expected]


def test_cdna_to_protein_positive_strand_cDNA_protein_end(ensembl100):
    position = CdnaPosition(
        _data=ensembl100,
        contig_id="13",
        start=10252,
        end=10254,
        strand="+",
        gene_id="ENSG00000139618",
        gene_name="BRCA2",
        transcript_id="ENST00000380152",
        transcript_name="BRCA2-201",
        protein_id="ENSP00000369497",
    )
    expected = ProteinPosition(
        _data=ensembl100,
        contig_id="13",
        start=3418,
        end=3418,
        strand="+",
        gene_id="ENSG00000139618",
        gene_name="BRCA2",
        transcript_id="ENST00000380152",
        transcript_name="BRCA2-201",
        protein_id="ENSP00000369497",
    )
    assert position.to_protein() == [expected]


def test_cdna_to_protein_positive_strand_stop_codon(ensembl100):
    position = CdnaPosition(
        _data=ensembl100,
        contig_id="13",
        start=10255,
        end=10257,
        strand="+",
        gene_id="ENSG00000139618",
        gene_name="BRCA2",
        transcript_id="ENST00000380152",
        transcript_name="BRCA2-201",
        protein_id="ENSP00000369497",
    )
    assert position.to_protein() == []


def test_cdna_to_protein_positive_strand_cDNA_protein_start(ensembl100):
    position = CdnaPosition(
        _data=ensembl100,
        contig_id="13",
        start=1,
        end=3,
        strand="+",
        gene_id="ENSG00000139618",
        gene_name="BRCA2",
        transcript_id="ENST00000380152",
        transcript_name="BRCA2-201",
        protein_id="ENSP00000369497",
    )
    expected = ProteinPosition(
        _data=ensembl100,
        contig_id="13",
        start=1,
        end=1,
        strand="+",
        gene_id="ENSG00000139618",
        gene_name="BRCA2",
        transcript_id="ENST00000380152",
        transcript_name="BRCA2-201",
        protein_id="ENSP00000369497",
    )
    assert position.to_protein() == [expected]


def test_cdna_to_protein_positive_strand_overlapping_genes_different_strands(ensembl100):
    position = CdnaPosition(
        _data=ensembl100,
        contig_id="6",
        start=580,
        end=582,
        strand="+",
        gene_id="ENSG00000137310",
        gene_name="TCF19",
        transcript_id="ENST00000542218",
        transcript_name="TCF19-204",
        protein_id="ENSP00000439397",
    )
    expected = ProteinPosition(
        _data=ensembl100,
        contig_id="6",
        start=194,
        end=194,
        strand="+",
        gene_id="ENSG00000137310",
        gene_name="TCF19",
        transcript_id="ENST00000542218",
        transcript_name="TCF19-204",
        protein_id="ENSP00000439397",
    )
    assert position.to_protein() == [expected]


def test_cdna_to_rna_negative_strand(ensembl100):
    position = CdnaPosition(
        _data=ensembl100,
        contig_id="5",
        start=124,
        end=126,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
        protein_id="ENSP00000309572",
    )
    expected = RnaPosition(
        _data=ensembl100,
        contig_id="5",
        start=203,
        end=205,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
    )
    assert position.to_rna() == [expected]


def test_cdna_to_rna_negative_strand_across_exon_boundary(ensembl100):
    position = CdnaPosition(
        _data=ensembl100,
        contig_id="5",
        start=1573,
        end=1575,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
        protein_id="ENSP00000309572",
    )
    expected = RnaPosition(
        _data=ensembl100,
        contig_id="5",
        start=1652,
        end=1654,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
    )
    assert position.to_rna() == [expected]


def test_cdna_to_rna_negative_strand_cDNA_protein_end(ensembl100):
    position = CdnaPosition(
        _data=ensembl100,
        contig_id="5",
        start=3397,
        end=3399,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
        protein_id="ENSP00000309572",
    )
    expected = RnaPosition(
        _data=ensembl100,
        contig_id="5",
        start=3476,
        end=3478,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
    )
    assert position.to_rna() == [expected]


def test_cdna_to_rna_negative_strand_cDNA_protein_start(ensembl100):
    position = CdnaPosition(
        _data=ensembl100,
        contig_id="5",
        start=1,
        end=3,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
        protein_id="ENSP00000309572",
    )
    expected = RnaPosition(
        _data=ensembl100,
        contig_id="5",
        start=80,
        end=82,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
    )
    assert position.to_rna() == [expected]


def test_cdna_to_rna_positive_strand(ensembl100):
    position = CdnaPosition(
        _data=ensembl100,
        contig_id="13",
        start=8566,
        end=8568,
        strand="+",
        gene_id="ENSG00000139618",
        gene_name="BRCA2",
        transcript_id="ENST00000380152",
        transcript_name="BRCA2-201",
        protein_id="ENSP00000369497",
    )
    expected = RnaPosition(
        _data=ensembl100,
        contig_id="13",
        start=8799,
        end=8801,
        strand="+",
        gene_id="ENSG00000139618",
        gene_name="BRCA2",
        transcript_id="ENST00000380152",
        transcript_name="BRCA2-201",
    )
    assert position.to_rna() == [expected]


def test_cdna_to_rna_positive_strand_across_exon_boundary(ensembl100):
    position = CdnaPosition(
        _data=ensembl100,
        contig_id="4",
        start=67,
        end=69,
        strand="+",
        gene_id="ENSG00000157404",
        gene_name="KIT",
        transcript_id="ENST00000288135",
        transcript_name="KIT-201",
        protein_id="ENSP00000288135",
    )
    expected = RnaPosition(
        _data=ensembl100,
        contig_id="4",
        start=125,
        end=127,
        strand="+",
        gene_id="ENSG00000157404",
        gene_name="KIT",
        transcript_id="ENST00000288135",
        transcript_name="KIT-201",
    )
    assert position.to_rna() == [expected]


def test_cdna_to_rna_positive_strand_cDNA_protein_end(ensembl100):
    position = CdnaPosition(
        _data=ensembl100,
        contig_id="13",
        start=10255,
        end=10257,
        strand="+",
        gene_id="ENSG00000139618",
        gene_name="BRCA2",
        transcript_id="ENST00000380152",
        transcript_name="BRCA2-201",
        protein_id="ENSP00000369497",
    )
    expected = RnaPosition(
        _data=ensembl100,
        contig_id="13",
        start=10488,
        end=10490,
        strand="+",
        gene_id="ENSG00000139618",
        gene_name="BRCA2",
        transcript_id="ENST00000380152",
        transcript_name="BRCA2-201",
    )
    assert position.to_rna() == [expected]


def test_cdna_to_rna_positive_strand_cDNA_protein_start(ensembl100):
    position = CdnaPosition(
        _data=ensembl100,
        contig_id="13",
        start=1,
        end=3,
        strand="+",
        gene_id="ENSG00000139618",
        gene_name="BRCA2",
        transcript_id="ENST00000380152",
        transcript_name="BRCA2-201",
        protein_id="ENSP00000369497",
    )
    expected = RnaPosition(
        _data=ensembl100,
        contig_id="13",
        start=234,
        end=236,
        strand="+",
        gene_id="ENSG00000139618",
        gene_name="BRCA2",
        transcript_id="ENST00000380152",
        transcript_name="BRCA2-201",
    )
    assert position.to_rna() == [expected]


def test_cdna_to_rna_positive_strand_overlapping_genes_different_strands(ensembl100):
    position = CdnaPosition(
        _data=ensembl100,
        contig_id="6",
        start=580,
        end=582,
        strand="+",
        gene_id="ENSG00000137310",
        gene_name="TCF19",
        transcript_id="ENST00000542218",
        transcript_name="TCF19-204",
        protein_id="ENSP00000439397",
    )
    expected = RnaPosition(
        _data=ensembl100,
        contig_id="6",
        start=580,
        end=582,
        strand="+",
        gene_id="ENSG00000137310",
        gene_name="TCF19",
        transcript_id="ENST00000542218",
        transcript_name="TCF19-204",
    )
    assert position.to_rna() == [expected]


def test_dna_to_cdna_negative_strand(ensembl100):
    position = DnaPosition(_data=ensembl100, contig_id="5", start=1294864, end=1294866, strand="-")
    expected = CdnaPosition(
        _data=ensembl100,
        contig_id="5",
        start=124,
        end=126,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000334602",
        transcript_name="TERT-202",
        protein_id="ENSP00000334346",
    )
    assert expected in position.to_cdna()


def test_dna_to_cdna_negative_strand_across_exon_boundary(ensembl100):
    position = DnaPosition(_data=ensembl100, contig_id="5", start=1282623, end=1293313, strand="-")
    expected = CdnaPosition(
        _data=ensembl100,
        contig_id="5",
        start=1573,
        end=1575,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
        protein_id="ENSP00000309572",
    )
    assert expected in position.to_cdna()


def test_dna_to_cdna_negative_strand_cDNA_protein_end(ensembl100):
    position = DnaPosition(_data=ensembl100, contig_id="5", start=1253728, end=1253730, strand="-")
    expected = CdnaPosition(
        _data=ensembl100,
        contig_id="5",
        start=3208,
        end=3210,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000334602",
        transcript_name="TERT-202",
        protein_id="ENSP00000334346",
    )
    assert expected in position.to_cdna()


def test_dna_to_cdna_negative_strand_cDNA_protein_start(ensembl100):
    position = DnaPosition(_data=ensembl100, contig_id="5", start=1294987, end=1294989, strand="-")
    expected = CdnaPosition(
        _data=ensembl100,
        contig_id="5",
        start=1,
        end=3,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000334602",
        transcript_name="TERT-202",
        protein_id="ENSP00000334346",
    )
    assert expected in position.to_cdna()


def test_dna_to_cdna_positive_strand(ensembl100):
    position = DnaPosition(
        _data=ensembl100, contig_id="13", start=32371034, end=32371036, strand="+"
    )
    expected = CdnaPosition(
        _data=ensembl100,
        contig_id="13",
        start=8566,
        end=8568,
        strand="+",
        gene_id="ENSG00000139618",
        gene_name="BRCA2",
        transcript_id="ENST00000380152",
        transcript_name="BRCA2-201",
        protein_id="ENSP00000369497",
    )
    assert expected in position.to_cdna()


def test_dna_to_cdna_positive_strand_across_exon_boundary(ensembl100):
    position = DnaPosition(
        _data=ensembl100, contig_id="4", start=54658081, end=54695513, strand="+"
    )
    expected = CdnaPosition(
        _data=ensembl100,
        contig_id="4",
        start=67,
        end=69,
        strand="+",
        gene_id="ENSG00000157404",
        gene_name="KIT",
        transcript_id="ENST00000288135",
        transcript_name="KIT-201",
        protein_id="ENSP00000288135",
    )
    assert expected in position.to_cdna()


def test_dna_to_cdna_positive_strand_cDNA_protein_end(ensembl100):
    position = DnaPosition(
        _data=ensembl100, contig_id="13", start=32398768, end=32398770, strand="+"
    )
    expected = CdnaPosition(
        _data=ensembl100,
        contig_id="13",
        start=10255,
        end=10257,
        strand="+",
        gene_id="ENSG00000139618",
        gene_name="BRCA2",
        transcript_id="ENST00000380152",
        transcript_name="BRCA2-201",
        protein_id="ENSP00000369497",
    )
    assert expected in position.to_cdna()


def test_dna_to_cdna_positive_strand_cDNA_protein_start(ensembl100):
    position = DnaPosition(
        _data=ensembl100, contig_id="13", start=32316461, end=32316463, strand="+"
    )
    expected = CdnaPosition(
        _data=ensembl100,
        contig_id="13",
        start=1,
        end=3,
        strand="+",
        gene_id="ENSG00000139618",
        gene_name="BRCA2",
        transcript_id="ENST00000380152",
        transcript_name="BRCA2-201",
        protein_id="ENSP00000369497",
    )
    assert expected in position.to_cdna()


def test_dna_to_cdna_positive_strand_overlapping_genes_different_strands(ensembl100):
    position = DnaPosition(
        _data=ensembl100, contig_id="6", start=31166302, end=31166304, strand="+"
    )
    expected = CdnaPosition(
        _data=ensembl100,
        contig_id="6",
        start=580,
        end=582,
        strand="+",
        gene_id="ENSG00000137310",
        gene_name="TCF19",
        transcript_id="ENST00000542218",
        transcript_name="TCF19-204",
        protein_id="ENSP00000439397",
    )
    assert position.to_cdna() == [expected]


def test_dna_to_dna_negative_strand(ensembl100):
    position = DnaPosition(_data=ensembl100, contig_id="5", start=1294864, end=1294866, strand="-")
    assert position.to_dna() == [position]


def test_dna_to_dna_negative_strand_across_exon_boundary(ensembl100):
    position = DnaPosition(_data=ensembl100, contig_id="5", start=1282623, end=1293313, strand="-")
    assert position.to_dna() == [position]


def test_dna_to_dna_negative_strand_cDNA_protein_end(ensembl100):
    position = DnaPosition(_data=ensembl100, contig_id="5", start=1253728, end=1253730, strand="-")
    assert position.to_dna() == [position]


def test_dna_to_dna_negative_strand_cDNA_protein_start(ensembl100):
    position = DnaPosition(_data=ensembl100, contig_id="5", start=1294987, end=1294989, strand="-")
    assert position.to_dna() == [position]


def test_dna_to_dna_negative_strand_overlapping_genes_different_strands(ensembl100):
    position = DnaPosition(
        _data=ensembl100, contig_id="6", start=31166302, end=31166304, strand="-"
    )
    assert position.to_dna() == [position]


def test_dna_to_dna_negative_strand_transcript_end(ensembl100):
    position = DnaPosition(_data=ensembl100, contig_id="5", start=1253167, end=1253169, strand="-")
    assert position.to_dna() == [position]


def test_dna_to_dna_negative_strand_transcript_start(ensembl100):
    position = DnaPosition(_data=ensembl100, contig_id="5", start=1295066, end=1295068, strand="-")
    assert position.to_dna() == [position]


def test_dna_to_dna_positive_strand(ensembl100):
    position = DnaPosition(
        _data=ensembl100, contig_id="13", start=32371034, end=32371036, strand="+"
    )
    assert position.to_dna() == [position]


def test_dna_to_dna_positive_strand_across_exon_boundary(ensembl100):
    position = DnaPosition(
        _data=ensembl100, contig_id="4", start=54658081, end=54695513, strand="+"
    )
    assert position.to_dna() == [position]


def test_dna_to_dna_positive_strand_cDNA_protein_end(ensembl100):
    position = DnaPosition(
        _data=ensembl100, contig_id="13", start=32398768, end=32398770, strand="+"
    )
    assert position.to_dna() == [position]


def test_dna_to_dna_positive_strand_cDNA_protein_start(ensembl100):
    position = DnaPosition(
        _data=ensembl100, contig_id="13", start=32316461, end=32316463, strand="+"
    )
    assert position.to_dna() == [position]


def test_dna_to_dna_positive_strand_overlapping_genes_different_strands(ensembl100):
    position = DnaPosition(
        _data=ensembl100, contig_id="6", start=31166302, end=31166304, strand="+"
    )
    assert position.to_dna() == [position]


def test_dna_to_dna_positive_strand_transcript_end(ensembl100):
    position = DnaPosition(
        _data=ensembl100, contig_id="13", start=32400264, end=32400266, strand="+"
    )
    assert position.to_dna() == [position]


def test_dna_to_dna_positive_strand_transcript_start(ensembl100):
    position = DnaPosition(
        _data=ensembl100, contig_id="13", start=32315474, end=32315476, strand="+"
    )
    assert position.to_dna() == [position]


def test_dna_to_exon_negative_strand(ensembl100):
    position = DnaPosition(_data=ensembl100, contig_id="5", start=1294864, end=1294866, strand="-")
    expected = ExonPosition(
        _data=ensembl100,
        contig_id="5",
        start=1,
        end=1,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
        exon_id="ENSE00003896691",
    )
    assert expected in position.to_exon()


def test_dna_to_exon_negative_strand_cDNA_protein_end(ensembl100):
    position = DnaPosition(_data=ensembl100, contig_id="5", start=1253728, end=1253730, strand="-")
    expected = ExonPosition(
        _data=ensembl100,
        contig_id="5",
        start=16,
        end=16,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
        exon_id="ENSE00001863787",
    )
    assert expected in position.to_exon()


def test_dna_to_exon_negative_strand_cDNA_protein_start(ensembl100):
    position = DnaPosition(_data=ensembl100, contig_id="5", start=1294987, end=1294989, strand="-")
    expected = ExonPosition(
        _data=ensembl100,
        contig_id="5",
        start=1,
        end=1,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
        exon_id="ENSE00003896691",
    )
    assert expected in position.to_exon()


def test_dna_to_exon_negative_strand_overlapping_genes_different_strands(ensembl100):
    position = DnaPosition(
        _data=ensembl100, contig_id="6", start=31166302, end=31166304, strand="-"
    )
    expected = [
        ExonPosition(
            _data=ensembl100,
            contig_id="6",
            start=1,
            end=1,
            strand="-",
            gene_id="ENSG00000204531",
            gene_name="POU5F1",
            transcript_id="ENST00000471529",
            transcript_name="POU5F1-204",
            exon_id="ENSE00002568331",
        ),
        ExonPosition(
            _data=ensembl100,
            contig_id="6",
            start=1,
            end=1,
            strand="-",
            gene_id="ENSG00000204531",
            gene_name="POU5F1",
            transcript_id="ENST00000513407",
            transcript_name="POU5F1-206",
            exon_id="ENSE00002033137",
        ),
    ]
    assert position.to_exon() == expected


def test_dna_to_exon_negative_strand_transcript_end(ensembl100):
    position = DnaPosition(_data=ensembl100, contig_id="5", start=1253167, end=1253169, strand="-")
    expected = ExonPosition(
        _data=ensembl100,
        contig_id="5",
        start=16,
        end=16,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
        exon_id="ENSE00001863787",
    )
    assert expected in position.to_exon()


def test_dna_to_exon_negative_strand_transcript_start(ensembl100):
    position = DnaPosition(_data=ensembl100, contig_id="5", start=1295066, end=1295068, strand="-")
    expected = ExonPosition(
        _data=ensembl100,
        contig_id="5",
        start=1,
        end=1,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
        exon_id="ENSE00003896691",
    )
    assert expected in position.to_exon()


def test_dna_to_exon_positive_strand(ensembl100):
    position = DnaPosition(
        _data=ensembl100, contig_id="13", start=32371034, end=32371036, strand="+"
    )
    expected = ExonPosition(
        _data=ensembl100,
        contig_id="13",
        start=20,
        end=20,
        strand="+",
        gene_id="ENSG00000139618",
        gene_name="BRCA2",
        transcript_id="ENST00000380152",
        transcript_name="BRCA2-201",
        exon_id="ENSE00000939180",
    )
    assert expected in position.to_exon()


def test_dna_to_exon_positive_strand_cDNA_protein_end(ensembl100):
    position = DnaPosition(
        _data=ensembl100, contig_id="13", start=32398768, end=32398770, strand="+"
    )
    expected = ExonPosition(
        _data=ensembl100,
        contig_id="13",
        start=27,
        end=27,
        strand="+",
        gene_id="ENSG00000139618",
        gene_name="BRCA2",
        transcript_id="ENST00000380152",
        transcript_name="BRCA2-201",
        exon_id="ENSE00003717596",
    )
    assert expected in position.to_exon()


def test_dna_to_exon_positive_strand_cDNA_protein_start(ensembl100):
    position = DnaPosition(
        _data=ensembl100, contig_id="13", start=32316461, end=32316463, strand="+"
    )
    expected = ExonPosition(
        _data=ensembl100,
        contig_id="13",
        start=2,
        end=2,
        strand="+",
        gene_id="ENSG00000139618",
        gene_name="BRCA2",
        transcript_id="ENST00000380152",
        transcript_name="BRCA2-201",
        exon_id="ENSE00001484009",
    )
    assert expected in position.to_exon()


def test_dna_to_exon_positive_strand_overlapping_genes_different_strands(ensembl100):
    position = DnaPosition(
        _data=ensembl100, contig_id="6", start=31166302, end=31166304, strand="+"
    )
    expected = ExonPosition(
        _data=ensembl100,
        contig_id="6",
        start=2,
        end=2,
        strand="+",
        gene_id="ENSG00000137310",
        gene_name="TCF19",
        transcript_id="ENST00000542218",
        transcript_name="TCF19-204",
        exon_id="ENSE00002283659",
    )
    assert position.to_exon() == [expected]


def test_dna_to_exon_positive_strand_transcript_end(ensembl100):
    position = DnaPosition(
        _data=ensembl100, contig_id="13", start=32400264, end=32400266, strand="+"
    )
    expected = ExonPosition(
        _data=ensembl100,
        contig_id="13",
        start=27,
        end=27,
        strand="+",
        gene_id="ENSG00000139618",
        gene_name="BRCA2",
        transcript_id="ENST00000380152",
        transcript_name="BRCA2-201",
        exon_id="ENSE00003717596",
    )
    assert expected in position.to_exon()


def test_dna_to_exon_positive_strand_transcript_start(ensembl100):
    position = DnaPosition(
        _data=ensembl100, contig_id="13", start=32315474, end=32315476, strand="+"
    )
    expected = ExonPosition(
        _data=ensembl100,
        contig_id="13",
        start=1,
        end=1,
        strand="+",
        gene_id="ENSG00000139618",
        gene_name="BRCA2",
        transcript_id="ENST00000380152",
        transcript_name="BRCA2-201",
        exon_id="ENSE00001184784",
    )
    assert expected in position.to_exon()


def test_dna_to_protein_negative_strand(ensembl100):
    position = DnaPosition(_data=ensembl100, contig_id="5", start=1294864, end=1294866, strand="-")
    expected = ProteinPosition(
        _data=ensembl100,
        contig_id="5",
        start=42,
        end=42,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
        protein_id="ENSP00000309572",
    )
    assert expected in position.to_protein()


def test_dna_to_protein_negative_strand_across_exon_boundary(ensembl100):
    position = DnaPosition(_data=ensembl100, contig_id="5", start=1282623, end=1293313, strand="-")
    expected = ProteinPosition(
        _data=ensembl100,
        contig_id="5",
        start=525,
        end=525,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
        protein_id="ENSP00000309572",
    )
    assert expected in position.to_protein()


def test_dna_to_protein_negative_strand_cDNA_protein_end(ensembl100):
    position = DnaPosition(_data=ensembl100, contig_id="5", start=1253731, end=1253733, strand="-")
    expected = ProteinPosition(
        _data=ensembl100,
        contig_id="5",
        start=1132,
        end=1132,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
        protein_id="ENSP00000309572",
    )
    assert expected in position.to_protein()


def test_dna_to_protein_negative_strand_stop_codong(ensembl100):
    position = DnaPosition(_data=ensembl100, contig_id="5", start=1253728, end=1253730, strand="-")
    assert position.to_protein() == []


def test_dna_to_protein_negative_strand_cDNA_protein_start(ensembl100):
    position = DnaPosition(_data=ensembl100, contig_id="5", start=1294987, end=1294989, strand="-")
    expected = ProteinPosition(
        _data=ensembl100,
        contig_id="5",
        start=1,
        end=1,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
        protein_id="ENSP00000309572",
    )
    assert expected in position.to_protein()


def test_dna_to_protein_positive_strand(ensembl100):
    position = DnaPosition(
        _data=ensembl100, contig_id="13", start=32371034, end=32371036, strand="+"
    )
    expected = ProteinPosition(
        _data=ensembl100,
        contig_id="13",
        start=2856,
        end=2856,
        strand="+",
        gene_id="ENSG00000139618",
        gene_name="BRCA2",
        transcript_id="ENST00000380152",
        transcript_name="BRCA2-201",
        protein_id="ENSP00000369497",
    )
    assert expected in position.to_protein()


def test_dna_to_protein_positive_strand_across_exon_boundary(ensembl100):
    position = DnaPosition(
        _data=ensembl100, contig_id="4", start=54658081, end=54695513, strand="+"
    )
    expected = ProteinPosition(
        _data=ensembl100,
        contig_id="4",
        start=23,
        end=23,
        strand="+",
        gene_id="ENSG00000157404",
        gene_name="KIT",
        transcript_id="ENST00000288135",
        transcript_name="KIT-201",
        protein_id="ENSP00000288135",
    )
    assert expected in position.to_protein()


def test_dna_to_protein_positive_strand_cDNA_protein_end(ensembl100):
    position = DnaPosition(
        _data=ensembl100, contig_id="13", start=32398765, end=32398767, strand="+"
    )
    expected = ProteinPosition(
        _data=ensembl100,
        contig_id="13",
        start=3418,
        end=3418,
        strand="+",
        gene_id="ENSG00000139618",
        gene_name="BRCA2",
        transcript_id="ENST00000380152",
        transcript_name="BRCA2-201",
        protein_id="ENSP00000369497",
    )
    assert expected in position.to_protein()


def test_dna_to_protein_positive_strand_stop_codon(ensembl100):
    position = DnaPosition(
        _data=ensembl100, contig_id="13", start=32398768, end=32398770, strand="+"
    )
    assert position.to_protein() == []


def test_dna_to_protein_positive_strand_cDNA_protein_start(ensembl100):
    position = DnaPosition(
        _data=ensembl100, contig_id="13", start=32316461, end=32316463, strand="+"
    )
    expected = ProteinPosition(
        _data=ensembl100,
        contig_id="13",
        start=1,
        end=1,
        strand="+",
        gene_id="ENSG00000139618",
        gene_name="BRCA2",
        transcript_id="ENST00000380152",
        transcript_name="BRCA2-201",
        protein_id="ENSP00000369497",
    )
    assert expected in position.to_protein()


def test_dna_to_protein_positive_strand_overlapping_genes_different_strands(ensembl100):
    position = DnaPosition(
        _data=ensembl100, contig_id="6", start=31166302, end=31166304, strand="+"
    )
    expected = ProteinPosition(
        _data=ensembl100,
        contig_id="6",
        start=194,
        end=194,
        strand="+",
        gene_id="ENSG00000137310",
        gene_name="TCF19",
        transcript_id="ENST00000542218",
        transcript_name="TCF19-204",
        protein_id="ENSP00000439397",
    )
    assert position.to_protein() == [expected]


def test_dna_to_rna_negative_strand(ensembl100):
    position = DnaPosition(_data=ensembl100, contig_id="5", start=1294864, end=1294866, strand="-")
    expected = RnaPosition(
        _data=ensembl100,
        contig_id="5",
        start=203,
        end=205,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
    )
    assert expected in position.to_rna()


def test_dna_to_rna_negative_strand_across_exon_boundary(ensembl100):
    position = DnaPosition(_data=ensembl100, contig_id="5", start=1282623, end=1293313, strand="-")
    expected = RnaPosition(
        _data=ensembl100,
        contig_id="5",
        start=1652,
        end=1654,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
    )
    assert expected in position.to_rna()


def test_dna_to_rna_negative_strand_cDNA_protein_end(ensembl100):
    position = DnaPosition(_data=ensembl100, contig_id="5", start=1253728, end=1253730, strand="-")
    expected = RnaPosition(
        _data=ensembl100,
        contig_id="5",
        start=3476,
        end=3478,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
    )
    assert expected in position.to_rna()


def test_dna_to_rna_negative_strand_cDNA_protein_start(ensembl100):
    position = DnaPosition(_data=ensembl100, contig_id="5", start=1294987, end=1294989, strand="-")
    expected = RnaPosition(
        _data=ensembl100,
        contig_id="5",
        start=80,
        end=82,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
    )
    assert expected in position.to_rna()


def test_dna_to_rna_negative_strand_overlapping_genes_different_strands(ensembl100):
    position = DnaPosition(
        _data=ensembl100, contig_id="6", start=31166302, end=31166304, strand="-"
    )
    expected = [
        RnaPosition(
            _data=ensembl100,
            contig_id="6",
            start=535,
            end=537,
            strand="-",
            gene_id="ENSG00000204531",
            gene_name="POU5F1",
            transcript_id="ENST00000471529",
            transcript_name="POU5F1-204",
        ),
        RnaPosition(
            _data=ensembl100,
            contig_id="6",
            start=867,
            end=869,
            strand="-",
            gene_id="ENSG00000204531",
            gene_name="POU5F1",
            transcript_id="ENST00000513407",
            transcript_name="POU5F1-206",
        ),
    ]
    assert position.to_rna() == expected


def test_dna_to_rna_negative_strand_transcript_end(ensembl100):
    position = DnaPosition(_data=ensembl100, contig_id="5", start=1253167, end=1253169, strand="-")
    expected = RnaPosition(
        _data=ensembl100,
        contig_id="5",
        start=4037,
        end=4039,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
    )
    assert expected in position.to_rna()


def test_dna_to_rna_negative_strand_transcript_start(ensembl100):
    position = DnaPosition(_data=ensembl100, contig_id="5", start=1295066, end=1295068, strand="-")
    expected = RnaPosition(
        _data=ensembl100,
        contig_id="5",
        start=1,
        end=3,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
    )
    assert expected in position.to_rna()


def test_dna_to_rna_positive_strand(ensembl100):
    position = DnaPosition(
        _data=ensembl100, contig_id="13", start=32371034, end=32371036, strand="+"
    )
    expected = RnaPosition(
        _data=ensembl100,
        contig_id="13",
        start=8799,
        end=8801,
        strand="+",
        gene_id="ENSG00000139618",
        gene_name="BRCA2",
        transcript_id="ENST00000380152",
        transcript_name="BRCA2-201",
    )
    assert expected in position.to_rna()


def test_dna_to_rna_positive_strand_across_exon_boundary(ensembl100):
    position = DnaPosition(
        _data=ensembl100, contig_id="4", start=54658081, end=54695513, strand="+"
    )
    expected = RnaPosition(
        _data=ensembl100,
        contig_id="4",
        start=125,
        end=127,
        strand="+",
        gene_id="ENSG00000157404",
        gene_name="KIT",
        transcript_id="ENST00000288135",
        transcript_name="KIT-201",
    )
    assert expected in position.to_rna()


def test_dna_to_rna_positive_strand_cDNA_protein_end(ensembl100):
    position = DnaPosition(
        _data=ensembl100, contig_id="13", start=32398768, end=32398770, strand="+"
    )
    expected = RnaPosition(
        _data=ensembl100,
        contig_id="13",
        start=10488,
        end=10490,
        strand="+",
        gene_id="ENSG00000139618",
        gene_name="BRCA2",
        transcript_id="ENST00000380152",
        transcript_name="BRCA2-201",
    )
    assert expected in position.to_rna()


def test_dna_to_rna_positive_strand_cDNA_protein_start(ensembl100):
    position = DnaPosition(
        _data=ensembl100, contig_id="13", start=32316461, end=32316463, strand="+"
    )
    expected = RnaPosition(
        _data=ensembl100,
        contig_id="13",
        start=234,
        end=236,
        strand="+",
        gene_id="ENSG00000139618",
        gene_name="BRCA2",
        transcript_id="ENST00000380152",
        transcript_name="BRCA2-201",
    )
    assert expected in position.to_rna()


def test_dna_to_rna_positive_strand_overlapping_genes_different_strands(ensembl100):
    position = DnaPosition(
        _data=ensembl100, contig_id="6", start=31166302, end=31166304, strand="+"
    )
    expected = RnaPosition(
        _data=ensembl100,
        contig_id="6",
        start=580,
        end=582,
        strand="+",
        gene_id="ENSG00000137310",
        gene_name="TCF19",
        transcript_id="ENST00000542218",
        transcript_name="TCF19-204",
    )
    assert position.to_rna() == [expected]


def test_dna_to_rna_positive_strand_transcript_end(ensembl100):
    position = DnaPosition(
        _data=ensembl100, contig_id="13", start=32400264, end=32400266, strand="+"
    )
    expected = RnaPosition(
        _data=ensembl100,
        contig_id="13",
        start=11984,
        end=11986,
        strand="+",
        gene_id="ENSG00000139618",
        gene_name="BRCA2",
        transcript_id="ENST00000380152",
        transcript_name="BRCA2-201",
    )
    assert expected in position.to_rna()


def test_dna_to_rna_positive_strand_transcript_start(ensembl100):
    position = DnaPosition(
        _data=ensembl100, contig_id="13", start=32315474, end=32315476, strand="+"
    )
    expected = RnaPosition(
        _data=ensembl100,
        contig_id="13",
        start=1,
        end=3,
        strand="+",
        gene_id="ENSG00000139618",
        gene_name="BRCA2",
        transcript_id="ENST00000380152",
        transcript_name="BRCA2-201",
    )
    assert expected in position.to_rna()


def test_exon_to_cdna_negative_strand(ensembl100):
    position = ExonPosition(
        _data=ensembl100,
        contig_id="5",
        start=2,
        end=2,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
        exon_id="ENSE00001197112",
    )
    expected = CdnaPosition(
        _data=ensembl100,
        contig_id="5",
        start=220,
        end=1573,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
        protein_id="ENSP00000309572",
    )
    assert position.to_cdna() == [expected]


def test_exon_to_cdna_positive_strand(ensembl100):
    position = ExonPosition(
        _data=ensembl100,
        contig_id="13",
        start=3,
        end=3,
        strand="+",
        gene_id="ENSG00000139618",
        gene_name="BRCA2",
        transcript_id="ENST00000380152",
        transcript_name="BRCA2-201",
        exon_id="ENSE00003666217",
    )
    expected = CdnaPosition(
        _data=ensembl100,
        contig_id="13",
        start=68,
        end=316,
        strand="+",
        gene_id="ENSG00000139618",
        gene_name="BRCA2",
        transcript_id="ENST00000380152",
        transcript_name="BRCA2-201",
        protein_id="ENSP00000369497",
    )
    assert position.to_cdna() == [expected]


def test_exon_to_dna_negative_strand(ensembl100):
    position = ExonPosition(
        _data=ensembl100,
        contig_id="5",
        start=2,
        end=2,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
        exon_id="ENSE00001197112",
    )
    expected = DnaPosition(_data=ensembl100, contig_id="5", start=1293313, end=1294666, strand="-")
    assert position.to_dna() == [expected]


def test_exon_to_dna_positive_strand(ensembl100):
    position = ExonPosition(
        _data=ensembl100,
        contig_id="13",
        start=3,
        end=3,
        strand="+",
        gene_id="ENSG00000139618",
        gene_name="BRCA2",
        transcript_id="ENST00000380152",
        transcript_name="BRCA2-201",
        exon_id="ENSE00003666217",
    )
    expected = DnaPosition(
        _data=ensembl100, contig_id="13", start=32319077, end=32319325, strand="+"
    )
    assert position.to_dna() == [expected]


def test_exon_to_exon_negative_strand(ensembl100):
    position = ExonPosition(
        _data=ensembl100,
        contig_id="5",
        start=2,
        end=2,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
        exon_id="ENSE00001197112",
    )
    assert position.to_exon() == [position]


def test_exon_to_exon_negative_strand_cDNA_protein_end(ensembl100):
    position = ExonPosition(
        _data=ensembl100,
        contig_id="5",
        start=16,
        end=16,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
        exon_id="ENSE00001863787",
    )
    assert position.to_exon() == [position]


def test_exon_to_exon_negative_strand_cDNA_protein_start(ensembl100):
    position = ExonPosition(
        _data=ensembl100,
        contig_id="5",
        start=1,
        end=1,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
        exon_id="ENSE00003896691",
    )
    assert position.to_exon() == [position]


def test_exon_to_exon_negative_strand_overlapping_genes_different_strands(ensembl100):
    position = ExonPosition(
        _data=ensembl100,
        contig_id="6",
        start=1,
        end=1,
        strand="-",
        gene_id="ENSG00000204531",
        gene_name="POU5F1",
        transcript_id="ENST00000471529",
        transcript_name="POU5F1-204",
        exon_id="ENSE00002568331",
    )
    assert position.to_exon() == [position]


def test_exon_to_exon_negative_strand_transcript_end(ensembl100):
    position = ExonPosition(
        _data=ensembl100,
        contig_id="5",
        start=16,
        end=16,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
        exon_id="ENSE00001863787",
    )
    assert position.to_exon() == [position]


def test_exon_to_exon_negative_strand_transcript_start(ensembl100):
    position = ExonPosition(
        _data=ensembl100,
        contig_id="5",
        start=1,
        end=1,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
        exon_id="ENSE00003896691",
    )
    assert position.to_exon() == [position]


def test_exon_to_exon_positive_strand(ensembl100):
    position = ExonPosition(
        _data=ensembl100,
        contig_id="13",
        start=20,
        end=20,
        strand="+",
        gene_id="ENSG00000139618",
        gene_name="BRCA2",
        transcript_id="ENST00000380152",
        transcript_name="BRCA2-201",
        exon_id="ENSE00000939180",
    )
    assert position.to_exon() == [position]


def test_exon_to_exon_positive_strand_cDNA_protein_end(ensembl100):
    position = ExonPosition(
        _data=ensembl100,
        contig_id="13",
        start=27,
        end=27,
        strand="+",
        gene_id="ENSG00000139618",
        gene_name="BRCA2",
        transcript_id="ENST00000380152",
        transcript_name="BRCA2-201",
        exon_id="ENSE00003717596",
    )
    assert position.to_exon() == [position]


def test_exon_to_exon_positive_strand_cDNA_protein_start(ensembl100):
    position = ExonPosition(
        _data=ensembl100,
        contig_id="13",
        start=2,
        end=2,
        strand="+",
        gene_id="ENSG00000139618",
        gene_name="BRCA2",
        transcript_id="ENST00000380152",
        transcript_name="BRCA2-201",
        exon_id="ENSE00001484009",
    )
    assert position.to_exon() == [position]


def test_exon_to_exon_positive_strand_overlapping_genes_different_strands(ensembl100):
    position = ExonPosition(
        _data=ensembl100,
        contig_id="6",
        start=2,
        end=2,
        strand="+",
        gene_id="ENSG00000137310",
        gene_name="TCF19",
        transcript_id="ENST00000542218",
        transcript_name="TCF19-204",
        exon_id="ENSE00002283659",
    )
    assert position.to_exon() == [position]


def test_exon_to_exon_positive_strand_transcript_end(ensembl100):
    position = ExonPosition(
        _data=ensembl100,
        contig_id="13",
        start=27,
        end=27,
        strand="+",
        gene_id="ENSG00000139618",
        gene_name="BRCA2",
        transcript_id="ENST00000380152",
        transcript_name="BRCA2-201",
        exon_id="ENSE00003717596",
    )
    assert position.to_exon() == [position]


def test_exon_to_exon_positive_strand_transcript_start(ensembl100):
    position = ExonPosition(
        _data=ensembl100,
        contig_id="13",
        start=1,
        end=1,
        strand="+",
        gene_id="ENSG00000139618",
        gene_name="BRCA2",
        transcript_id="ENST00000380152",
        transcript_name="BRCA2-201",
        exon_id="ENSE00001184784",
    )
    assert position.to_exon() == [position]


def test_exon_to_protein_negative_strand(ensembl100):
    position = ExonPosition(
        _data=ensembl100,
        contig_id="5",
        start=2,
        end=2,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
        exon_id="ENSE00001197112",
    )
    expected = ProteinPosition(
        _data=ensembl100,
        contig_id="5",
        start=74,
        end=525,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
        protein_id="ENSP00000309572",
    )
    assert position.to_protein() == [expected]


def test_exon_to_protein_positive_strand(ensembl100):
    position = ExonPosition(
        _data=ensembl100,
        contig_id="13",
        start=3,
        end=3,
        strand="+",
        gene_id="ENSG00000139618",
        gene_name="BRCA2",
        transcript_id="ENST00000380152",
        transcript_name="BRCA2-201",
        exon_id="ENSE00003666217",
    )
    expected = ProteinPosition(
        _data=ensembl100,
        contig_id="13",
        start=23,
        end=106,
        strand="+",
        gene_id="ENSG00000139618",
        gene_name="BRCA2",
        transcript_id="ENST00000380152",
        transcript_name="BRCA2-201",
        protein_id="ENSP00000369497",
    )
    assert position.to_protein() == [expected]


def test_exon_to_rna_negative_strand(ensembl100):
    position = ExonPosition(
        _data=ensembl100,
        contig_id="5",
        start=2,
        end=2,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
        exon_id="ENSE00001197112",
    )
    expected = RnaPosition(
        _data=ensembl100,
        contig_id="5",
        start=299,
        end=1652,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
    )
    assert position.to_rna() == [expected]


def test_exon_to_rna_positive_strand(ensembl100):
    position = ExonPosition(
        _data=ensembl100,
        contig_id="13",
        start=3,
        end=3,
        strand="+",
        gene_id="ENSG00000139618",
        gene_name="BRCA2",
        transcript_id="ENST00000380152",
        transcript_name="BRCA2-201",
        exon_id="ENSE00003666217",
    )
    expected = RnaPosition(
        _data=ensembl100,
        contig_id="13",
        start=301,
        end=549,
        strand="+",
        gene_id="ENSG00000139618",
        gene_name="BRCA2",
        transcript_id="ENST00000380152",
        transcript_name="BRCA2-201",
    )
    assert position.to_rna() == [expected]


def test_protein_to_cdna_negative_strand(ensembl100):
    position = ProteinPosition(
        _data=ensembl100,
        contig_id="5",
        start=42,
        end=42,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
        protein_id="ENSP00000309572",
    )
    expected = CdnaPosition(
        _data=ensembl100,
        contig_id="5",
        start=124,
        end=126,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
        protein_id="ENSP00000309572",
    )
    assert position.to_cdna() == [expected]


def test_protein_to_cdna_negative_strand_across_exon_boundary(ensembl100):
    position = ProteinPosition(
        _data=ensembl100,
        contig_id="5",
        start=525,
        end=525,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
        protein_id="ENSP00000309572",
    )
    expected = CdnaPosition(
        _data=ensembl100,
        contig_id="5",
        start=1573,
        end=1575,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
        protein_id="ENSP00000309572",
    )
    assert position.to_cdna() == [expected]


def test_protein_to_cdna_negative_strand_cDNA_protein_end(ensembl100):
    position = ProteinPosition(
        _data=ensembl100,
        contig_id="5",
        start=1132,
        end=1132,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
        protein_id="ENSP00000309572",
    )
    expected = CdnaPosition(
        _data=ensembl100,
        contig_id="5",
        start=3394,
        end=3396,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
        protein_id="ENSP00000309572",
    )
    assert position.to_cdna() == [expected]


def test_protein_to_cdna_negative_strand_cDNA_protein_start(ensembl100):
    position = ProteinPosition(
        _data=ensembl100,
        contig_id="5",
        start=1,
        end=1,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
        protein_id="ENSP00000309572",
    )
    expected = CdnaPosition(
        _data=ensembl100,
        contig_id="5",
        start=1,
        end=3,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
        protein_id="ENSP00000309572",
    )
    assert position.to_cdna() == [expected]


def test_protein_to_cdna_positive_strand(ensembl100):
    position = ProteinPosition(
        _data=ensembl100,
        contig_id="13",
        start=2856,
        end=2856,
        strand="+",
        gene_id="ENSG00000139618",
        gene_name="BRCA2",
        transcript_id="ENST00000380152",
        transcript_name="BRCA2-201",
        protein_id="ENSP00000369497",
    )
    expected = CdnaPosition(
        _data=ensembl100,
        contig_id="13",
        start=8566,
        end=8568,
        strand="+",
        gene_id="ENSG00000139618",
        gene_name="BRCA2",
        transcript_id="ENST00000380152",
        transcript_name="BRCA2-201",
        protein_id="ENSP00000369497",
    )
    assert position.to_cdna() == [expected]


def test_protein_to_cdna_positive_strand_across_exon_boundary(ensembl100):
    position = ProteinPosition(
        _data=ensembl100,
        contig_id="4",
        start=23,
        end=23,
        strand="+",
        gene_id="ENSG00000157404",
        gene_name="KIT",
        transcript_id="ENST00000288135",
        transcript_name="KIT-201",
        protein_id="ENSP00000288135",
    )
    expected = CdnaPosition(
        _data=ensembl100,
        contig_id="4",
        start=67,
        end=69,
        strand="+",
        gene_id="ENSG00000157404",
        gene_name="KIT",
        transcript_id="ENST00000288135",
        transcript_name="KIT-201",
        protein_id="ENSP00000288135",
    )
    assert position.to_cdna() == [expected]


def test_protein_to_cdna_positive_strand_cDNA_protein_end(ensembl100):
    position = ProteinPosition(
        _data=ensembl100,
        contig_id="13",
        start=3418,
        end=3418,
        strand="+",
        gene_id="ENSG00000139618",
        gene_name="BRCA2",
        transcript_id="ENST00000380152",
        transcript_name="BRCA2-201",
        protein_id="ENSP00000369497",
    )
    expected = CdnaPosition(
        _data=ensembl100,
        contig_id="13",
        start=10252,
        end=10254,
        strand="+",
        gene_id="ENSG00000139618",
        gene_name="BRCA2",
        transcript_id="ENST00000380152",
        transcript_name="BRCA2-201",
        protein_id="ENSP00000369497",
    )
    assert position.to_cdna() == [expected]


def test_protein_to_cdna_positive_strand_cDNA_protein_start(ensembl100):
    position = ProteinPosition(
        _data=ensembl100,
        contig_id="13",
        start=1,
        end=1,
        strand="+",
        gene_id="ENSG00000139618",
        gene_name="BRCA2",
        transcript_id="ENST00000380152",
        transcript_name="BRCA2-201",
        protein_id="ENSP00000369497",
    )
    expected = CdnaPosition(
        _data=ensembl100,
        contig_id="13",
        start=1,
        end=3,
        strand="+",
        gene_id="ENSG00000139618",
        gene_name="BRCA2",
        transcript_id="ENST00000380152",
        transcript_name="BRCA2-201",
        protein_id="ENSP00000369497",
    )
    assert position.to_cdna() == [expected]


def test_protein_to_cdna_positive_strand_overlapping_genes_different_strands(ensembl100):
    position = ProteinPosition(
        _data=ensembl100,
        contig_id="6",
        start=194,
        end=194,
        strand="+",
        gene_id="ENSG00000137310",
        gene_name="TCF19",
        transcript_id="ENST00000542218",
        transcript_name="TCF19-204",
        protein_id="ENSP00000439397",
    )
    expected = CdnaPosition(
        _data=ensembl100,
        contig_id="6",
        start=580,
        end=582,
        strand="+",
        gene_id="ENSG00000137310",
        gene_name="TCF19",
        transcript_id="ENST00000542218",
        transcript_name="TCF19-204",
        protein_id="ENSP00000439397",
    )
    assert position.to_cdna() == [expected]


def test_protein_to_dna_negative_strand(ensembl100):
    position = ProteinPosition(
        _data=ensembl100,
        contig_id="5",
        start=42,
        end=42,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
        protein_id="ENSP00000309572",
    )
    expected = DnaPosition(_data=ensembl100, contig_id="5", start=1294864, end=1294866, strand="-")
    assert position.to_dna() == [expected]


def test_protein_to_dna_negative_strand_across_exon_boundary(ensembl100):
    position = ProteinPosition(
        _data=ensembl100,
        contig_id="5",
        start=525,
        end=525,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
        protein_id="ENSP00000309572",
    )
    expected = DnaPosition(_data=ensembl100, contig_id="5", start=1282623, end=1293313, strand="-")
    assert position.to_dna() == [expected]


def test_protein_to_dna_negative_strand_cDNA_protein_end(ensembl100):
    position = ProteinPosition(
        _data=ensembl100,
        contig_id="5",
        start=1132,
        end=1132,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
        protein_id="ENSP00000309572",
    )
    expected = DnaPosition(_data=ensembl100, contig_id="5", start=1253731, end=1253733, strand="-")
    assert position.to_dna() == [expected]


def test_protein_to_dna_negative_strand_cDNA_protein_start(ensembl100):
    position = ProteinPosition(
        _data=ensembl100,
        contig_id="5",
        start=1,
        end=1,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
        protein_id="ENSP00000309572",
    )
    expected = DnaPosition(_data=ensembl100, contig_id="5", start=1294987, end=1294989, strand="-")
    assert position.to_dna() == [expected]


def test_protein_to_dna_positive_strand(ensembl100):
    position = ProteinPosition(
        _data=ensembl100,
        contig_id="13",
        start=2856,
        end=2856,
        strand="+",
        gene_id="ENSG00000139618",
        gene_name="BRCA2",
        transcript_id="ENST00000380152",
        transcript_name="BRCA2-201",
        protein_id="ENSP00000369497",
    )
    expected = DnaPosition(
        _data=ensembl100, contig_id="13", start=32371034, end=32371036, strand="+"
    )
    assert position.to_dna() == [expected]


def test_protein_to_dna_positive_strand_across_exon_boundary(ensembl100):
    position = ProteinPosition(
        _data=ensembl100,
        contig_id="4",
        start=23,
        end=23,
        strand="+",
        gene_id="ENSG00000157404",
        gene_name="KIT",
        transcript_id="ENST00000288135",
        transcript_name="KIT-201",
        protein_id="ENSP00000288135",
    )
    expected = DnaPosition(
        _data=ensembl100, contig_id="4", start=54658081, end=54695513, strand="+"
    )
    assert position.to_dna() == [expected]


def test_protein_to_dna_positive_strand_cDNA_protein_end(ensembl100):
    position = ProteinPosition(
        _data=ensembl100,
        contig_id="13",
        start=3418,
        end=3418,
        strand="+",
        gene_id="ENSG00000139618",
        gene_name="BRCA2",
        transcript_id="ENST00000380152",
        transcript_name="BRCA2-201",
        protein_id="ENSP00000369497",
    )
    expected = DnaPosition(
        _data=ensembl100, contig_id="13", start=32398765, end=32398767, strand="+"
    )
    assert position.to_dna() == [expected]


def test_protein_to_dna_positive_strand_cDNA_protein_start(ensembl100):
    position = ProteinPosition(
        _data=ensembl100,
        contig_id="13",
        start=1,
        end=1,
        strand="+",
        gene_id="ENSG00000139618",
        gene_name="BRCA2",
        transcript_id="ENST00000380152",
        transcript_name="BRCA2-201",
        protein_id="ENSP00000369497",
    )
    expected = DnaPosition(
        _data=ensembl100, contig_id="13", start=32316461, end=32316463, strand="+"
    )
    assert position.to_dna() == [expected]


def test_protein_to_dna_positive_strand_overlapping_genes_different_strands(ensembl100):
    position = ProteinPosition(
        _data=ensembl100,
        contig_id="6",
        start=194,
        end=194,
        strand="+",
        gene_id="ENSG00000137310",
        gene_name="TCF19",
        transcript_id="ENST00000542218",
        transcript_name="TCF19-204",
        protein_id="ENSP00000439397",
    )
    expected = DnaPosition(
        _data=ensembl100, contig_id="6", start=31166302, end=31166304, strand="+"
    )
    assert position.to_dna() == [expected]


def test_protein_to_exon_negative_strand(ensembl100):
    position = ProteinPosition(
        _data=ensembl100,
        contig_id="5",
        start=42,
        end=42,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
        protein_id="ENSP00000309572",
    )
    expected = ExonPosition(
        _data=ensembl100,
        contig_id="5",
        start=1,
        end=1,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
        exon_id="ENSE00003896691",
    )
    assert position.to_exon() == [expected]


def test_protein_to_exon_negative_strand_cDNA_protein_end(ensembl100):
    position = ProteinPosition(
        _data=ensembl100,
        contig_id="5",
        start=1132,
        end=1132,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
        protein_id="ENSP00000309572",
    )
    expected = ExonPosition(
        _data=ensembl100,
        contig_id="5",
        start=16,
        end=16,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
        exon_id="ENSE00001863787",
    )
    assert position.to_exon() == [expected]


def test_protein_to_exon_negative_strand_cDNA_protein_start(ensembl100):
    position = ProteinPosition(
        _data=ensembl100,
        contig_id="5",
        start=1,
        end=1,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
        protein_id="ENSP00000309572",
    )
    expected = ExonPosition(
        _data=ensembl100,
        contig_id="5",
        start=1,
        end=1,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
        exon_id="ENSE00003896691",
    )
    assert position.to_exon() == [expected]


def test_protein_to_exon_positive_strand(ensembl100):
    position = ProteinPosition(
        _data=ensembl100,
        contig_id="13",
        start=2856,
        end=2856,
        strand="+",
        gene_id="ENSG00000139618",
        gene_name="BRCA2",
        transcript_id="ENST00000380152",
        transcript_name="BRCA2-201",
        protein_id="ENSP00000369497",
    )
    expected = ExonPosition(
        _data=ensembl100,
        contig_id="13",
        start=20,
        end=20,
        strand="+",
        gene_id="ENSG00000139618",
        gene_name="BRCA2",
        transcript_id="ENST00000380152",
        transcript_name="BRCA2-201",
        exon_id="ENSE00000939180",
    )
    assert position.to_exon() == [expected]


def test_protein_to_exon_positive_strand_cDNA_protein_end(ensembl100):
    position = ProteinPosition(
        _data=ensembl100,
        contig_id="13",
        start=3418,
        end=3418,
        strand="+",
        gene_id="ENSG00000139618",
        gene_name="BRCA2",
        transcript_id="ENST00000380152",
        transcript_name="BRCA2-201",
        protein_id="ENSP00000369497",
    )
    expected = ExonPosition(
        _data=ensembl100,
        contig_id="13",
        start=27,
        end=27,
        strand="+",
        gene_id="ENSG00000139618",
        gene_name="BRCA2",
        transcript_id="ENST00000380152",
        transcript_name="BRCA2-201",
        exon_id="ENSE00003717596",
    )
    assert position.to_exon() == [expected]


def test_protein_to_exon_positive_strand_cDNA_protein_start(ensembl100):
    position = ProteinPosition(
        _data=ensembl100,
        contig_id="13",
        start=1,
        end=1,
        strand="+",
        gene_id="ENSG00000139618",
        gene_name="BRCA2",
        transcript_id="ENST00000380152",
        transcript_name="BRCA2-201",
        protein_id="ENSP00000369497",
    )
    expected = ExonPosition(
        _data=ensembl100,
        contig_id="13",
        start=2,
        end=2,
        strand="+",
        gene_id="ENSG00000139618",
        gene_name="BRCA2",
        transcript_id="ENST00000380152",
        transcript_name="BRCA2-201",
        exon_id="ENSE00001484009",
    )
    assert position.to_exon() == [expected]


def test_protein_to_exon_positive_strand_overlapping_genes_different_strands(ensembl100):
    position = ProteinPosition(
        _data=ensembl100,
        contig_id="6",
        start=194,
        end=194,
        strand="+",
        gene_id="ENSG00000137310",
        gene_name="TCF19",
        transcript_id="ENST00000542218",
        transcript_name="TCF19-204",
        protein_id="ENSP00000439397",
    )
    expected = ExonPosition(
        _data=ensembl100,
        contig_id="6",
        start=2,
        end=2,
        strand="+",
        gene_id="ENSG00000137310",
        gene_name="TCF19",
        transcript_id="ENST00000542218",
        transcript_name="TCF19-204",
        exon_id="ENSE00002283659",
    )
    assert position.to_exon() == [expected]


def test_protein_to_protein_negative_strand(ensembl100):
    position = ProteinPosition(
        _data=ensembl100,
        contig_id="5",
        start=42,
        end=42,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
        protein_id="ENSP00000309572",
    )
    assert position.to_protein() == [position]


def test_protein_to_protein_negative_strand_across_exon_boundary(ensembl100):
    position = ProteinPosition(
        _data=ensembl100,
        contig_id="5",
        start=525,
        end=525,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
        protein_id="ENSP00000309572",
    )
    assert position.to_protein() == [position]


def test_protein_to_protein_negative_strand_cDNA_protein_end(ensembl100):
    position = ProteinPosition(
        _data=ensembl100,
        contig_id="5",
        start=1132,
        end=1132,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
        protein_id="ENSP00000309572",
    )
    assert position.to_protein() == [position]


def test_protein_to_protein_negative_strand_stop_codon(ensembl100):
    position = ProteinPosition(
        _data=ensembl100,
        contig_id="5",
        start=1133,
        end=1133,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
        protein_id="ENSP00000309572",
    )
    assert position.to_protein() == []


def test_protein_to_protein_negative_strand_cDNA_protein_start(ensembl100):
    position = ProteinPosition(
        _data=ensembl100,
        contig_id="5",
        start=1,
        end=1,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
        protein_id="ENSP00000309572",
    )
    assert position.to_protein() == [position]


def test_protein_to_protein_positive_strand(ensembl100):
    position = ProteinPosition(
        _data=ensembl100,
        contig_id="13",
        start=2856,
        end=2856,
        strand="+",
        gene_id="ENSG00000139618",
        gene_name="BRCA2",
        transcript_id="ENST00000380152",
        transcript_name="BRCA2-201",
        protein_id="ENSP00000369497",
    )
    assert position.to_protein() == [position]


def test_protein_to_protein_positive_strand_across_exon_boundary(ensembl100):
    position = ProteinPosition(
        _data=ensembl100,
        contig_id="4",
        start=23,
        end=23,
        strand="+",
        gene_id="ENSG00000157404",
        gene_name="KIT",
        transcript_id="ENST00000288135",
        transcript_name="KIT-201",
        protein_id="ENSP00000288135",
    )
    assert position.to_protein() == [position]


def test_protein_to_protein_positive_strand_cDNA_protein_end(ensembl100):
    position = ProteinPosition(
        _data=ensembl100,
        contig_id="13",
        start=3418,
        end=3418,
        strand="+",
        gene_id="ENSG00000139618",
        gene_name="BRCA2",
        transcript_id="ENST00000380152",
        transcript_name="BRCA2-201",
        protein_id="ENSP00000369497",
    )
    assert position.to_protein() == [position]


def test_protein_to_protein_positive_strand_stop_codon(ensembl100):
    position = ProteinPosition(
        _data=ensembl100,
        contig_id="13",
        start=3419,
        end=3419,
        strand="+",
        gene_id="ENSG00000139618",
        gene_name="BRCA2",
        transcript_id="ENST00000380152",
        transcript_name="BRCA2-201",
        protein_id="ENSP00000369497",
    )
    assert position.to_protein() == []


def test_protein_to_protein_positive_strand_cDNA_protein_start(ensembl100):
    position = ProteinPosition(
        _data=ensembl100,
        contig_id="13",
        start=1,
        end=1,
        strand="+",
        gene_id="ENSG00000139618",
        gene_name="BRCA2",
        transcript_id="ENST00000380152",
        transcript_name="BRCA2-201",
        protein_id="ENSP00000369497",
    )
    assert position.to_protein() == [position]


def test_protein_to_protein_positive_strand_overlapping_genes_different_strands(ensembl100):
    position = ProteinPosition(
        _data=ensembl100,
        contig_id="6",
        start=194,
        end=194,
        strand="+",
        gene_id="ENSG00000137310",
        gene_name="TCF19",
        transcript_id="ENST00000542218",
        transcript_name="TCF19-204",
        protein_id="ENSP00000439397",
    )
    assert position.to_protein() == [position]


def test_protein_to_rna_negative_strand(ensembl100):
    position = ProteinPosition(
        _data=ensembl100,
        contig_id="5",
        start=42,
        end=42,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
        protein_id="ENSP00000309572",
    )
    expected = RnaPosition(
        _data=ensembl100,
        contig_id="5",
        start=203,
        end=205,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
    )
    assert position.to_rna() == [expected]


def test_protein_to_rna_negative_strand_across_exon_boundary(ensembl100):
    position = ProteinPosition(
        _data=ensembl100,
        contig_id="5",
        start=525,
        end=525,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
        protein_id="ENSP00000309572",
    )
    expected = RnaPosition(
        _data=ensembl100,
        contig_id="5",
        start=1652,
        end=1654,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
    )
    assert position.to_rna() == [expected]


def test_protein_to_rna_negative_strand_cDNA_protein_end(ensembl100):
    position = ProteinPosition(
        _data=ensembl100,
        contig_id="5",
        start=1132,
        end=1132,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
        protein_id="ENSP00000309572",
    )
    expected = RnaPosition(
        _data=ensembl100,
        contig_id="5",
        start=3473,
        end=3475,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
    )
    assert position.to_rna() == [expected]


def test_protein_to_rna_negative_strand_cDNA_protein_start(ensembl100):
    position = ProteinPosition(
        _data=ensembl100,
        contig_id="5",
        start=1,
        end=1,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
        protein_id="ENSP00000309572",
    )
    expected = RnaPosition(
        _data=ensembl100,
        contig_id="5",
        start=80,
        end=82,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
    )
    assert position.to_rna() == [expected]


def test_protein_to_rna_positive_strand(ensembl100):
    position = ProteinPosition(
        _data=ensembl100,
        contig_id="13",
        start=2856,
        end=2856,
        strand="+",
        gene_id="ENSG00000139618",
        gene_name="BRCA2",
        transcript_id="ENST00000380152",
        transcript_name="BRCA2-201",
        protein_id="ENSP00000369497",
    )
    expected = RnaPosition(
        _data=ensembl100,
        contig_id="13",
        start=8799,
        end=8801,
        strand="+",
        gene_id="ENSG00000139618",
        gene_name="BRCA2",
        transcript_id="ENST00000380152",
        transcript_name="BRCA2-201",
    )
    assert position.to_rna() == [expected]


def test_protein_to_rna_positive_strand_across_exon_boundary(ensembl100):
    position = ProteinPosition(
        _data=ensembl100,
        contig_id="4",
        start=23,
        end=23,
        strand="+",
        gene_id="ENSG00000157404",
        gene_name="KIT",
        transcript_id="ENST00000288135",
        transcript_name="KIT-201",
        protein_id="ENSP00000288135",
    )
    expected = RnaPosition(
        _data=ensembl100,
        contig_id="4",
        start=125,
        end=127,
        strand="+",
        gene_id="ENSG00000157404",
        gene_name="KIT",
        transcript_id="ENST00000288135",
        transcript_name="KIT-201",
    )
    assert position.to_rna() == [expected]


def test_protein_to_rna_positive_strand_cDNA_protein_end(ensembl100):
    position = ProteinPosition(
        _data=ensembl100,
        contig_id="13",
        start=3418,
        end=3418,
        strand="+",
        gene_id="ENSG00000139618",
        gene_name="BRCA2",
        transcript_id="ENST00000380152",
        transcript_name="BRCA2-201",
        protein_id="ENSP00000369497",
    )
    expected = RnaPosition(
        _data=ensembl100,
        contig_id="13",
        start=10485,
        end=10487,
        strand="+",
        gene_id="ENSG00000139618",
        gene_name="BRCA2",
        transcript_id="ENST00000380152",
        transcript_name="BRCA2-201",
    )
    assert position.to_rna() == [expected]


def test_protein_to_rna_positive_strand_cDNA_protein_start(ensembl100):
    position = ProteinPosition(
        _data=ensembl100,
        contig_id="13",
        start=1,
        end=1,
        strand="+",
        gene_id="ENSG00000139618",
        gene_name="BRCA2",
        transcript_id="ENST00000380152",
        transcript_name="BRCA2-201",
        protein_id="ENSP00000369497",
    )
    expected = RnaPosition(
        _data=ensembl100,
        contig_id="13",
        start=234,
        end=236,
        strand="+",
        gene_id="ENSG00000139618",
        gene_name="BRCA2",
        transcript_id="ENST00000380152",
        transcript_name="BRCA2-201",
    )
    assert position.to_rna() == [expected]


def test_protein_to_rna_positive_strand_overlapping_genes_different_strands(ensembl100):
    position = ProteinPosition(
        _data=ensembl100,
        contig_id="6",
        start=194,
        end=194,
        strand="+",
        gene_id="ENSG00000137310",
        gene_name="TCF19",
        transcript_id="ENST00000542218",
        transcript_name="TCF19-204",
        protein_id="ENSP00000439397",
    )
    expected = RnaPosition(
        _data=ensembl100,
        contig_id="6",
        start=580,
        end=582,
        strand="+",
        gene_id="ENSG00000137310",
        gene_name="TCF19",
        transcript_id="ENST00000542218",
        transcript_name="TCF19-204",
    )
    assert position.to_rna() == [expected]


def test_rna_to_cdna_negative_strand(ensembl100):
    position = RnaPosition(
        _data=ensembl100,
        contig_id="5",
        start=203,
        end=205,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
    )
    expected = CdnaPosition(
        _data=ensembl100,
        contig_id="5",
        start=124,
        end=126,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
        protein_id="ENSP00000309572",
    )
    assert position.to_cdna() == [expected]


def test_rna_to_cdna_negative_strand_across_exon_boundary(ensembl100):
    position = RnaPosition(
        _data=ensembl100,
        contig_id="5",
        start=1652,
        end=1654,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
    )
    expected = CdnaPosition(
        _data=ensembl100,
        contig_id="5",
        start=1573,
        end=1575,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
        protein_id="ENSP00000309572",
    )
    assert position.to_cdna() == [expected]


def test_rna_to_cdna_negative_strand_cDNA_protein_end(ensembl100):
    position = RnaPosition(
        _data=ensembl100,
        contig_id="5",
        start=3476,
        end=3478,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
    )
    expected = CdnaPosition(
        _data=ensembl100,
        contig_id="5",
        start=3397,
        end=3399,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
        protein_id="ENSP00000309572",
    )
    assert position.to_cdna() == [expected]


def test_rna_to_cdna_negative_strand_cDNA_protein_start(ensembl100):
    position = RnaPosition(
        _data=ensembl100,
        contig_id="5",
        start=80,
        end=82,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
    )
    expected = CdnaPosition(
        _data=ensembl100,
        contig_id="5",
        start=1,
        end=3,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
        protein_id="ENSP00000309572",
    )
    assert position.to_cdna() == [expected]


def test_rna_to_cdna_positive_strand(ensembl100):
    position = RnaPosition(
        _data=ensembl100,
        contig_id="13",
        start=8799,
        end=8801,
        strand="+",
        gene_id="ENSG00000139618",
        gene_name="BRCA2",
        transcript_id="ENST00000380152",
        transcript_name="BRCA2-201",
    )
    expected = CdnaPosition(
        _data=ensembl100,
        contig_id="13",
        start=8566,
        end=8568,
        strand="+",
        gene_id="ENSG00000139618",
        gene_name="BRCA2",
        transcript_id="ENST00000380152",
        transcript_name="BRCA2-201",
        protein_id="ENSP00000369497",
    )
    assert position.to_cdna() == [expected]


def test_rna_to_cdna_positive_strand_across_exon_boundary(ensembl100):
    position = RnaPosition(
        _data=ensembl100,
        contig_id="4",
        start=125,
        end=127,
        strand="+",
        gene_id="ENSG00000157404",
        gene_name="KIT",
        transcript_id="ENST00000288135",
        transcript_name="KIT-201",
    )
    expected = CdnaPosition(
        _data=ensembl100,
        contig_id="4",
        start=67,
        end=69,
        strand="+",
        gene_id="ENSG00000157404",
        gene_name="KIT",
        transcript_id="ENST00000288135",
        transcript_name="KIT-201",
        protein_id="ENSP00000288135",
    )
    assert position.to_cdna() == [expected]


def test_rna_to_cdna_positive_strand_cDNA_protein_end(ensembl100):
    position = RnaPosition(
        _data=ensembl100,
        contig_id="13",
        start=10488,
        end=10490,
        strand="+",
        gene_id="ENSG00000139618",
        gene_name="BRCA2",
        transcript_id="ENST00000380152",
        transcript_name="BRCA2-201",
    )
    expected = CdnaPosition(
        _data=ensembl100,
        contig_id="13",
        start=10255,
        end=10257,
        strand="+",
        gene_id="ENSG00000139618",
        gene_name="BRCA2",
        transcript_id="ENST00000380152",
        transcript_name="BRCA2-201",
        protein_id="ENSP00000369497",
    )
    assert position.to_cdna() == [expected]


def test_rna_to_cdna_positive_strand_cDNA_protein_start(ensembl100):
    position = RnaPosition(
        _data=ensembl100,
        contig_id="13",
        start=234,
        end=236,
        strand="+",
        gene_id="ENSG00000139618",
        gene_name="BRCA2",
        transcript_id="ENST00000380152",
        transcript_name="BRCA2-201",
    )
    expected = CdnaPosition(
        _data=ensembl100,
        contig_id="13",
        start=1,
        end=3,
        strand="+",
        gene_id="ENSG00000139618",
        gene_name="BRCA2",
        transcript_id="ENST00000380152",
        transcript_name="BRCA2-201",
        protein_id="ENSP00000369497",
    )
    assert position.to_cdna() == [expected]


def test_rna_to_cdna_positive_strand_overlapping_genes_different_strands(ensembl100):
    position = RnaPosition(
        _data=ensembl100,
        contig_id="6",
        start=580,
        end=582,
        strand="+",
        gene_id="ENSG00000137310",
        gene_name="TCF19",
        transcript_id="ENST00000542218",
        transcript_name="TCF19-204",
    )
    expected = CdnaPosition(
        _data=ensembl100,
        contig_id="6",
        start=580,
        end=582,
        strand="+",
        gene_id="ENSG00000137310",
        gene_name="TCF19",
        transcript_id="ENST00000542218",
        transcript_name="TCF19-204",
        protein_id="ENSP00000439397",
    )
    assert position.to_cdna() == [expected]


def test_rna_to_dna_negative_strand(ensembl100):
    position = RnaPosition(
        _data=ensembl100,
        contig_id="5",
        start=203,
        end=205,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
    )
    expected = DnaPosition(_data=ensembl100, contig_id="5", start=1294864, end=1294866, strand="-")
    assert position.to_dna() == [expected]


def test_rna_to_dna_negative_strand_across_exon_boundary(ensembl100):
    position = RnaPosition(
        _data=ensembl100,
        contig_id="5",
        start=1652,
        end=1654,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
    )
    expected = DnaPosition(_data=ensembl100, contig_id="5", start=1282623, end=1293313, strand="-")
    assert position.to_dna() == [expected]


def test_rna_to_dna_negative_strand_cDNA_protein_end(ensembl100):
    position = RnaPosition(
        _data=ensembl100,
        contig_id="5",
        start=3476,
        end=3478,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
    )
    expected = DnaPosition(_data=ensembl100, contig_id="5", start=1253728, end=1253730, strand="-")
    assert position.to_dna() == [expected]


def test_rna_to_dna_negative_strand_cDNA_protein_start(ensembl100):
    position = RnaPosition(
        _data=ensembl100,
        contig_id="5",
        start=80,
        end=82,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
    )
    expected = DnaPosition(_data=ensembl100, contig_id="5", start=1294987, end=1294989, strand="-")
    assert position.to_dna() == [expected]


def test_rna_to_dna_negative_strand_overlapping_genes_different_strands(ensembl100):
    position = RnaPosition(
        _data=ensembl100,
        contig_id="6",
        start=535,
        end=537,
        strand="-",
        gene_id="ENSG00000204531",
        gene_name="POU5F1",
        transcript_id="ENST00000471529",
        transcript_name="POU5F1-204",
    )
    expected = DnaPosition(
        _data=ensembl100, contig_id="6", start=31166302, end=31166304, strand="-"
    )
    assert position.to_dna() == [expected]


def test_rna_to_dna_negative_strand_transcript_end(ensembl100):
    position = RnaPosition(
        _data=ensembl100,
        contig_id="5",
        start=4037,
        end=4039,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
    )
    expected = DnaPosition(_data=ensembl100, contig_id="5", start=1253167, end=1253169, strand="-")
    assert position.to_dna() == [expected]


def test_rna_to_dna_negative_strand_transcript_start(ensembl100):
    position = RnaPosition(
        _data=ensembl100,
        contig_id="5",
        start=1,
        end=3,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
    )
    expected = DnaPosition(_data=ensembl100, contig_id="5", start=1295066, end=1295068, strand="-")
    assert position.to_dna() == [expected]


def test_rna_to_dna_positive_strand(ensembl100):
    position = RnaPosition(
        _data=ensembl100,
        contig_id="13",
        start=8799,
        end=8801,
        strand="+",
        gene_id="ENSG00000139618",
        gene_name="BRCA2",
        transcript_id="ENST00000380152",
        transcript_name="BRCA2-201",
    )
    expected = DnaPosition(
        _data=ensembl100, contig_id="13", start=32371034, end=32371036, strand="+"
    )
    assert position.to_dna() == [expected]


def test_rna_to_dna_positive_strand_across_exon_boundary(ensembl100):
    position = RnaPosition(
        _data=ensembl100,
        contig_id="4",
        start=125,
        end=127,
        strand="+",
        gene_id="ENSG00000157404",
        gene_name="KIT",
        transcript_id="ENST00000288135",
        transcript_name="KIT-201",
    )
    expected = DnaPosition(
        _data=ensembl100, contig_id="4", start=54658081, end=54695513, strand="+"
    )
    assert position.to_dna() == [expected]


def test_rna_to_dna_positive_strand_cDNA_protein_end(ensembl100):
    position = RnaPosition(
        _data=ensembl100,
        contig_id="13",
        start=10488,
        end=10490,
        strand="+",
        gene_id="ENSG00000139618",
        gene_name="BRCA2",
        transcript_id="ENST00000380152",
        transcript_name="BRCA2-201",
    )
    expected = DnaPosition(
        _data=ensembl100, contig_id="13", start=32398768, end=32398770, strand="+"
    )
    assert position.to_dna() == [expected]


def test_rna_to_dna_positive_strand_cDNA_protein_start(ensembl100):
    position = RnaPosition(
        _data=ensembl100,
        contig_id="13",
        start=234,
        end=236,
        strand="+",
        gene_id="ENSG00000139618",
        gene_name="BRCA2",
        transcript_id="ENST00000380152",
        transcript_name="BRCA2-201",
    )
    expected = DnaPosition(
        _data=ensembl100, contig_id="13", start=32316461, end=32316463, strand="+"
    )
    assert position.to_dna() == [expected]


def test_rna_to_dna_positive_strand_overlapping_genes_different_strands(ensembl100):
    position = RnaPosition(
        _data=ensembl100,
        contig_id="6",
        start=580,
        end=582,
        strand="+",
        gene_id="ENSG00000137310",
        gene_name="TCF19",
        transcript_id="ENST00000542218",
        transcript_name="TCF19-204",
    )
    expected = DnaPosition(
        _data=ensembl100, contig_id="6", start=31166302, end=31166304, strand="+"
    )
    assert position.to_dna() == [expected]


def test_rna_to_dna_positive_strand_transcript_end(ensembl100):
    position = RnaPosition(
        _data=ensembl100,
        contig_id="13",
        start=11984,
        end=11986,
        strand="+",
        gene_id="ENSG00000139618",
        gene_name="BRCA2",
        transcript_id="ENST00000380152",
        transcript_name="BRCA2-201",
    )
    expected = DnaPosition(
        _data=ensembl100, contig_id="13", start=32400264, end=32400266, strand="+"
    )
    assert position.to_dna() == [expected]


def test_rna_to_dna_positive_strand_transcript_start(ensembl100):
    position = RnaPosition(
        _data=ensembl100,
        contig_id="13",
        start=1,
        end=3,
        strand="+",
        gene_id="ENSG00000139618",
        gene_name="BRCA2",
        transcript_id="ENST00000380152",
        transcript_name="BRCA2-201",
    )
    expected = DnaPosition(
        _data=ensembl100, contig_id="13", start=32315474, end=32315476, strand="+"
    )
    assert position.to_dna() == [expected]


def test_rna_to_exon_negative_strand(ensembl100):
    position = RnaPosition(
        _data=ensembl100,
        contig_id="5",
        start=203,
        end=205,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
    )
    expected = ExonPosition(
        _data=ensembl100,
        contig_id="5",
        start=1,
        end=1,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
        exon_id="ENSE00003896691",
    )
    assert position.to_exon() == [expected]


def test_rna_to_exon_negative_strand_cDNA_protein_end(ensembl100):
    position = RnaPosition(
        _data=ensembl100,
        contig_id="5",
        start=3476,
        end=3478,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
    )
    expected = ExonPosition(
        _data=ensembl100,
        contig_id="5",
        start=16,
        end=16,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
        exon_id="ENSE00001863787",
    )
    assert position.to_exon() == [expected]


def test_rna_to_exon_negative_strand_cDNA_protein_start(ensembl100):
    position = RnaPosition(
        _data=ensembl100,
        contig_id="5",
        start=80,
        end=82,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
    )
    expected = ExonPosition(
        _data=ensembl100,
        contig_id="5",
        start=1,
        end=1,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
        exon_id="ENSE00003896691",
    )
    assert position.to_exon() == [expected]


def test_rna_to_exon_negative_strand_overlapping_genes_different_strands(ensembl100):
    position = RnaPosition(
        _data=ensembl100,
        contig_id="6",
        start=535,
        end=537,
        strand="-",
        gene_id="ENSG00000204531",
        gene_name="POU5F1",
        transcript_id="ENST00000471529",
        transcript_name="POU5F1-204",
    )
    expected = ExonPosition(
        _data=ensembl100,
        contig_id="6",
        start=1,
        end=1,
        strand="-",
        gene_id="ENSG00000204531",
        gene_name="POU5F1",
        transcript_id="ENST00000471529",
        transcript_name="POU5F1-204",
        exon_id="ENSE00002568331",
    )
    assert position.to_exon() == [expected]


def test_rna_to_exon_negative_strand_transcript_end(ensembl100):
    position = RnaPosition(
        _data=ensembl100,
        contig_id="5",
        start=4037,
        end=4039,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
    )
    expected = ExonPosition(
        _data=ensembl100,
        contig_id="5",
        start=16,
        end=16,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
        exon_id="ENSE00001863787",
    )
    assert position.to_exon() == [expected]


def test_rna_to_exon_negative_strand_transcript_start(ensembl100):
    position = RnaPosition(
        _data=ensembl100,
        contig_id="5",
        start=1,
        end=3,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
    )
    expected = ExonPosition(
        _data=ensembl100,
        contig_id="5",
        start=1,
        end=1,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
        exon_id="ENSE00003896691",
    )
    assert position.to_exon() == [expected]


def test_rna_to_exon_positive_strand(ensembl100):
    position = RnaPosition(
        _data=ensembl100,
        contig_id="13",
        start=8799,
        end=8801,
        strand="+",
        gene_id="ENSG00000139618",
        gene_name="BRCA2",
        transcript_id="ENST00000380152",
        transcript_name="BRCA2-201",
    )
    expected = ExonPosition(
        _data=ensembl100,
        contig_id="13",
        start=20,
        end=20,
        strand="+",
        gene_id="ENSG00000139618",
        gene_name="BRCA2",
        transcript_id="ENST00000380152",
        transcript_name="BRCA2-201",
        exon_id="ENSE00000939180",
    )
    assert position.to_exon() == [expected]


def test_rna_to_exon_positive_strand_cDNA_protein_end(ensembl100):
    position = RnaPosition(
        _data=ensembl100,
        contig_id="13",
        start=10488,
        end=10490,
        strand="+",
        gene_id="ENSG00000139618",
        gene_name="BRCA2",
        transcript_id="ENST00000380152",
        transcript_name="BRCA2-201",
    )
    expected = ExonPosition(
        _data=ensembl100,
        contig_id="13",
        start=27,
        end=27,
        strand="+",
        gene_id="ENSG00000139618",
        gene_name="BRCA2",
        transcript_id="ENST00000380152",
        transcript_name="BRCA2-201",
        exon_id="ENSE00003717596",
    )
    assert position.to_exon() == [expected]


def test_rna_to_exon_positive_strand_cDNA_protein_start(ensembl100):
    position = RnaPosition(
        _data=ensembl100,
        contig_id="13",
        start=234,
        end=236,
        strand="+",
        gene_id="ENSG00000139618",
        gene_name="BRCA2",
        transcript_id="ENST00000380152",
        transcript_name="BRCA2-201",
    )
    expected = ExonPosition(
        _data=ensembl100,
        contig_id="13",
        start=2,
        end=2,
        strand="+",
        gene_id="ENSG00000139618",
        gene_name="BRCA2",
        transcript_id="ENST00000380152",
        transcript_name="BRCA2-201",
        exon_id="ENSE00001484009",
    )
    assert position.to_exon() == [expected]


def test_rna_to_exon_positive_strand_overlapping_genes_different_strands(ensembl100):
    position = RnaPosition(
        _data=ensembl100,
        contig_id="6",
        start=580,
        end=582,
        strand="+",
        gene_id="ENSG00000137310",
        gene_name="TCF19",
        transcript_id="ENST00000542218",
        transcript_name="TCF19-204",
    )
    expected = ExonPosition(
        _data=ensembl100,
        contig_id="6",
        start=2,
        end=2,
        strand="+",
        gene_id="ENSG00000137310",
        gene_name="TCF19",
        transcript_id="ENST00000542218",
        transcript_name="TCF19-204",
        exon_id="ENSE00002283659",
    )
    assert position.to_exon() == [expected]


def test_rna_to_exon_positive_strand_transcript_end(ensembl100):
    position = RnaPosition(
        _data=ensembl100,
        contig_id="13",
        start=11984,
        end=11986,
        strand="+",
        gene_id="ENSG00000139618",
        gene_name="BRCA2",
        transcript_id="ENST00000380152",
        transcript_name="BRCA2-201",
    )
    expected = ExonPosition(
        _data=ensembl100,
        contig_id="13",
        start=27,
        end=27,
        strand="+",
        gene_id="ENSG00000139618",
        gene_name="BRCA2",
        transcript_id="ENST00000380152",
        transcript_name="BRCA2-201",
        exon_id="ENSE00003717596",
    )
    assert position.to_exon() == [expected]


def test_rna_to_exon_positive_strand_transcript_start(ensembl100):
    position = RnaPosition(
        _data=ensembl100,
        contig_id="13",
        start=1,
        end=3,
        strand="+",
        gene_id="ENSG00000139618",
        gene_name="BRCA2",
        transcript_id="ENST00000380152",
        transcript_name="BRCA2-201",
    )
    expected = ExonPosition(
        _data=ensembl100,
        contig_id="13",
        start=1,
        end=1,
        strand="+",
        gene_id="ENSG00000139618",
        gene_name="BRCA2",
        transcript_id="ENST00000380152",
        transcript_name="BRCA2-201",
        exon_id="ENSE00001184784",
    )
    assert position.to_exon() == [expected]


def test_rna_to_protein_negative_strand(ensembl100):
    position = RnaPosition(
        _data=ensembl100,
        contig_id="5",
        start=203,
        end=205,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
    )
    expected = ProteinPosition(
        _data=ensembl100,
        contig_id="5",
        start=42,
        end=42,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
        protein_id="ENSP00000309572",
    )
    assert position.to_protein() == [expected]


def test_rna_to_protein_negative_strand_across_exon_boundary(ensembl100):
    position = RnaPosition(
        _data=ensembl100,
        contig_id="5",
        start=1652,
        end=1654,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
    )
    expected = ProteinPosition(
        _data=ensembl100,
        contig_id="5",
        start=525,
        end=525,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
        protein_id="ENSP00000309572",
    )
    assert position.to_protein() == [expected]


def test_rna_to_protein_negative_strand_cDNA_protein_end(ensembl100):
    position = RnaPosition(
        _data=ensembl100,
        contig_id="5",
        start=3473,
        end=3475,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
    )
    expected = ProteinPosition(
        _data=ensembl100,
        contig_id="5",
        start=1132,
        end=1132,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
        protein_id="ENSP00000309572",
    )
    assert position.to_protein() == [expected]


def test_rna_to_protein_negative_strand_stop_codon(ensembl100):
    position = RnaPosition(
        _data=ensembl100,
        contig_id="5",
        start=3476,
        end=3478,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
    )
    assert position.to_protein() == []


def test_rna_to_protein_negative_strand_cDNA_protein_start(ensembl100):
    position = RnaPosition(
        _data=ensembl100,
        contig_id="5",
        start=80,
        end=82,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
    )
    expected = ProteinPosition(
        _data=ensembl100,
        contig_id="5",
        start=1,
        end=1,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
        protein_id="ENSP00000309572",
    )
    assert position.to_protein() == [expected]


def test_rna_to_protein_positive_strand(ensembl100):
    position = RnaPosition(
        _data=ensembl100,
        contig_id="13",
        start=8799,
        end=8801,
        strand="+",
        gene_id="ENSG00000139618",
        gene_name="BRCA2",
        transcript_id="ENST00000380152",
        transcript_name="BRCA2-201",
    )
    expected = ProteinPosition(
        _data=ensembl100,
        contig_id="13",
        start=2856,
        end=2856,
        strand="+",
        gene_id="ENSG00000139618",
        gene_name="BRCA2",
        transcript_id="ENST00000380152",
        transcript_name="BRCA2-201",
        protein_id="ENSP00000369497",
    )
    assert position.to_protein() == [expected]


def test_rna_to_protein_positive_strand_across_exon_boundary(ensembl100):
    position = RnaPosition(
        _data=ensembl100,
        contig_id="4",
        start=125,
        end=127,
        strand="+",
        gene_id="ENSG00000157404",
        gene_name="KIT",
        transcript_id="ENST00000288135",
        transcript_name="KIT-201",
    )
    expected = ProteinPosition(
        _data=ensembl100,
        contig_id="4",
        start=23,
        end=23,
        strand="+",
        gene_id="ENSG00000157404",
        gene_name="KIT",
        transcript_id="ENST00000288135",
        transcript_name="KIT-201",
        protein_id="ENSP00000288135",
    )
    assert position.to_protein() == [expected]


def test_rna_to_protein_positive_strand_cDNA_protein_end(ensembl100):
    position = RnaPosition(
        _data=ensembl100,
        contig_id="13",
        start=10485,
        end=10487,
        strand="+",
        gene_id="ENSG00000139618",
        gene_name="BRCA2",
        transcript_id="ENST00000380152",
        transcript_name="BRCA2-201",
    )
    expected = ProteinPosition(
        _data=ensembl100,
        contig_id="13",
        start=3418,
        end=3418,
        strand="+",
        gene_id="ENSG00000139618",
        gene_name="BRCA2",
        transcript_id="ENST00000380152",
        transcript_name="BRCA2-201",
        protein_id="ENSP00000369497",
    )
    assert position.to_protein() == [expected]


def test_rna_to_protein_positive_strand_stop_codon(ensembl100):
    position = RnaPosition(
        _data=ensembl100,
        contig_id="13",
        start=10488,
        end=10490,
        strand="+",
        gene_id="ENSG00000139618",
        gene_name="BRCA2",
        transcript_id="ENST00000380152",
        transcript_name="BRCA2-201",
    )
    assert position.to_protein() == []


def test_rna_to_protein_positive_strand_cDNA_protein_start(ensembl100):
    position = RnaPosition(
        _data=ensembl100,
        contig_id="13",
        start=234,
        end=236,
        strand="+",
        gene_id="ENSG00000139618",
        gene_name="BRCA2",
        transcript_id="ENST00000380152",
        transcript_name="BRCA2-201",
    )
    expected = ProteinPosition(
        _data=ensembl100,
        contig_id="13",
        start=1,
        end=1,
        strand="+",
        gene_id="ENSG00000139618",
        gene_name="BRCA2",
        transcript_id="ENST00000380152",
        transcript_name="BRCA2-201",
        protein_id="ENSP00000369497",
    )
    assert position.to_protein() == [expected]


def test_rna_to_protein_positive_strand_overlapping_genes_different_strands(ensembl100):
    position = RnaPosition(
        _data=ensembl100,
        contig_id="6",
        start=580,
        end=582,
        strand="+",
        gene_id="ENSG00000137310",
        gene_name="TCF19",
        transcript_id="ENST00000542218",
        transcript_name="TCF19-204",
    )
    expected = ProteinPosition(
        _data=ensembl100,
        contig_id="6",
        start=194,
        end=194,
        strand="+",
        gene_id="ENSG00000137310",
        gene_name="TCF19",
        transcript_id="ENST00000542218",
        transcript_name="TCF19-204",
        protein_id="ENSP00000439397",
    )
    assert position.to_protein() == [expected]


def test_rna_to_rna_negative_strand(ensembl100):
    position = RnaPosition(
        _data=ensembl100,
        contig_id="5",
        start=203,
        end=205,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
    )
    assert position.to_rna() == [position]


def test_rna_to_rna_negative_strand_across_exon_boundary(ensembl100):
    position = RnaPosition(
        _data=ensembl100,
        contig_id="5",
        start=1652,
        end=1654,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
    )
    assert position.to_rna() == [position]


def test_rna_to_rna_negative_strand_cDNA_protein_end(ensembl100):
    position = RnaPosition(
        _data=ensembl100,
        contig_id="5",
        start=3476,
        end=3478,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
    )
    assert position.to_rna() == [position]


def test_rna_to_rna_negative_strand_cDNA_protein_start(ensembl100):
    position = RnaPosition(
        _data=ensembl100,
        contig_id="5",
        start=80,
        end=82,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
    )
    assert position.to_rna() == [position]


def test_rna_to_rna_negative_strand_overlapping_genes_different_strands(ensembl100):
    position = RnaPosition(
        _data=ensembl100,
        contig_id="6",
        start=535,
        end=537,
        strand="-",
        gene_id="ENSG00000204531",
        gene_name="POU5F1",
        transcript_id="ENST00000471529",
        transcript_name="POU5F1-204",
    )
    assert position.to_rna() == [position]


def test_rna_to_rna_negative_strand_transcript_end(ensembl100):
    position = RnaPosition(
        _data=ensembl100,
        contig_id="5",
        start=4037,
        end=4039,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
    )
    assert position.to_rna() == [position]


def test_rna_to_rna_negative_strand_transcript_start(ensembl100):
    position = RnaPosition(
        _data=ensembl100,
        contig_id="5",
        start=1,
        end=3,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
    )
    assert position.to_rna() == [position]


def test_rna_to_rna_positive_strand(ensembl100):
    position = RnaPosition(
        _data=ensembl100,
        contig_id="13",
        start=8799,
        end=8801,
        strand="+",
        gene_id="ENSG00000139618",
        gene_name="BRCA2",
        transcript_id="ENST00000380152",
        transcript_name="BRCA2-201",
    )
    assert position.to_rna() == [position]


def test_rna_to_rna_positive_strand_across_exon_boundary(ensembl100):
    position = RnaPosition(
        _data=ensembl100,
        contig_id="4",
        start=125,
        end=127,
        strand="+",
        gene_id="ENSG00000157404",
        gene_name="KIT",
        transcript_id="ENST00000288135",
        transcript_name="KIT-201",
    )
    assert position.to_rna() == [position]


def test_rna_to_rna_positive_strand_cDNA_protein_end(ensembl100):
    position = RnaPosition(
        _data=ensembl100,
        contig_id="13",
        start=10488,
        end=10490,
        strand="+",
        gene_id="ENSG00000139618",
        gene_name="BRCA2",
        transcript_id="ENST00000380152",
        transcript_name="BRCA2-201",
    )
    assert position.to_rna() == [position]


def test_rna_to_rna_positive_strand_cDNA_protein_start(ensembl100):
    position = RnaPosition(
        _data=ensembl100,
        contig_id="13",
        start=234,
        end=236,
        strand="+",
        gene_id="ENSG00000139618",
        gene_name="BRCA2",
        transcript_id="ENST00000380152",
        transcript_name="BRCA2-201",
    )
    assert position.to_rna() == [position]


def test_rna_to_rna_positive_strand_overlapping_genes_different_strands(ensembl100):
    position = RnaPosition(
        _data=ensembl100,
        contig_id="6",
        start=580,
        end=582,
        strand="+",
        gene_id="ENSG00000137310",
        gene_name="TCF19",
        transcript_id="ENST00000542218",
        transcript_name="TCF19-204",
    )
    assert position.to_rna() == [position]


def test_rna_to_rna_positive_strand_transcript_end(ensembl100):
    position = RnaPosition(
        _data=ensembl100,
        contig_id="13",
        start=11984,
        end=11986,
        strand="+",
        gene_id="ENSG00000139618",
        gene_name="BRCA2",
        transcript_id="ENST00000380152",
        transcript_name="BRCA2-201",
    )
    assert position.to_rna() == [position]


def test_rna_to_rna_positive_strand_transcript_start(ensembl100):
    position = RnaPosition(
        _data=ensembl100,
        contig_id="13",
        start=1,
        end=3,
        strand="+",
        gene_id="ENSG00000139618",
        gene_name="BRCA2",
        transcript_id="ENST00000380152",
        transcript_name="BRCA2-201",
    )
    assert position.to_rna() == [position]
