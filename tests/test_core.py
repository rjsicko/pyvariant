import pandas as pd
import pytest
from constants import (
    CANONICAL_TRANSCRIPT,
    CONTIG_ALIAS,
    EXON_ALIAS,
    GENE_ALIAS,
    PROTEIN_ALIAS,
    TEST_ENS100_CDNA_FASTA,
    TEST_ENS100_DNA_FASTA,
    TEST_ENS100_GTF,
    TEST_ENS100_NCRNA_FASTA,
    TEST_ENS100_PEP_FASTA,
    TRANSCRIPT_ALIAS,
)
from pyfaidx import Fasta

from variant_map.constants import (
    CONTIG_ID,
    EXON_ID,
    GENE_ID,
    GENE_NAME,
    PROTEIN_ID,
    TRANSCRIPT_ID,
    TRANSCRIPT_NAME,
)
from variant_map.core import Core
from variant_map.positions import (
    CdnaPosition,
    CdnaSubstitution,
    DnaPosition,
    DnaSubstitution,
    ExonFusion,
    ExonPosition,
    ProteinPosition,
    RnaPosition,
)


# -------------------------------------------------------------------------------------------------
# test init
# -------------------------------------------------------------------------------------------------
def test_init():
    obj = Core(
        gtf=TEST_ENS100_GTF,
        cds=[TEST_ENS100_CDNA_FASTA],
        dna=[TEST_ENS100_DNA_FASTA],
        peptide=[TEST_ENS100_PEP_FASTA],
        rna=[TEST_ENS100_NCRNA_FASTA],
        canonical_transcript=CANONICAL_TRANSCRIPT,
        contig_alias=CONTIG_ALIAS,
        exon_alias=EXON_ALIAS,
        gene_alias=GENE_ALIAS,
        protein_alias=PROTEIN_ALIAS,
        transcript_alias=TRANSCRIPT_ALIAS,
    )
    assert isinstance(obj.df, pd.DataFrame)
    assert isinstance(obj.cds_fasta, list)
    assert isinstance(obj.cds_fasta[0], Fasta)
    assert isinstance(obj.dna_fasta, list)
    assert isinstance(obj.dna_fasta[0], Fasta)
    assert isinstance(obj.protein_fasta, list)
    assert isinstance(obj.protein_fasta[0], Fasta)
    assert isinstance(obj.rna_fasta, list)
    assert isinstance(obj.rna_fasta[0], Fasta)
    assert isinstance(obj._canonical_transcript, list)
    assert isinstance(obj._contig_alias, dict)
    assert isinstance(obj._exon_alias, dict)
    assert isinstance(obj._gene_alias, dict)
    assert isinstance(obj._protein_alias, dict)
    assert isinstance(obj._transcript_alias, dict)


# -------------------------------------------------------------------------------------------------
# test all_<feature_symbol>s
# -------------------------------------------------------------------------------------------------
def test_all_contig_ids(ensembl100):
    result = ensembl100.all_contig_ids()
    assert isinstance(result, list)
    assert "16" in result


def test_all_exon_ids(ensembl100):
    result = ensembl100.all_exon_ids()
    assert isinstance(result, list)
    assert "ENSE00003826864" in result


def test_all_gene_ids(ensembl100):
    result = ensembl100.all_gene_ids()
    assert isinstance(result, list)
    assert "ENSG00000149925" in result


def test_all_gene_names(ensembl100):
    result = ensembl100.all_gene_names()
    assert isinstance(result, list)
    assert "ALDOA" in result


def test_all_protein_ids(ensembl100):
    result = ensembl100.all_protein_ids()
    assert isinstance(result, list)
    assert "ENSP00000494188" in result


def test_all_transcript_ids(ensembl100):
    result = ensembl100.all_transcript_ids()
    assert isinstance(result, list)
    assert "ENST00000643777" in result


def test_all_transcript_names(ensembl100):
    result = ensembl100.all_transcript_names()
    assert isinstance(result, list)
    assert "ALDOA-219" in result


# -------------------------------------------------------------------------------------------------
# test get_<feature_symbol>s
# -------------------------------------------------------------------------------------------------
def test_contig_ids_from_contig_id(ensembl100):
    result = ensembl100.contig_ids("16")
    assert isinstance(result, list)
    assert "16" in result


def test_contig_ids_from_contig_with_chr(ensembl100):
    ret_by_id = ensembl100.contig_ids("chr4")
    assert "4" in ret_by_id


def test_contig_ids_from_exon_id(ensembl100):
    result = ensembl100.contig_ids("ENSE00003826864")
    assert isinstance(result, list)
    assert "16" in result


def test_contig_ids_from_gene(ensembl100):
    ret_by_id = ensembl100.contig_ids("ENSG00000149925")
    ret_by_name = ensembl100.contig_ids("ALDOA")
    assert "16" in ret_by_id
    assert "16" in ret_by_name
    assert ret_by_id == ret_by_name


def test_contig_ids_from_protein_id(ensembl100):
    result = ensembl100.contig_ids("ENSP00000494188")
    assert isinstance(result, list)
    assert "16" in result


def test_contig_ids_from_transcript(ensembl100):
    ret_by_id = ensembl100.contig_ids("ENST00000643777")
    ret_by_name = ensembl100.contig_ids("ALDOA-219")
    assert "16" in ret_by_id
    assert "16" in ret_by_name
    assert ret_by_id == ret_by_name


def test_exon_ids_from_contig_id(ensembl100):
    result = ensembl100.exon_ids("16")
    assert isinstance(result, list)
    assert "ENSE00003826864" in result


def test_exon_ids_from_exon_id(ensembl100):
    result = ensembl100.exon_ids("ENSE00003826864")
    assert isinstance(result, list)
    assert "ENSE00003826864" in result


def test_exon_ids_from_gene(ensembl100):
    ret_by_id = ensembl100.exon_ids("ENSG00000149925")
    ret_by_name = ensembl100.exon_ids("ALDOA")
    assert "ENSE00003826864" in ret_by_id
    assert "ENSE00003826864" in ret_by_name
    assert ret_by_id == ret_by_name


def test_exon_ids_from_protein_id(ensembl100):
    result = ensembl100.exon_ids("ENSP00000494188")
    assert isinstance(result, list)
    assert "ENSE00003826864" in result


def test_exon_ids_from_transcript(ensembl100):
    ret_by_id = ensembl100.exon_ids("ENST00000643777")
    ret_by_name = ensembl100.exon_ids("ALDOA-219")
    assert "ENSE00003826864" in ret_by_id
    assert "ENSE00003826864" in ret_by_name
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
    result = ensembl100.gene_ids("16")
    assert isinstance(result, list)
    assert "ENSG00000149925" in result


def test_gene_ids_from_exon_id(ensembl100):
    result = ensembl100.gene_ids("ENSE00003826864")
    assert isinstance(result, list)
    assert "ENSG00000149925" in result


def test_gene_ids_from_gene(ensembl100):
    ret_by_id = ensembl100.gene_ids("ENSG00000149925")
    ret_by_name = ensembl100.gene_ids("ALDOA")
    assert "ENSG00000149925" in ret_by_id
    assert "ENSG00000149925" in ret_by_name
    assert ret_by_id == ret_by_name


def test_gene_ids_from_protein_id(ensembl100):
    result = ensembl100.gene_ids("ENSP00000494188")
    assert isinstance(result, list)
    assert "ENSG00000149925" in result


def test_gene_ids_from_transcript(ensembl100):
    ret_by_id = ensembl100.gene_ids("ENST00000643777")
    ret_by_name = ensembl100.gene_ids("ALDOA-219")
    assert "ENSG00000149925" in ret_by_id
    assert "ENSG00000149925" in ret_by_name
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
    result = ensembl100.gene_names("16")
    assert isinstance(result, list)
    assert "ALDOA" in result


def test_gene_names_from_exon_id(ensembl100):
    result = ensembl100.gene_names("ENSE00003826864")
    assert isinstance(result, list)
    assert "ALDOA" in result


def test_gene_names_from_gene(ensembl100):
    ret_by_id = ensembl100.gene_names("ENSG00000149925")
    ret_by_name = ensembl100.gene_names("ALDOA")
    assert "ALDOA" in ret_by_id
    assert "ALDOA" in ret_by_name
    assert ret_by_id == ret_by_name


def test_gene_names_from_protein_id(ensembl100):
    result = ensembl100.gene_names("ENSP00000494188")
    assert isinstance(result, list)
    assert "ALDOA" in result


def test_gene_names_from_refseq_transcript(ensembl100):
    ret_by_id = ensembl100.gene_names("NM_001370259")
    assert "MEN1" in ret_by_id
    ret_by_id = ensembl100.gene_names("NM_007194")
    assert "CHEK2" in ret_by_id
    ret_by_id = ensembl100.gene_names("NM_003072")
    assert "SMARCA4" in ret_by_id
    ret_by_id = ensembl100.gene_names("NM_000314")
    assert "PTEN" in ret_by_id


def test_gene_names_from_transcript(ensembl100):
    ret_by_id = ensembl100.gene_names("ENST00000643777")
    ret_by_name = ensembl100.gene_names("ALDOA-219")
    assert "ALDOA" in ret_by_id
    assert "ALDOA" in ret_by_name
    assert ret_by_id == ret_by_name


def test_protein_ids_from_contig_id(ensembl100):
    result = ensembl100.protein_ids("16")
    assert isinstance(result, list)
    assert "ENSP00000494188" in result


def test_protein_ids_from_exon_id(ensembl100):
    result = ensembl100.protein_ids("ENSE00003826864")
    assert isinstance(result, list)
    assert "ENSP00000494188" in result


def test_protein_ids_from_gene(ensembl100):
    ret_by_id = ensembl100.protein_ids("ENSG00000149925")
    ret_by_name = ensembl100.protein_ids("ALDOA")
    assert "ENSP00000494188" in ret_by_id
    assert "ENSP00000494188" in ret_by_name
    assert ret_by_id == ret_by_name


def test_protein_ids_from_protein_id(ensembl100):
    result = ensembl100.protein_ids("ENSP00000494188")
    assert isinstance(result, list)
    assert "ENSP00000494188" in result


def test_protein_ids_from_transcript(ensembl100):
    ret_by_id = ensembl100.protein_ids("ENST00000643777")
    ret_by_name = ensembl100.protein_ids("ALDOA-219")
    assert "ENSP00000494188" in ret_by_id
    assert "ENSP00000494188" in ret_by_name
    assert ret_by_id == ret_by_name


def test_transcript_ids_from_contig_id(ensembl100):
    result = ensembl100.transcript_ids("16")
    assert isinstance(result, list)
    assert "ENST00000643777" in result


def test_transcript_ids_from_exon_id(ensembl100):
    result = ensembl100.transcript_ids("ENSE00003826864")
    assert isinstance(result, list)
    assert "ENST00000643777" in result


def test_transcript_ids_from_gene(ensembl100):
    ret_by_id = ensembl100.transcript_ids("ENSG00000149925")
    ret_by_name = ensembl100.transcript_ids("ALDOA")
    assert "ENST00000643777" in ret_by_id
    assert "ENST00000643777" in ret_by_name
    assert ret_by_id == ret_by_name


def test_transcript_ids_from_protein_id(ensembl100):
    result = ensembl100.transcript_ids("ENSP00000494188")
    assert isinstance(result, list)
    assert "ENST00000643777" in result


def test_transcript_ids_from_transcript(ensembl100):
    ret_by_id = ensembl100.transcript_ids("ENST00000643777")
    ret_by_name = ensembl100.transcript_ids("ALDOA-219")
    assert "ENST00000643777" in ret_by_id
    assert "ENST00000643777" in ret_by_name
    assert ret_by_id == ret_by_name


def test_transcript_names_from_contig_id(ensembl100):
    result = ensembl100.transcript_names("16")
    assert isinstance(result, list)
    assert "ALDOA-219" in result


def test_transcript_names_from_exon_id(ensembl100):
    result = ensembl100.transcript_names("ENSE00003826864")
    assert isinstance(result, list)
    assert "ALDOA-219" in result


def test_transcript_names_from_gene(ensembl100):
    ret_by_id = ensembl100.transcript_names("ENSG00000149925")
    ret_by_name = ensembl100.transcript_names("ALDOA")
    assert "ALDOA-219" in ret_by_id
    assert "ALDOA-219" in ret_by_name
    assert ret_by_id == ret_by_name


def test_transcript_names_from_protein_id(ensembl100):
    result = ensembl100.transcript_names("ENSP00000494188")
    assert isinstance(result, list)
    assert "ALDOA-219" in result


def test_transcript_names_from_transcript(ensembl100):
    ret_by_id = ensembl100.transcript_names("ENST00000643777")
    ret_by_name = ensembl100.transcript_names("ALDOA-219")
    assert "ALDOA-219" in ret_by_id
    assert "ALDOA-219" in ret_by_name
    assert ret_by_id == ret_by_name


# -------------------------------------------------------------------------------------------------
# test is_canonical_transcript
# -------------------------------------------------------------------------------------------------
def test_is_canonical_transcript(ensembl100):
    assert ensembl100.is_canonical_transcript("ENST00000000233") is True
    assert ensembl100.is_canonical_transcript("spam") is False


# -------------------------------------------------------------------------------------------------
# test is_<feature>
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
    assert ensembl100.is_exon("ENSE00003826864") is True


def test_is_gene_false(ensembl100):
    assert ensembl100.is_gene("spam") is False


def test_is_gene_true_from_id(ensembl100):
    assert ensembl100.is_gene("ENSG00000149925") is True


def test_is_gene_true_from_name(ensembl100):
    assert ensembl100.is_gene("ALDOA") is True


def test_is_protein_false(ensembl100):
    assert ensembl100.is_protein("spam") is False


def test_is_protein_true(ensembl100):
    assert ensembl100.is_protein("ENSP00000494188") is True


def test_is_transcript_false(ensembl100):
    assert ensembl100.is_transcript("spam") is False


def test_is_transcript_true_from_id(ensembl100):
    assert ensembl100.is_transcript("ENST00000643777") is True


def test_is_transcript_true_from_name(ensembl100):
    assert ensembl100.is_transcript("ALDOA-219") is True


# -------------------------------------------------------------------------------------------------
# test <feature>_sequence
# -------------------------------------------------------------------------------------------------
@pytest.fixture
def test_cds_sequence_pos():
    return CdnaPosition(
        "4",
        7,
        0,
        9,
        0,
        "+",
        "ENSG00000157404",
        "KIT",
        "ENST00000288135",
        "KIT-201",
        "ENSP00000288135",
    )


def test_cds_sequence(ensembl100, test_cds_sequence_pos):
    assert ensembl100.sequence(test_cds_sequence_pos) == "GGC"


@pytest.fixture
def test_dna_sequence_pos_plus():
    return DnaPosition("4", 54695512, 0, 54695514, 0, "+")


@pytest.fixture
def test_dna_sequence_pos_minus():
    return DnaPosition("4", 54695512, 0, 54695514, 0, "-")


def test_dna_sequence(ensembl100, test_dna_sequence_pos_plus, test_dna_sequence_pos_minus):
    assert ensembl100.sequence(test_dna_sequence_pos_plus) == "GCT"
    assert ensembl100.sequence(test_dna_sequence_pos_minus) == "AGC"


@pytest.fixture
def test_protein_sequence_pos():
    return ProteinPosition(
        "4",
        3,
        0,
        4,
        0,
        "+",
        "ENSG00000157404",
        "KIT",
        "ENST00000288135",
        "KIT-201",
        "ENSP00000288135",
    )


def test_protein_sequence(ensembl100, test_protein_sequence_pos):
    assert ensembl100.sequence(test_protein_sequence_pos) == "GA"


@pytest.fixture
def test_rna_sequence_pos():
    return RnaPosition(
        "4", 126, 0, 128, 0, "+", "ENSG00000157404", "KIT", "ENST00000288135", "KIT-201"
    )


def test_rna_sequence(ensembl100, test_rna_sequence_pos):
    assert ensembl100.sequence(test_rna_sequence_pos) == "GCT"


# -------------------------------------------------------------------------------------------------
# test normalize_id(feature, feature_type)
# -------------------------------------------------------------------------------------------------
def test_normalize_id_contig_id(ensembl100):
    assert ensembl100.normalize_id("4") == [("4", CONTIG_ID)]


def test_normalize_id_exon_id(ensembl100):
    assert ensembl100.normalize_id("ENSE00003826864") == [("ENSE00003826864", EXON_ID)]


def test_normalize_id_gene_id(ensembl100):
    assert ensembl100.normalize_id("ENSG00000149925") == [("ENSG00000149925", GENE_ID)]


def test_normalize_id_gene_name(ensembl100):
    assert ensembl100.normalize_id("ALDOA") == [("ALDOA", GENE_NAME)]


def test_normalize_id_protein_id(ensembl100):
    assert ensembl100.normalize_id("ENSP00000494188") == [("ENSP00000494188", PROTEIN_ID)]


def test_normalize_id_refseq_transcript_id(ensembl100):
    assert ensembl100.normalize_id("NM_000314.4") == [("ENST00000371953", TRANSCRIPT_ID)]


def test_normalize_id_transcript_id(ensembl100):
    assert ensembl100.normalize_id("ENST00000643777") == [("ENST00000643777", TRANSCRIPT_ID)]


def test_normalize_id_transcript_name(ensembl100):
    assert ensembl100.normalize_id("ALDOA-219") == [("ALDOA-219", TRANSCRIPT_NAME)]


# -------------------------------------------------------------------------------------------------
# test str
# -------------------------------------------------------------------------------------------------
def test_cdna_to_str():
    position1 = CdnaPosition(
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


def test_dna_to_str():
    position1 = DnaPosition(
        contig_id="5", start=1234, start_offset=0, end=1234, end_offset=0, strand="-"
    )
    position2 = DnaPosition(
        contig_id="5", start=1, start_offset=0, end=1234, end_offset=0, strand="-"
    )
    assert str(position1) == "5:g.1234"
    assert str(position2) == "5:g.1_1234"


def test_exon_to_str():
    position1 = ExonPosition(
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


def test_protein_to_str():
    position1 = ProteinPosition(
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


def test_transcript_to_str():
    position1 = RnaPosition(
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
# test <feature>
# -------------------------------------------------------------------------------------------------
@pytest.fixture
def expected_cdna():
    return CdnaPosition(
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
        ("5", 4289),
        ("ENSE00003896691", 1),
        ("ENSG00000164362", 4),
        ("ENSP00000309572", 1),
        ("ENST00000310581", 1),
        ("TERT", 4),
        ("TERT-201", 1),
    ],
)
def test_cdna(ensembl100, expected_cdna, feature, num_results):
    results = ensembl100.cdna(feature)
    assert len(results) == num_results
    assert expected_cdna in results


@pytest.fixture
def expected_dna():
    return DnaPosition(
        contig_id="5", start=1, start_offset=0, end=181538259, end_offset=0, strand="-"
    )


@pytest.mark.parametrize(
    "feature,num_results",
    [
        ("5", 2),
        ("ENSE00003896691", 1),
        ("ENSG00000164362", 1),
        ("ENSP00000309572", 1),
        ("ENST00000310581", 1),
        ("TERT", 1),
        ("TERT-201", 1),
    ],
)
def test_dna(ensembl100, expected_dna, feature, num_results):
    results = ensembl100.dna(feature)
    assert len(results) == num_results
    assert expected_dna in results


@pytest.fixture
def expected_exon():
    return ExonPosition(
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
        ("5", 62673),
        ("ENSE00003896691", 1),
        ("ENSG00000164362", 79),
        ("ENSP00000309572", 33),
        ("ENST00000310581", 33),
        ("TERT", 79),
        ("TERT-201", 33),
    ],
)
def test_exon(ensembl100, expected_exon, feature, num_results):
    results = ensembl100.exon(feature)
    assert len(results) == num_results
    assert expected_exon in results


@pytest.fixture
def expected_gene():
    return DnaPosition(
        contig_id="5", start=1253147, start_offset=0, end=1295068, end_offset=0, strand="-"
    )


@pytest.mark.parametrize(
    "feature,num_results",
    [
        ("5", 2987),
        ("ENSE00003896691", 1),
        ("ENSG00000164362", 1),
        ("ENSP00000309572", 1),
        ("ENST00000310581", 1),
        ("TERT", 1),
        ("TERT-201", 1),
    ],
)
def test_gene(ensembl100, expected_gene, feature, num_results):
    results = ensembl100.gene(feature)
    assert len(results) == num_results
    assert expected_gene in results


@pytest.fixture
def expected_protein():
    return ProteinPosition(
        contig_id="5",
        start=1,
        start_offset=0,
        end=1133,
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
        ("5", 4289),
        ("ENSE00003896691", 1),
        ("ENSG00000164362", 4),
        ("ENSP00000309572", 1),
        ("ENST00000310581", 1),
        ("TERT", 4),
        ("TERT-201", 1),
    ],
)
def test_protein(ensembl100, expected_protein, feature, num_results):
    results = ensembl100.protein(feature)
    assert len(results) == num_results
    assert expected_protein in results


@pytest.fixture
def expected_rna():
    return RnaPosition(
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
        ("5", 10806),
        ("ENSE00003896691", 1),
        ("ENSG00000164362", 7),
        ("ENSP00000309572", 1),
        ("ENST00000310581", 1),
        ("TERT", 7),
        ("TERT-201", 1),
    ],
)
def test_transcript(ensembl100, expected_rna, feature, num_results):
    results = ensembl100.rna(feature)
    assert len(results) == num_results
    assert expected_rna in results


# -------------------------------------------------------------------------------------------------
# test <feature> canonical only
# -------------------------------------------------------------------------------------------------
def test_cdna_canonical(ensembl100):
    results = ensembl100.cdna("ARF5", canonical=True)
    assert len(results) == 1
    assert results[0].transcript_id == "ENST00000000233"


def test_exon_canonical(ensembl100):
    results = ensembl100.exon("ARF5", canonical=True)
    assert len(results) == 6
    for exon in results:
        assert exon.transcript_id == "ENST00000000233"


def test_protein_canonical(ensembl100):
    results = ensembl100.protein("ARF5", canonical=True)
    assert len(results) == 1
    assert results[0].transcript_id == "ENST00000000233"


def test_rna_canonical(ensembl100):
    results = ensembl100.rna("ARF5", canonical=True)
    assert len(results) == 1
    assert results[0].transcript_id == "ENST00000000233"


# -------------------------------------------------------------------------------------------------
# test translate_cdna_variant
# -------------------------------------------------------------------------------------------------
@pytest.fixture
def test_cds_sequence_pos_1():
    return CdnaPosition(
        "4",
        38,
        0,
        38,
        0,
        "+",
        "ENSG00000157404",
        "KIT",
        "ENST00000288135",
        "KIT-201",
        "ENSP00000288135",
    )


def test_translate_cdna_variant_1(ensembl100, test_cds_sequence_pos_1):
    # GTT -> GAT
    assert ensembl100.translate_cdna_variant(test_cds_sequence_pos_1, "A") == ["D"]


def test_translate_cdna_variant_3(ensembl100, test_cds_sequence_pos_1):
    # GTT -> GAT or GCT
    assert ensembl100.translate_cdna_variant(test_cds_sequence_pos_1, "M") == ["A", "D"]


@pytest.fixture
def test_cds_sequence_pos_2():
    return CdnaPosition(
        "4",
        38,
        0,
        40,
        0,
        "+",
        "ENSG00000157404",
        "KIT",
        "ENST00000288135",
        "KIT-201",
        "ENSP00000288135",
    )


def test_translate_cdna_variant_2(ensembl100, test_cds_sequence_pos_2):
    # GTTCTG -> GAAATG
    assert ensembl100.translate_cdna_variant(test_cds_sequence_pos_2, "AAA") == ["EM"]


@pytest.fixture
def test_cds_sequence_pos_3():
    return CdnaPosition(
        "4",
        1674,
        0,
        1675,
        0,
        "+",
        "ENSG00000157404",
        "KIT",
        "ENST00000288135",
        "KIT-201",
        "ENSP00000288135",
    )


def test_translate_cdna_variant_4(ensembl100, test_cds_sequence_pos_3):
    # GG -> GTTCG
    assert ensembl100.translate_cdna_variant(test_cds_sequence_pos_3, "GTTCG") == ["KFV"]


# -------------------------------------------------------------------------------------------------
# test load
# -------------------------------------------------------------------------------------------------
def test_load_cdna_minimal(ensembl100):
    result = ensembl100.variant(position_type="cdna", feature="ENST00000375562", start=123)
    assert result == [
        CdnaPosition(
            contig_id="6",
            start=123,
            start_offset=0,
            end=123,
            end_offset=0,
            strand="-",
            gene_id="ENSG00000204385",
            gene_name="SLC44A4",
            transcript_id="ENST00000375562",
            transcript_name="SLC44A4-203",
            protein_id="ENSP00000364712",
        )
    ]


def test_load_cdna_position(ensembl100):
    result = ensembl100.variant(
        position_type="cdna",
        feature="ENST00000288135",
        start=68,
        start_offset=-1,
        end=68,
        end_offset=-1,
        strand="+",
    )
    assert result == [
        CdnaPosition(
            contig_id="4",
            start=68,
            start_offset=-1,
            end=68,
            end_offset=-1,
            strand="+",
            gene_id="ENSG00000157404",
            gene_name="KIT",
            transcript_id="ENST00000288135",
            transcript_name="KIT-201",
            protein_id="ENSP00000288135",
        )
    ]


def test_load_cdna_substitution(ensembl100):
    result = ensembl100.variant(
        position_type="cdna",
        feature="ENST00000256078",
        start=38,
        start_offset=0,
        end=38,
        end_offset=0,
        strand="-",
        refseq="G",
        altseq="A",
        variant_type="substitution",
    )
    assert result == [
        CdnaSubstitution(
            contig_id="12",
            start=38,
            start_offset=0,
            end=38,
            end_offset=0,
            strand="-",
            gene_id="ENSG00000133703",
            gene_name="KRAS",
            transcript_id="ENST00000256078",
            transcript_name="KRAS-201",
            protein_id="ENSP00000256078",
            refseq="G",
            altseq="A",
        )
    ]


def test_load_dna_minimal(ensembl100):
    result = ensembl100.variant(position_type="dna", feature="5", start=1282623)
    assert result == [
        DnaPosition(
            contig_id="5", start=1282623, start_offset=0, end=1282623, end_offset=0, strand="+"
        ),
        DnaPosition(
            contig_id="5", start=1282623, start_offset=0, end=1282623, end_offset=0, strand="-"
        ),
    ]


def test_load_dna_position(ensembl100):
    result = ensembl100.variant(
        position_type="dna", feature="5", start=1282623, end=1282626, strand="-"
    )
    assert result == [
        DnaPosition(
            contig_id="5", start=1282623, start_offset=0, end=1282626, end_offset=0, strand="-"
        )
    ]


def test_load_dna_substitution(ensembl100):
    result = ensembl100.variant(
        position_type="dna",
        feature="5",
        start=1282623,
        end=1282626,
        strand="-",
        refseq="A",
        altseq="G",
        variant_type="substitution",
    )
    assert result == [
        DnaSubstitution(
            contig_id="5",
            start=1282623,
            start_offset=0,
            end=1282626,
            end_offset=0,
            strand="-",
            refseq="A",
            altseq="G",
        )
    ]


def test_load_exon_fusion(ensembl100):
    result = ensembl100.variant(
        position_type="exon",
        feature="ENST00000315869",
        start=3,
        end=4,
        strand="-",
        variant_type="fusion",
        position_type2="exon",
        feature2="ENST00000271526",
        start2=2,
        end2=2,
        strand2="+",
    )
    assert result == [
        ExonFusion(
            ExonPosition(
                contig_id="X",
                start=3,
                start_offset=0,
                end=4,
                end_offset=0,
                strand="-",
                gene_id="ENSG00000068323",
                gene_name="TFE3",
                transcript_id="ENST00000315869",
                transcript_name="TFE3-201",
                exon_id="ENSE00003528623",
            ),
            ExonPosition(
                contig_id="1",
                start=2,
                start_offset=0,
                end=2,
                end_offset=0,
                strand="+",
                gene_id="ENSG00000143294",
                gene_name="PRCC",
                transcript_id="ENST00000271526",
                transcript_name="PRCC-201",
                exon_id="ENSE00003525835",
            ),
        )
    ]
