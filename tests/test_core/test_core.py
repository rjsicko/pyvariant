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
from variant_map.core import (
    CdnaMappablePosition,
    Core,
    DnaMappablePosition,
    ExonMappablePosition,
    ProteinMappablePosition,
    RnaMappablePosition,
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
def test_cds_sequence(ensembl100):
    assert ensembl100.cds_sequence("ENST00000288135", 7, 9) == "GGC"


def test_dna_sequence(ensembl100):
    assert ensembl100.dna_sequence("4", 54695512, 54695514, "+") == "GCT"
    assert ensembl100.dna_sequence("4", 54695512, 54695514, "-") == "AGC"


def test_protein_sequence(ensembl100):
    assert ensembl100.protein_sequence("ENSP00000288135", 3, 4) == "GA"


def test_rna_sequence(ensembl100):
    assert ensembl100.rna_sequence("ENST00000288135", 126, 128) == "GCT"


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
def test_cdna_to_str(ensembl100):
    position1 = CdnaMappablePosition(
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
    position2 = CdnaMappablePosition(
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
    position1 = DnaMappablePosition(
        _data=ensembl100,
        contig_id="5",
        start=1234,
        start_offset=0,
        end=1234,
        end_offset=0,
        strand="-",
    )
    position2 = DnaMappablePosition(
        _data=ensembl100, contig_id="5", start=1, start_offset=0, end=1234, end_offset=0, strand="-"
    )
    assert str(position1) == "5:g.1234"
    assert str(position2) == "5:g.1_1234"


def test_exon_to_str(ensembl100):
    position1 = ExonMappablePosition(
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
    position2 = ExonMappablePosition(
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
    position1 = ProteinMappablePosition(
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
    position2 = ProteinMappablePosition(
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
    position1 = RnaMappablePosition(
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
    position2 = RnaMappablePosition(
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
# test get_<feature>
# -------------------------------------------------------------------------------------------------
@pytest.fixture
def test_cdna(ensembl100):
    return CdnaMappablePosition(
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
        ("5", 4289),
        ("ENSE00003896691", 1),
        ("ENSG00000164362", 4),
        ("ENSP00000309572", 1),
        ("ENST00000310581", 1),
        ("TERT", 4),
        ("TERT-201", 1),
    ],
)
def test_get_cdna(ensembl100, test_cdna, feature, num_results):
    results = ensembl100.cdna(feature)
    assert len(results) == num_results
    assert test_cdna in results


@pytest.fixture
def test_dna(ensembl100):
    return DnaMappablePosition(
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
        ("5", 2),
        ("ENSE00003896691", 1),
        ("ENSG00000164362", 1),
        ("ENSP00000309572", 1),
        ("ENST00000310581", 1),
        ("TERT", 1),
        ("TERT-201", 1),
    ],
)
def test_get_dna(ensembl100, test_dna, feature, num_results):
    results = ensembl100.dna(feature)
    assert len(results) == num_results
    assert test_dna in results


@pytest.fixture
def test_exon(ensembl100):
    return ExonMappablePosition(
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
        ("5", 62673),
        ("ENSE00003896691", 1),
        ("ENSG00000164362", 79),
        ("ENSP00000309572", 33),
        ("ENST00000310581", 33),
        ("TERT", 79),
        ("TERT-201", 33),
    ],
)
def test_get_exon(ensembl100, test_exon, feature, num_results):
    results = ensembl100.exon(feature)
    assert len(results) == num_results
    assert test_exon in results


@pytest.fixture
def test_gene(ensembl100):
    return DnaMappablePosition(
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
        ("5", 2987),
        ("ENSE00003896691", 1),
        ("ENSG00000164362", 1),
        ("ENSP00000309572", 1),
        ("ENST00000310581", 1),
        ("TERT", 1),
        ("TERT-201", 1),
    ],
)
def test_get_gene(ensembl100, test_gene, feature, num_results):
    results = ensembl100.gene(feature)
    assert len(results) == num_results
    assert test_gene in results


@pytest.fixture
def test_transcript(ensembl100):
    return RnaMappablePosition(
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
        ("5", 10806),
        ("ENSE00003896691", 1),
        ("ENSG00000164362", 7),
        ("ENSP00000309572", 1),
        ("ENST00000310581", 1),
        ("TERT", 7),
        ("TERT-201", 1),
    ],
)
def test_get_transcript(ensembl100, test_transcript, feature, num_results):
    results = ensembl100.rna(feature)
    assert len(results) == num_results
    assert test_transcript in results


# -------------------------------------------------------------------------------------------------
# test <feature>_to_<feature>
# -------------------------------------------------------------------------------------------------
def test_cdna_to_cdna(ensembl100):
    result = ensembl100.cdna_to_cdna(
        "ENST00000310581", 100, end=100, start_offset=-1, end_offset=-1, strand="-"
    )
    expected = [
        CdnaMappablePosition(
            _data=ensembl100,
            contig_id="5",
            start=99,
            start_offset=0,
            end=99,
            end_offset=0,
            strand="-",
            gene_id="ENSG00000164362",
            gene_name="TERT",
            transcript_id="ENST00000310581",
            transcript_name="TERT-201",
            protein_id="ENSP00000309572",
        )
    ]
    assert result == expected


def test_cdna_to_dna(ensembl100):
    result = ensembl100.cdna_to_dna(
        "ENST00000288135", 68, end=68, start_offset=-1, end_offset=-1, strand="+"
    )
    expected = [
        DnaMappablePosition(
            _data=ensembl100,
            contig_id="4",
            start=54695511,
            start_offset=0,
            end=54695511,
            end_offset=0,
            strand="+",
        )
    ]
    assert result == expected


def test_cdna_to_exon(ensembl100):
    result = ensembl100.cdna_to_exon(
        "ENST00000380152", 10252, end=10254, start_offset=0, end_offset=0, strand="+"
    )
    expected = [
        ExonMappablePosition(
            _data=ensembl100,
            contig_id="13",
            start=27,
            start_offset=0,
            end=27,
            end_offset=0,
            strand="+",
            gene_id="ENSG00000139618",
            gene_name="BRCA2",
            transcript_id="ENST00000380152",
            transcript_name="BRCA2-201",
            exon_id="ENSE00003717596",
        )
    ]
    assert result == expected


def test_cdna_to_protein(ensembl100):
    result = ensembl100.cdna_to_protein(
        "ENST00000310581", 100, end=100, start_offset=-1, end_offset=-1, strand="-"
    )
    expected = [
        ProteinMappablePosition(
            _data=ensembl100,
            contig_id="5",
            start=33,
            start_offset=0,
            end=33,
            end_offset=0,
            strand="-",
            gene_id="ENSG00000164362",
            gene_name="TERT",
            transcript_id="ENST00000310581",
            transcript_name="TERT-201",
            protein_id="ENSP00000309572",
        )
    ]
    assert result == expected


def test_cdna_to_rna(ensembl100):
    result = ensembl100.cdna_to_rna(
        "ENST00000288135", 68, end=68, start_offset=-1, end_offset=-1, strand="+"
    )
    expected = [
        RnaMappablePosition(
            _data=ensembl100,
            contig_id="4",
            start=126,
            start_offset=-1,
            end=126,
            end_offset=-1,
            strand="+",
            gene_id="ENSG00000157404",
            gene_name="KIT",
            transcript_id="ENST00000288135",
            transcript_name="KIT-201",
        )
    ]
    assert result == expected


def test_dna_to_cdna(ensembl100):
    result = ensembl100.dna_to_cdna(
        "13", 32316526, end=32319077, start_offset=1, end_offset=1, strand="+"
    )
    expected = [
        CdnaMappablePosition(
            _data=ensembl100,
            contig_id="13",
            start=67,
            start_offset=0,
            end=69,
            end_offset=0,
            strand="+",
            gene_id="ENSG00000139618",
            gene_name="BRCA2",
            transcript_id="ENST00000380152",
            transcript_name="BRCA2-201",
            protein_id="ENSP00000369497",
        ),
        CdnaMappablePosition(
            _data=ensembl100,
            contig_id="13",
            start=67,
            start_offset=0,
            end=69,
            end_offset=0,
            strand="+",
            gene_id="ENSG00000139618",
            gene_name="BRCA2",
            transcript_id="ENST00000544455",
            transcript_name="BRCA2-206",
            protein_id="ENSP00000439902",
        ),
    ]
    assert result == expected


def test_dna_to_dna(ensembl100):
    result = ensembl100.dna_to_dna(
        "5", 1282624, end=1293314, start_offset=1, end_offset=1, strand="-"
    )
    expected = [
        DnaMappablePosition(
            _data=ensembl100,
            contig_id="5",
            start=1282623,
            start_offset=0,
            end=1293313,
            end_offset=0,
            strand="-",
        )
    ]
    assert result == expected


def test_dna_to_exon(ensembl100):
    result = ensembl100.dna_to_exon(
        "6", 31166302, end=31166304, start_offset=0, end_offset=0, strand="-"
    )
    expected = [
        ExonMappablePosition(
            _data=ensembl100,
            contig_id="6",
            start=1,
            start_offset=0,
            end=1,
            end_offset=0,
            strand="-",
            gene_id="ENSG00000204531",
            gene_name="POU5F1",
            transcript_id="ENST00000471529",
            transcript_name="POU5F1-204",
            exon_id="ENSE00002568331",
        ),
        ExonMappablePosition(
            _data=ensembl100,
            contig_id="6",
            start=1,
            start_offset=0,
            end=1,
            end_offset=0,
            strand="-",
            gene_id="ENSG00000204531",
            gene_name="POU5F1",
            transcript_id="ENST00000513407",
            transcript_name="POU5F1-206",
            exon_id="ENSE00002033137",
        ),
    ]
    assert result == expected


def test_dna_to_protein(ensembl100):
    result = ensembl100.dna_to_protein(
        "13", 32316526, end=32319077, start_offset=1, end_offset=1, strand="+"
    )
    expected = [
        ProteinMappablePosition(
            _data=ensembl100,
            contig_id="13",
            start=23,
            start_offset=0,
            end=23,
            end_offset=0,
            strand="+",
            gene_id="ENSG00000139618",
            gene_name="BRCA2",
            transcript_id="ENST00000380152",
            transcript_name="BRCA2-201",
            protein_id="ENSP00000369497",
        ),
        ProteinMappablePosition(
            _data=ensembl100,
            contig_id="13",
            start=23,
            start_offset=0,
            end=23,
            end_offset=0,
            strand="+",
            gene_id="ENSG00000139618",
            gene_name="BRCA2",
            transcript_id="ENST00000544455",
            transcript_name="BRCA2-206",
            protein_id="ENSP00000439902",
        ),
    ]
    assert result == expected


def test_dna_to_rna(ensembl100):
    result = ensembl100.dna_to_rna(
        "5", 1253148, end=1253150, start_offset=1, end_offset=1, strand="-"
    )
    expected = [
        RnaMappablePosition(
            _data=ensembl100,
            contig_id="5",
            start=2420,
            start_offset=0,
            end=2422,
            end_offset=0,
            strand="-",
            gene_id="ENSG00000164362",
            gene_name="TERT",
            transcript_id="ENST00000484238",
            transcript_name="TERT-204",
        )
    ]
    assert result == expected


def test_exon_to_cdna(ensembl100):
    result = ensembl100.exon_to_cdna(
        "ENST00000380152", 27, end=27, start_offset=0, end_offset=0, strand="+"
    )
    expected = [
        CdnaMappablePosition(
            _data=ensembl100,
            contig_id="13",
            start=10255,
            start_offset=0,
            end=10257,
            end_offset=0,
            strand="+",
            gene_id="ENSG00000139618",
            gene_name="BRCA2",
            transcript_id="ENST00000380152",
            transcript_name="BRCA2-201",
            protein_id="ENSP00000369497",
        ),
        CdnaMappablePosition(
            _data=ensembl100,
            contig_id="13",
            start=9649,
            start_offset=0,
            end=10254,
            end_offset=0,
            strand="+",
            gene_id="ENSG00000139618",
            gene_name="BRCA2",
            transcript_id="ENST00000380152",
            transcript_name="BRCA2-201",
            protein_id="ENSP00000369497",
        ),
    ]
    assert result == expected


def test_exon_to_dna(ensembl100):
    result = ensembl100.exon_to_dna(
        "ENST00000310581", 16, end=16, start_offset=0, end_offset=0, strand="-"
    )
    expected = [
        DnaMappablePosition(
            _data=ensembl100,
            contig_id="5",
            start=1253167,
            start_offset=0,
            end=1253831,
            end_offset=0,
            strand="-",
        )
    ]
    assert result == expected


def test_exon_to_exon(ensembl100):
    result = ensembl100.exon_to_exon(
        "ENST00000380152", 27, end=27, start_offset=0, end_offset=0, strand="+"
    )
    expected = [
        ExonMappablePosition(
            _data=ensembl100,
            contig_id="13",
            start=27,
            start_offset=0,
            end=27,
            end_offset=0,
            strand="+",
            gene_id="ENSG00000139618",
            gene_name="BRCA2",
            transcript_id="ENST00000380152",
            transcript_name="BRCA2-201",
            exon_id="ENSE00003717596",
        )
    ]
    assert result == expected


def test_exon_to_protein(ensembl100):
    result = ensembl100.exon_to_protein(
        "ENST00000310581", 16, end=16, start_offset=0, end_offset=0, strand="-"
    )
    expected = [
        ProteinMappablePosition(
            _data=ensembl100,
            contig_id="5",
            start=1099,
            start_offset=0,
            end=1132,
            end_offset=0,
            strand="-",
            gene_id="ENSG00000164362",
            gene_name="TERT",
            transcript_id="ENST00000310581",
            transcript_name="TERT-201",
            protein_id="ENSP00000309572",
        )
    ]
    assert result == expected


def test_exon_to_rna(ensembl100):
    result = ensembl100.exon_to_rna(
        "ENST00000380152", 27, end=27, start_offset=0, end_offset=0, strand="+"
    )
    expected = [
        RnaMappablePosition(
            _data=ensembl100,
            contig_id="13",
            start=9882,
            start_offset=0,
            end=11986,
            end_offset=0,
            strand="+",
            gene_id="ENSG00000139618",
            gene_name="BRCA2",
            transcript_id="ENST00000380152",
            transcript_name="BRCA2-201",
        )
    ]
    assert result == expected


def test_protein_to_cdna(ensembl100):
    result = ensembl100.protein_to_cdna(
        "ENST00000310581", 525, end=525, start_offset=0, end_offset=0, strand="-"
    )
    expected = [
        CdnaMappablePosition(
            _data=ensembl100,
            contig_id="5",
            start=1573,
            start_offset=0,
            end=1575,
            end_offset=0,
            strand="-",
            gene_id="ENSG00000164362",
            gene_name="TERT",
            transcript_id="ENST00000310581",
            transcript_name="TERT-201",
            protein_id="ENSP00000309572",
        )
    ]
    assert result == expected


def test_protein_to_dna(ensembl100):
    result = ensembl100.protein_to_dna(
        "ENST00000310581", 1132, end=1132, start_offset=0, end_offset=0, strand="-"
    )
    expected = [
        DnaMappablePosition(
            _data=ensembl100,
            contig_id="5",
            start=1253731,
            start_offset=0,
            end=1253733,
            end_offset=0,
            strand="-",
        )
    ]
    assert result == expected


def test_protein_to_exon(ensembl100):
    result = ensembl100.protein_to_exon(
        "ENST00000380152", 3418, end=3418, start_offset=0, end_offset=0, strand="+"
    )
    to_positions = [
        ExonMappablePosition(
            _data=ensembl100,
            contig_id="13",
            start=27,
            start_offset=0,
            end=27,
            end_offset=0,
            strand="+",
            gene_id="ENSG00000139618",
            gene_name="BRCA2",
            transcript_id="ENST00000380152",
            transcript_name="BRCA2-201",
            exon_id="ENSE00003717596",
        )
    ]
    assert result == to_positions


def test_protein_to_protein(ensembl100):
    result = ensembl100.protein_to_protein(
        "ENSP00000309572", 525, end=525, start_offset=0, end_offset=0, strand="-"
    )
    to_positions = [
        ProteinMappablePosition(
            _data=ensembl100,
            contig_id="5",
            start=525,
            start_offset=0,
            end=525,
            end_offset=0,
            strand="-",
            gene_id="ENSG00000164362",
            gene_name="TERT",
            transcript_id="ENST00000310581",
            transcript_name="TERT-201",
            protein_id="ENSP00000309572",
        )
    ]
    assert result == to_positions


def test_protein_to_rna(ensembl100):
    result = ensembl100.protein_to_rna(
        "ENSP00000369497", 3418, end=3418, start_offset=0, end_offset=0, strand="+"
    )
    to_positions = [
        RnaMappablePosition(
            _data=ensembl100,
            contig_id="13",
            start=10485,
            start_offset=0,
            end=10487,
            end_offset=0,
            strand="+",
            gene_id="ENSG00000139618",
            gene_name="BRCA2",
            transcript_id="ENST00000380152",
            transcript_name="BRCA2-201",
        )
    ]
    assert result == to_positions


def test_rna_to_cdna(ensembl100):
    result = ensembl100.rna_to_cdna(
        "ENST00000310581", 1652, end=1654, start_offset=0, end_offset=0, strand="-"
    )
    to_positions = [
        CdnaMappablePosition(
            _data=ensembl100,
            contig_id="5",
            start=1573,
            start_offset=0,
            end=1575,
            end_offset=0,
            strand="-",
            gene_id="ENSG00000164362",
            gene_name="TERT",
            transcript_id="ENST00000310581",
            transcript_name="TERT-201",
            protein_id="ENSP00000309572",
        )
    ]
    assert result == to_positions


def test_rna_to_dna(ensembl100):
    result = ensembl100.rna_to_dna(
        "ENST00000310581", 1652, end=1654, start_offset=0, end_offset=0, strand="-"
    )
    to_positions = [
        DnaMappablePosition(
            _data=ensembl100,
            contig_id="5",
            start=1282623,
            start_offset=0,
            end=1293313,
            end_offset=0,
            strand="-",
        )
    ]
    assert result == to_positions


def test_rna_to_exon(ensembl100):
    result = ensembl100.rna_to_exon(
        "ENST00000310581", 4037, end=4039, start_offset=0, end_offset=0, strand="-"
    )
    to_positions = [
        ExonMappablePosition(
            _data=ensembl100,
            contig_id="5",
            start=16,
            start_offset=0,
            end=16,
            end_offset=0,
            strand="-",
            gene_id="ENSG00000164362",
            gene_name="TERT",
            transcript_id="ENST00000310581",
            transcript_name="TERT-201",
            exon_id="ENSE00001863787",
        )
    ]
    assert result == to_positions


def test_rna_to_protein(ensembl100):
    result = ensembl100.rna_to_protein(
        "ENST00000310581", 1652, end=1654, start_offset=0, end_offset=0, strand="-"
    )
    to_positions = [
        ProteinMappablePosition(
            _data=ensembl100,
            contig_id="5",
            start=525,
            start_offset=0,
            end=525,
            end_offset=0,
            strand="-",
            gene_id="ENSG00000164362",
            gene_name="TERT",
            transcript_id="ENST00000310581",
            transcript_name="TERT-201",
            protein_id="ENSP00000309572",
        )
    ]
    assert result == to_positions


def test_rna_to_rna(ensembl100):
    result = ensembl100.rna_to_rna(
        "ENST00000310581", 4037, end=4039, start_offset=0, end_offset=0, strand="-"
    )
    to_positions = [
        RnaMappablePosition(
            _data=ensembl100,
            contig_id="5",
            start=4037,
            start_offset=0,
            end=4039,
            end_offset=0,
            strand="-",
            gene_id="ENSG00000164362",
            gene_name="TERT",
            transcript_id="ENST00000310581",
            transcript_name="TERT-201",
        )
    ]
    assert result == to_positions


# -------------------------------------------------------------------------------------------------
# test translate_cds_variant
# -------------------------------------------------------------------------------------------------
def test_translate_cds_variant(ensembl100):
    # GTT -> GAT
    assert ensembl100.translate_cds_variant("ENST00000288135", 38, 38, "A") == ["D"]
    # GTTCTG -> GAAATG
    assert ensembl100.translate_cds_variant("ENST00000288135", 38, 40, "AAA") == ["EM"]
    # GTT -> GAT or GCT
    assert ensembl100.translate_cds_variant("ENST00000288135", 38, 38, "M") == ["A", "D"]
    # GG -> GTTCG
    assert ensembl100.translate_cds_variant("ENST00000288135", 1674, 1675, "GTTCG") == ["KFV"]
