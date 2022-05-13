import pytest

from ensembl_map.core import EnsemblRelease, instance
from ensembl_map.constants import CDS, CONTIG, EXON, GENE, PROTEIN, TRANSCRIPT

CACHE_DIR = "/home/matt/Downloads/ensembl_map_data"  # DEBUG

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
        species="homo_sapiens", reference="GRCh38", release=100, cache_dir=CACHE_DIR
    )


def test_contig_ids(ensembl100):
    result = ensembl100.contig_ids()
    assert isinstance(result, list)
    assert CONTIG_ID in result


def test_exon_ids(ensembl100):
    result = ensembl100.exon_ids()
    assert isinstance(result, list)
    assert EXON_ID in result


def test_gene_ids(ensembl100):
    result = ensembl100.gene_ids()
    assert isinstance(result, list)
    assert GENE_ID in result


def test_gene_names(ensembl100):
    result = ensembl100.gene_names()
    assert isinstance(result, list)
    assert GENE_NAME in result


def test_protein_ids(ensembl100):
    result = ensembl100.protein_ids()
    assert isinstance(result, list)
    assert PROTEIN_ID in result


def test_transcript_ids(ensembl100):
    result = ensembl100.transcript_ids()
    assert isinstance(result, list)
    assert TRANSCRIPT_ID in result


def test_transcript_names(ensembl100):
    result = ensembl100.transcript_names()
    assert isinstance(result, list)
    assert TRANSCRIPT_NAME in result


def test_contig_ids_of_contig_id(ensembl100):
    result = ensembl100.contig_ids_of_contig_id(CONTIG_ID)
    assert isinstance(result, list)
    assert CONTIG_ID in result


def test_contig_ids_of_exon_id(ensembl100):
    result = ensembl100.contig_ids_of_exon_id(EXON_ID)
    assert isinstance(result, list)
    assert CONTIG_ID in result


def test_contig_ids_of_gene_id(ensembl100):
    result = ensembl100.contig_ids_of_gene_id(GENE_ID)
    assert isinstance(result, list)
    assert CONTIG_ID in result


def test_contig_ids_of_gene_name(ensembl100):
    result = ensembl100.contig_ids_of_gene_name(GENE_NAME)
    assert isinstance(result, list)
    assert CONTIG_ID in result


def test_contig_ids_of_protein_id(ensembl100):
    result = ensembl100.contig_ids_of_protein_id(PROTEIN_ID)
    assert isinstance(result, list)
    assert CONTIG_ID in result


def test_contig_ids_of_transcript_id(ensembl100):
    result = ensembl100.contig_ids_of_transcript_id(TRANSCRIPT_ID)
    assert isinstance(result, list)
    assert CONTIG_ID in result


def test_contig_ids_of_transcript_name(ensembl100):
    result = ensembl100.contig_ids_of_transcript_name(TRANSCRIPT_NAME)
    assert isinstance(result, list)
    assert CONTIG_ID in result


def test_exon_ids_of_contig_id(ensembl100):
    result = ensembl100.exon_ids_of_contig_id(CONTIG_ID)
    assert isinstance(result, list)
    assert EXON_ID in result


def test_exon_ids_of_exon_id(ensembl100):
    result = ensembl100.exon_ids_of_exon_id(EXON_ID)
    assert isinstance(result, list)
    assert EXON_ID in result


def test_exon_ids_of_gene_id(ensembl100):
    result = ensembl100.exon_ids_of_gene_id(GENE_ID)
    assert isinstance(result, list)
    assert EXON_ID in result


def test_exon_ids_of_gene_name(ensembl100):
    result = ensembl100.exon_ids_of_gene_name(GENE_NAME)
    assert isinstance(result, list)
    assert EXON_ID in result


def test_exon_ids_of_protein_id(ensembl100):
    result = ensembl100.exon_ids_of_protein_id(PROTEIN_ID)
    assert isinstance(result, list)
    assert EXON_ID in result


def test_exon_ids_of_transcript_id(ensembl100):
    result = ensembl100.exon_ids_of_transcript_id(TRANSCRIPT_ID)
    assert isinstance(result, list)
    assert EXON_ID in result


def test_exon_ids_of_transcript_name(ensembl100):
    result = ensembl100.exon_ids_of_transcript_name(TRANSCRIPT_NAME)
    assert isinstance(result, list)
    assert EXON_ID in result


def test_gene_ids_of_contig_id(ensembl100):
    result = ensembl100.gene_ids_of_contig_id(CONTIG_ID)
    assert isinstance(result, list)
    assert GENE_ID in result


def test_gene_ids_of_exon_id(ensembl100):
    result = ensembl100.gene_ids_of_exon_id(EXON_ID)
    assert isinstance(result, list)
    assert GENE_ID in result


def test_gene_ids_of_gene_id(ensembl100):
    result = ensembl100.gene_ids_of_gene_id(GENE_ID)
    assert isinstance(result, list)
    assert GENE_ID in result


def test_gene_ids_of_gene_name(ensembl100):
    result = ensembl100.gene_ids_of_gene_name(GENE_NAME)
    assert isinstance(result, list)
    assert GENE_ID in result


def test_gene_ids_of_protein_id(ensembl100):
    result = ensembl100.gene_ids_of_protein_id(PROTEIN_ID)
    assert isinstance(result, list)
    assert GENE_ID in result


def test_gene_ids_of_transcript_id(ensembl100):
    result = ensembl100.gene_ids_of_transcript_id(TRANSCRIPT_ID)
    assert isinstance(result, list)
    assert GENE_ID in result


def test_gene_ids_of_transcript_name(ensembl100):
    result = ensembl100.gene_ids_of_transcript_name(TRANSCRIPT_NAME)
    assert isinstance(result, list)
    assert GENE_ID in result


def test_gene_names_of_contig_id(ensembl100):
    result = ensembl100.gene_names_of_contig_id(CONTIG_ID)
    assert isinstance(result, list)
    assert GENE_NAME in result


def test_gene_names_of_exon_id(ensembl100):
    result = ensembl100.gene_names_of_exon_id(EXON_ID)
    assert isinstance(result, list)
    assert GENE_NAME in result


def test_gene_names_of_gene_id(ensembl100):
    result = ensembl100.gene_names_of_gene_id(GENE_ID)
    assert isinstance(result, list)
    assert GENE_NAME in result


def test_gene_names_of_gene_name(ensembl100):
    result = ensembl100.gene_names_of_gene_name(GENE_NAME)
    assert isinstance(result, list)
    assert GENE_NAME in result


def test_gene_names_of_protein_id(ensembl100):
    result = ensembl100.gene_names_of_protein_id(PROTEIN_ID)
    assert isinstance(result, list)
    assert GENE_NAME in result


def test_gene_names_of_transcript_id(ensembl100):
    result = ensembl100.gene_names_of_transcript_id(TRANSCRIPT_ID)
    assert isinstance(result, list)
    assert GENE_NAME in result


def test_gene_names_of_transcript_name(ensembl100):
    result = ensembl100.gene_names_of_transcript_name(TRANSCRIPT_NAME)
    assert isinstance(result, list)
    assert GENE_NAME in result


def test_protein_ids_of_contig_id(ensembl100):
    result = ensembl100.protein_ids_of_contig_id(CONTIG_ID)
    assert isinstance(result, list)
    assert PROTEIN_ID in result


def test_protein_ids_of_exon_id(ensembl100):
    result = ensembl100.protein_ids_of_exon_id(EXON_ID)
    assert isinstance(result, list)
    assert PROTEIN_ID in result


def test_protein_ids_of_gene_id(ensembl100):
    result = ensembl100.protein_ids_of_gene_id(GENE_ID)
    assert isinstance(result, list)
    assert PROTEIN_ID in result


def test_protein_ids_of_gene_name(ensembl100):
    result = ensembl100.protein_ids_of_gene_name(GENE_NAME)
    assert isinstance(result, list)
    assert PROTEIN_ID in result


def test_protein_ids_of_protein_id(ensembl100):
    result = ensembl100.protein_ids_of_protein_id(PROTEIN_ID)
    assert isinstance(result, list)
    assert PROTEIN_ID in result


def test_protein_ids_of_transcript_id(ensembl100):
    result = ensembl100.protein_ids_of_transcript_id(TRANSCRIPT_ID)
    assert isinstance(result, list)
    assert PROTEIN_ID in result


def test_protein_ids_of_transcript_name(ensembl100):
    result = ensembl100.protein_ids_of_transcript_name(TRANSCRIPT_NAME)
    assert isinstance(result, list)
    assert PROTEIN_ID in result


def test_transcript_ids_of_contig_id(ensembl100):
    result = ensembl100.transcript_ids_of_contig_id(CONTIG_ID)
    assert isinstance(result, list)
    assert TRANSCRIPT_ID in result


def test_transcript_ids_of_exon_id(ensembl100):
    result = ensembl100.transcript_ids_of_exon_id(EXON_ID)
    assert isinstance(result, list)
    assert TRANSCRIPT_ID in result


def test_transcript_ids_of_gene_id(ensembl100):
    result = ensembl100.transcript_ids_of_gene_id(GENE_ID)
    assert isinstance(result, list)
    assert TRANSCRIPT_ID in result


def test_transcript_ids_of_gene_name(ensembl100):
    result = ensembl100.transcript_ids_of_gene_name(GENE_NAME)
    assert isinstance(result, list)
    assert TRANSCRIPT_ID in result


def test_transcript_ids_of_protein_id(ensembl100):
    result = ensembl100.transcript_ids_of_protein_id(PROTEIN_ID)
    assert isinstance(result, list)
    assert TRANSCRIPT_ID in result


def test_transcript_ids_of_transcript_id(ensembl100):
    result = ensembl100.transcript_ids_of_transcript_id(TRANSCRIPT_ID)
    assert isinstance(result, list)
    assert TRANSCRIPT_ID in result


def test_transcript_ids_of_transcript_name(ensembl100):
    result = ensembl100.transcript_ids_of_transcript_name(TRANSCRIPT_NAME)
    assert isinstance(result, list)
    assert TRANSCRIPT_ID in result


def test_transcript_names_of_contig_id(ensembl100):
    result = ensembl100.transcript_names_of_contig_id(CONTIG_ID)
    assert isinstance(result, list)
    assert TRANSCRIPT_NAME in result


def test_transcript_names_of_exon_id(ensembl100):
    result = ensembl100.transcript_names_of_exon_id(EXON_ID)
    assert isinstance(result, list)
    assert TRANSCRIPT_NAME in result


def test_transcript_names_of_gene_id(ensembl100):
    result = ensembl100.transcript_names_of_gene_id(GENE_ID)
    assert isinstance(result, list)
    assert TRANSCRIPT_NAME in result


def test_transcript_names_of_gene_name(ensembl100):
    result = ensembl100.transcript_names_of_gene_name(GENE_NAME)
    assert isinstance(result, list)
    assert TRANSCRIPT_NAME in result


def test_transcript_names_of_protein_id(ensembl100):
    result = ensembl100.transcript_names_of_protein_id(PROTEIN_ID)
    assert isinstance(result, list)
    assert TRANSCRIPT_NAME in result


def test_transcript_names_of_transcript_id(ensembl100):
    result = ensembl100.transcript_names_of_transcript_id(TRANSCRIPT_ID)
    assert isinstance(result, list)
    assert TRANSCRIPT_NAME in result


def test_transcript_names_of_transcript_name(ensembl100):
    result = ensembl100.transcript_names_of_transcript_name(TRANSCRIPT_NAME)
    assert isinstance(result, list)
    assert TRANSCRIPT_NAME in result


# def test_is_canonical_transcript_instance(ensembl100):
#     assert ensembl100.is_canonical_transcript("ENST00000000233") is True
#     assert ensembl100.is_canonical_transcript("spam") is False


# def test_is_canonical_transcript_class(ensembl100):
#     assert instance().is_canonical_transcript("ENST00000000233") is True
#     assert instance().is_canonical_transcript("spam") is False


# def test_normalize_contig_id_instance(ensembl100):
#     assert ensembl100.normalize_contig_id("chr1") == ["1"]
#     assert ensembl100.normalize_contig_id("spam") == []


# def test_normalize_contig_id_class(ensembl100):
#     assert instance().normalize_contig_id("chr1") == ["1"]
#     assert instance().normalize_contig_id("spam") == []


# def test_gene_to_gene_id_different_gene_names(ensembl69):
#     # Gene name changed bewteen 69 and 100
#     # ENSG00000118058 MLL
#     # ENSG00000118058 KMT2A
#     assert ensembl69.gene_to_gene_id("MLL") == "ENSG00000118058"
#     assert ensembl69.gene_to_gene_id("KMT2A") == "ENSG00000118058"


# def test_gene_to_gene_id_different_gene_ids(ensembl69, ensembl100):
#     # Gene ID changed bewteen 69 and 100
#     # ENSG00000166748 AGBL1
#     # ENSG00000273540 AGBL1
#     assert ensembl69.gene_to_gene_id("AGBL1") == "ENSG00000166748"
#     assert ensembl100.gene_to_gene_id("AGBL1") == "ENSG00000273540"


# def test_gene_to_gene_id_different_gene_ids_2(ensembl69, ensembl100):
#     # Gene ID changed bewteen 69 and 100
#     # ENSG00000204645 SSX4
#     # ENSG00000268009 SSX4
#     assert ensembl69.gene_to_gene_id("SSX4") == "ENSG00000204645"
#     assert ensembl100.gene_to_gene_id("SSX4") == "ENSG00000268009"


# def test_gene_to_gene_name_different_gene_names(ensembl69, ensembl100):
#     # Gene name changed bewteen 69 and 100
#     # ENSG00000118058 MLL
#     # ENSG00000118058 KMT2A
#     assert ensembl69.gene_to_gene_name("ENSG00000118058") == "MLL"
#     assert ensembl100.gene_to_gene_name("ENSG00000118058") == "KMT2A"


# def test_gene_to_gene_name_different_gene_ids(ensembl69, ensembl100):
#     # Gene ID changed bewteen 69 and 100
#     # ENSG00000204645 SSX4
#     # ENSG00000268009 SSX4
#     assert ensembl69.gene_to_gene_name("ENSG00000204645") == "SSX4"
#     assert ensembl100.gene_to_gene_name("ENSG00000268009") == "SSX4"


def test_get_contig_ids_from_cds(ensembl100):
    ret_by_id = ensembl100.get_contig_ids("ENST00000288135", CDS)
    ret_by_name = ensembl100.get_contig_ids("KIT-201", CDS)
    assert "4" in ret_by_id
    assert "4" in ret_by_name
    assert ret_by_id, ret_by_name


def test_get_contig_ids_from_contig(ensembl100):
    ret_by_id = ensembl100.get_contig_ids("12", CONTIG)
    assert "12" in ret_by_id


def test_get_contig_ids_from_contig_with_chr(ensembl100):
    ret_by_id = ensembl100.get_contig_ids("chr12", CONTIG)
    assert "12" in ret_by_id


def test_get_contig_ids_from_exon(ensembl100):
    ret_by_id = ensembl100.get_contig_ids("ENSE00001074448", EXON)
    assert "4" in ret_by_id


def test_get_contig_ids_from_gene(ensembl100):
    ret_by_id = ensembl100.get_contig_ids("ENSG00000157404", GENE)
    ret_by_name = ensembl100.get_contig_ids("KIT", GENE)
    assert "4" in ret_by_id
    assert "4" in ret_by_name
    assert ret_by_id, ret_by_name


def test_get_contig_ids_from_protein(ensembl100):
    ret_by_id = ensembl100.get_contig_ids("ENSP00000288135", PROTEIN)
    assert "4" in ret_by_id


def test_get_contig_ids_from_transcript(ensembl100):
    ret_by_id = ensembl100.get_contig_ids("ENST00000288135", TRANSCRIPT)
    ret_by_name = ensembl100.get_contig_ids("KIT-201", TRANSCRIPT)
    assert "4" in ret_by_id
    assert "4" in ret_by_name
    assert ret_by_id, ret_by_name


def test_get_exon_ids_from_cds(ensembl100):
    ret_by_id = ensembl100.get_exon_ids("ENST00000288135", CDS)
    ret_by_name = ensembl100.get_exon_ids("KIT-201", CDS)
    assert "ENSE00001074448" in ret_by_id
    assert "ENSE00001074448" in ret_by_name
    assert ret_by_id, ret_by_name


def test_get_exon_ids_from_exon(ensembl100):
    ret_by_id = ensembl100.get_exon_ids("ENSE00001074448", EXON)
    assert "ENSE00001074448" in ret_by_id


def test_get_exon_ids_from_gene(ensembl100):
    ret_by_id = ensembl100.get_exon_ids("ENSG00000157404", GENE)
    ret_by_name = ensembl100.get_exon_ids("KIT", GENE)
    assert "ENSE00001074448" in ret_by_id
    assert "ENSE00001074448" in ret_by_name
    assert ret_by_id, ret_by_name


def test_get_exon_ids_from_protein(ensembl100):
    ret_by_id = ensembl100.get_exon_ids("ENSP00000288135", PROTEIN)
    assert "ENSE00001074448" in ret_by_id


def test_get_exon_ids_from_transcript(ensembl100):
    ret_by_id = ensembl100.get_exon_ids("ENST00000288135", TRANSCRIPT)
    ret_by_name = ensembl100.get_exon_ids("KIT-201", TRANSCRIPT)
    assert "ENSE00001074448" in ret_by_id
    assert "ENSE00001074448" in ret_by_name
    assert ret_by_id, ret_by_name


def test_exon_info_from_cds(ensembl100):
    ret_by_id = ensembl100.exon_info("ENST00000288135", CDS)
    ret_by_name = ensembl100.exon_info("KIT-201", CDS)
    assert "ENSE00001074448" in [i.exon_id for i in ret_by_id]
    assert "ENSE00001074448" in [i.exon_id for i in ret_by_name]
    assert ret_by_id, ret_by_name


def test_exon_info_from_exon(ensembl100):
    ret_by_id = ensembl100.exon_info("ENSE00001074448", EXON)
    assert "ENSE00001074448" in [i.exon_id for i in ret_by_id]


def test_exon_info_from_gene(ensembl100):
    ret_by_id = ensembl100.exon_info("ENSG00000157404", GENE)
    ret_by_name = ensembl100.exon_info("KIT", GENE)
    assert "ENSE00001074448" in [i.exon_id for i in ret_by_id]
    assert "ENSE00001074448" in [i.exon_id for i in ret_by_name]
    assert ret_by_id, ret_by_name


def test_exon_info_from_protein(ensembl100):
    ret_by_id = ensembl100.exon_info("ENSP00000288135", PROTEIN)
    assert "ENSE00001074448" in [i.exon_id for i in ret_by_id]


def test_exon_info_from_transcript(ensembl100):
    ret_by_id = ensembl100.exon_info("ENST00000288135", TRANSCRIPT)
    ret_by_name = ensembl100.exon_info("KIT-201", TRANSCRIPT)
    assert "ENSE00001074448" in [i.exon_id for i in ret_by_id]
    assert "ENSE00001074448" in [i.exon_id for i in ret_by_name]
    assert ret_by_id, ret_by_name


def test_get_gene_ids_from_cds(ensembl100):
    ret_by_id = ensembl100.get_gene_ids("ENST00000288135", CDS)
    ret_by_name = ensembl100.get_gene_ids("KIT-201", CDS)
    assert "ENSG00000157404" in ret_by_id
    assert "ENSG00000157404" in ret_by_name
    assert ret_by_id, ret_by_name


def test_get_gene_ids_from_exon(ensembl100):
    ret_by_id = ensembl100.get_gene_ids("ENSE00001074448", EXON)
    assert "ENSG00000157404" in ret_by_id


def test_get_gene_ids_from_gene(ensembl100):
    ret_by_id = ensembl100.get_gene_ids("ENSG00000157404", GENE)
    ret_by_name = ensembl100.get_gene_ids("KIT", GENE)
    assert "ENSG00000157404" in ret_by_id
    assert "ENSG00000157404" in ret_by_name
    assert ret_by_id, ret_by_name


def test_get_gene_ids_from_protein(ensembl100):
    ret_by_id = ensembl100.get_gene_ids("ENSP00000288135", PROTEIN)
    assert "ENSG00000157404" in ret_by_id


def test_get_gene_ids_from_transcript(ensembl100):
    ret_by_id = ensembl100.get_gene_ids("ENST00000288135", TRANSCRIPT)
    ret_by_name = ensembl100.get_gene_ids("KIT-201", TRANSCRIPT)
    assert "ENSG00000157404" in ret_by_id
    assert "ENSG00000157404" in ret_by_name
    assert ret_by_id, ret_by_name


def test_get_names_from_cds(ensembl100):
    ret_by_id = ensembl100.get_gene_names("ENST00000288135", CDS)
    ret_by_name = ensembl100.get_gene_names("KIT-201", CDS)
    assert "KIT" in ret_by_id
    assert "KIT" in ret_by_name
    assert ret_by_id, ret_by_name


def test_get_names_from_exon(ensembl100):
    ret_by_id = ensembl100.get_gene_names("ENSE00001074448", EXON)
    assert "KIT" in ret_by_id


def test_get_names_from_gene(ensembl100):
    ret_by_id = ensembl100.get_gene_names("ENSG00000157404", GENE)
    ret_by_name = ensembl100.get_gene_names("KIT", GENE)
    assert "KIT" in ret_by_id
    assert "KIT" in ret_by_name
    assert ret_by_id, ret_by_name


def test_get_names_from_protein(ensembl100):
    ret_by_id = ensembl100.get_gene_names("ENSP00000288135", PROTEIN)
    assert "KIT" in ret_by_id


def test_get_names_from_transcript(ensembl100):
    ret_by_id = ensembl100.get_gene_names("ENST00000288135", TRANSCRIPT)
    ret_by_name = ensembl100.get_gene_names("KIT-201", TRANSCRIPT)
    assert "KIT" in ret_by_id
    assert "KIT" in ret_by_name
    assert ret_by_id, ret_by_name


def test_get_names_from_refseq_transcript(ensembl100):
    # from SDEV-2909:
    ret_by_id = ensembl100.get_gene_names("NM_000244.3", TRANSCRIPT)
    assert "MEN1" in ret_by_id
    ret_by_id = ensembl100.get_gene_names("NM_007194.3", TRANSCRIPT)
    assert "CHEK2" in ret_by_id
    ret_by_id = ensembl100.get_gene_names("NM_001128849.1", TRANSCRIPT)
    assert "SMARCA4" in ret_by_id
    ret_by_id = ensembl100.get_gene_names("NM_000314.4", TRANSCRIPT)
    assert "PTEN" in ret_by_id


def test_contig_info_of_contig_id(ensembl100):
    ret_by_id = ensembl100.contig_info_of_contig_id("1")
    assert ret_by_id == [{"contig_id": "1", "start": 1, "end": 1}]


def test_gene_info_from_cds(ensembl100):
    ret_by_id = ensembl100.gene_info("ENST00000288135", CDS)
    ret_by_name = ensembl100.gene_info("KIT-201", CDS)
    assert "ENSG00000157404" in [i.gene_id for i in ret_by_id]
    assert "ENSG00000157404" in [i.gene_id for i in ret_by_name]
    assert ret_by_id, ret_by_name


def test_gene_info_from_exon(ensembl100):
    ret_by_id = ensembl100.gene_info("ENSE00001074448", EXON)
    assert "ENSG00000157404" in [i.gene_id for i in ret_by_id]


def test_gene_info_from_gene(ensembl100):
    ret_by_id = ensembl100.gene_info("ENSG00000157404", GENE)
    ret_by_name = ensembl100.gene_info("KIT", GENE)
    assert "ENSG00000157404" in [i.gene_id for i in ret_by_id]
    assert "ENSG00000157404" in [i.gene_id for i in ret_by_name]
    assert ret_by_id, ret_by_name


def test_gene_info_from_gene_old_name(ensembl100):
    # FAM175A was renamed to ABRAXAS1 in HG38
    ret_by_id = ensembl100.gene_info("ENSG00000163322", GENE)
    ret_by_name = ensembl100.gene_info("FAM175A", GENE)
    assert "ABRAXAS1" in [i.gene_name for i in ret_by_id]
    assert "ABRAXAS1" in [i.gene_name for i in ret_by_name]
    assert ret_by_id, ret_by_name


def test_gene_info_from_protein(ensembl100):
    ret_by_id = ensembl100.gene_info("ENSP00000288135", PROTEIN)
    assert "ENSG00000157404" in [i.gene_id for i in ret_by_id]


def test_gene_info_from_transcript(ensembl100):
    ret_by_id = ensembl100.gene_info("ENST00000288135", TRANSCRIPT)
    ret_by_name = ensembl100.gene_info("KIT-201", TRANSCRIPT)
    assert "ENSG00000157404" in [i.gene_id for i in ret_by_id]
    assert "ENSG00000157404" in [i.gene_id for i in ret_by_name]
    assert ret_by_id, ret_by_name


def test_is_contig_false(ensembl100):
    assert ensembl100.is_contig("spam") is False


def test_is_contig_true(ensembl100):
    assert ensembl100.is_contig("2") is True


def test_is_contig_true_alias(ensembl100):
    assert ensembl100.is_contig("chr2") is True


def test_is_exon_false(ensembl100):
    assert ensembl100.is_exon("spam") is False


def test_is_exon_true(ensembl100):
    assert ensembl100.is_exon("ENSE00003659301") is True


def test_is_gene_false(ensembl100):
    assert ensembl100.is_gene("spam") is False


def test_is_gene_true_from_id(ensembl100):
    assert ensembl100.is_gene("ENSG00000157404") is True


def test_is_gene_true_from_name(ensembl100):
    assert ensembl100.is_gene("KIT") is True


def test_is_protein_false(ensembl100):
    assert ensembl100.is_protein("spam") is False


def test_is_protein_true(ensembl100):
    assert ensembl100.is_protein("ENSP00000288135") is True


def test_is_transcript_false(ensembl100):
    assert ensembl100.is_transcript("spam") is False


def test_is_transcript_true_from_id(ensembl100):
    assert ensembl100.is_transcript("ENST00000288135") is True


def test_is_transcript_true_from_name(ensembl100):
    assert ensembl100.is_transcript("KIT-201") is True


def test_normalize_feature_contig_id(ensembl100):
    assert ensembl100.normalize_feature("17") == [("17", CONTIG)]


def test_normalize_feature_exon_id(ensembl100):
    assert ensembl100.normalize_feature("ENSE00003659301") == [("ENSE00003659301", EXON)]


def test_normalize_feature_gene_id(ensembl100):
    assert ensembl100.normalize_feature("ENSG00000157404") == [("ENSG00000157404", GENE)]


def test_normalize_feature_hene_name(ensembl100):
    assert ensembl100.normalize_feature("KIT") == [("KIT", GENE)]


def test_normalize_feature_protein_id(ensembl100):
    assert ensembl100.normalize_feature("ENSP00000288135") == [("ENSP00000288135", PROTEIN)]


def test_normalize_feature_transcript_id(ensembl100):
    assert ensembl100.normalize_feature("ENST00000288135") == [("ENST00000288135", TRANSCRIPT)]


def test_normalize_feature_transcript_name(ensembl100):
    assert ensembl100.normalize_feature("KIT-201") == [("KIT-201", TRANSCRIPT)]


def test_normalize_feature_refseq_transcript_id(ensembl100):
    assert ensembl100.normalize_feature("NM_000314.4") == [
        ("ENST00000371953", TRANSCRIPT),
        ("ENST00000645317", TRANSCRIPT),
    ]


def test_get_protein_ids_from_cds(ensembl100):
    ret_by_id = ensembl100.get_protein_ids("ENST00000288135", CDS)
    ret_by_name = ensembl100.get_protein_ids("KIT-201", CDS)
    assert "ENSP00000288135" in ret_by_id
    assert "ENSP00000288135" in ret_by_name
    assert ret_by_id, ret_by_name


def test_get_protein_ids_from_exon(ensembl100):
    ret_by_id = ensembl100.get_protein_ids("ENSE00001074448", EXON)
    assert "ENSP00000288135" in ret_by_id


def test_get_protein_ids_from_gene(ensembl100):
    ret_by_id = ensembl100.get_protein_ids("ENSG00000157404", GENE)
    ret_by_name = ensembl100.get_protein_ids("KIT", GENE)
    assert "ENSP00000288135" in ret_by_id
    assert "ENSP00000288135" in ret_by_name
    assert ret_by_id, ret_by_name


def test_get_protein_ids_from_protein(ensembl100):
    ret_by_id = ensembl100.get_protein_ids("ENSP00000288135", PROTEIN)
    assert "ENSP00000288135" in ret_by_id


def test_get_protein_ids_from_transcript(ensembl100):
    ret_by_id = ensembl100.get_protein_ids("ENST00000288135", TRANSCRIPT)
    ret_by_name = ensembl100.get_protein_ids("KIT-201", TRANSCRIPT)
    assert "ENSP00000288135" in ret_by_id
    assert "ENSP00000288135" in ret_by_name
    assert ret_by_id, ret_by_name


def test_get_transcript_ids_from_cds(ensembl100):
    ret_by_id = ensembl100.get_transcript_ids("ENST00000288135", CDS)
    ret_by_name = ensembl100.get_transcript_ids("KIT-201", CDS)
    assert "ENST00000288135" in ret_by_id
    assert "ENST00000288135" in ret_by_name
    assert ret_by_id, ret_by_name


def test_get_transcript_ids_from_exon(ensembl100):
    ret_by_id = ensembl100.get_transcript_ids("ENSE00001074448", EXON)
    assert "ENST00000288135" in ret_by_id


def test_get_transcript_ids_from_gene(ensembl100):
    ret_by_id = ensembl100.get_transcript_ids("ENSG00000157404", GENE)
    ret_by_name = ensembl100.get_transcript_ids("KIT", GENE)
    assert "ENST00000288135" in ret_by_id
    assert "ENST00000288135" in ret_by_name
    assert ret_by_id, ret_by_name


def test_get_transcript_ids_from_protein(ensembl100):
    ret_by_id = ensembl100.get_transcript_ids("ENSP00000288135", PROTEIN)
    assert "ENST00000288135" in ret_by_id


def test_get_transcript_ids_from_transcript(ensembl100):
    ret_by_id = ensembl100.get_transcript_ids("ENST00000288135", TRANSCRIPT)
    ret_by_name = ensembl100.get_transcript_ids("KIT-201", TRANSCRIPT)
    assert "ENST00000288135" in ret_by_id
    assert "ENST00000288135" in ret_by_name
    assert ret_by_id, ret_by_name


def test_get_transcript_ids_from_gene_best_only(ensembl100):
    # assumes ENST00000288135 is the best transcript for ENSG00000157404 (KIT)
    ret_all = ensembl100.get_transcript_ids("ENSG00000157404", GENE, best_only=False)
    ret_best = ensembl100.get_transcript_ids("ENSG00000157404", GENE, best_only=True)
    assert "ENST00000412167" in ret_all
    assert "ENST00000412167" not in ret_best


def test_get_transcript_names_from_cds(ensembl100):
    ret_by_id = ensembl100.get_transcript_names("ENST00000288135", CDS)
    ret_by_name = ensembl100.get_transcript_names("KIT-201", CDS)
    assert "KIT-201" in ret_by_id
    assert "KIT-201" in ret_by_name
    assert ret_by_id, ret_by_name


def test_get_transcript_names_from_exon(ensembl100):
    ret_by_id = ensembl100.get_transcript_names("ENSE00001074448", EXON)
    assert "KIT-201" in ret_by_id


def test_get_transcript_names_from_gene(ensembl100):
    ret_by_id = ensembl100.get_transcript_names("ENSG00000157404", GENE)
    ret_by_name = ensembl100.get_transcript_names("KIT", GENE)
    assert "KIT-201" in ret_by_id
    assert "KIT-201" in ret_by_name
    assert ret_by_id, ret_by_name


def test_get_transcript_names_from_protein(ensembl100):
    ret_by_id = ensembl100.get_transcript_names("ENSP00000288135", PROTEIN)
    assert "KIT-201" in ret_by_id


def test_get_transcript_names_from_transcript(ensembl100):
    ret_by_id = ensembl100.get_transcript_names("ENST00000288135", TRANSCRIPT)
    ret_by_name = ensembl100.get_transcript_names("KIT-201", TRANSCRIPT)
    assert "KIT-201" in ret_by_id
    assert "KIT-201" in ret_by_name
    assert ret_by_id, ret_by_name


def test_transcript_info_from_cds(ensembl100):
    ret_by_id = ensembl100.transcript_info("ENST00000288135", CDS)
    ret_by_name = ensembl100.transcript_info("KIT-201", CDS)
    assert "ENST00000288135" in [i.transcript_id for i in ret_by_id]
    assert "ENST00000288135" in [i.transcript_id for i in ret_by_name]
    assert ret_by_id, ret_by_name


def test_transcript_info_from_exon(ensembl100):
    ret_by_id = ensembl100.transcript_info("ENSE00003659301", EXON)
    assert "ENST00000380152" in [i.transcript_id for i in ret_by_id]


def test_transcript_info_from_gene(ensembl100):
    ret_by_id = ensembl100.transcript_info("ENSG00000157404", GENE)
    ret_by_name = ensembl100.transcript_info("KIT", GENE)
    assert "ENST00000288135" in [i.transcript_id for i in ret_by_id]
    assert "ENST00000288135" in [i.transcript_id for i in ret_by_name]
    assert ret_by_id, ret_by_name


def test_transcript_info_from_protein(ensembl100):
    ret_by_id = ensembl100.transcript_info("ENSP00000288135", PROTEIN)
    assert "ENST00000288135" in [i.transcript_id for i in ret_by_id]


def test_transcript_info_from_transcript(ensembl100):
    ret_by_id = ensembl100.transcript_info("ENST00000288135", TRANSCRIPT)
    ret_by_name = ensembl100.transcript_info("KIT-201", TRANSCRIPT)
    assert "ENST00000288135" in [i.transcript_id for i in ret_by_id]
    assert "ENST00000288135" in [i.transcript_id for i in ret_by_name]
    assert ret_by_id, ret_by_name
