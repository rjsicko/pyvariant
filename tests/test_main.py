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
    CdnaMappablePosition,
    DnaMappablePosition,
    ExonMappablePosition,
    ProteinMappablePosition,
    RnaMappablePosition,
)
from ensembl_map.main import (
    all_contig_ids,
    all_exon_ids,
    all_gene_ids,
    all_gene_names,
    all_protein_ids,
    all_transcript_ids,
    all_transcript_names,
    cdna_to_cdna,
    cdna_to_dna,
    cdna_to_exon,
    cdna_to_protein,
    cdna_to_rna,
    cds_sequence,
    contig_ids,
    dna_sequence,
    dna_to_cdna,
    dna_to_dna,
    dna_to_exon,
    dna_to_protein,
    dna_to_rna,
    exon_ids,
    exon_to_cdna,
    exon_to_dna,
    exon_to_exon,
    exon_to_protein,
    exon_to_rna,
    gene_ids,
    gene_names,
    get_cdna,
    get_dna,
    get_exons,
    get_genes,
    get_transcripts,
    is_contig,
    is_exon,
    is_gene,
    is_protein,
    is_transcript,
    normalize_feature,
    peptide_sequence,
    protein_ids,
    protein_to_cdna,
    protein_to_dna,
    protein_to_exon,
    protein_to_protein,
    protein_to_rna,
    rna_sequence,
    rna_to_cdna,
    rna_to_dna,
    rna_to_exon,
    rna_to_protein,
    rna_to_rna,
    transcript_ids,
    transcript_names,
)


# -------------------------------------------------------------------------------------------------
# test all_<feature_symbol>s
# -------------------------------------------------------------------------------------------------
def test_all_contig_ids():
    result = all_contig_ids(release=100, species="homo_sapiens")
    assert isinstance(result, list)
    assert "16" in result


def test_all_exon_ids():
    result = all_exon_ids(release=100, species="homo_sapiens")
    assert isinstance(result, list)
    assert "ENSE00003826864" in result


def test_all_gene_ids():
    result = all_gene_ids(release=100, species="homo_sapiens")
    assert isinstance(result, list)
    assert "ENSG00000149925" in result


def test_all_gene_names():
    result = all_gene_names(release=100, species="homo_sapiens")
    assert isinstance(result, list)
    assert "ALDOA" in result


def test_all_protein_ids():
    result = all_protein_ids(release=100, species="homo_sapiens")
    assert isinstance(result, list)
    assert "ENSP00000494188" in result


def test_all_transcript_ids():
    result = all_transcript_ids(release=100, species="homo_sapiens")
    assert isinstance(result, list)
    assert "ENST00000643777" in result


def test_all_transcript_names():
    result = all_transcript_names(release=100, species="homo_sapiens")
    assert isinstance(result, list)
    assert "ALDOA-219" in result


# -------------------------------------------------------------------------------------------------
# test get_<feature_symbol>s
# -------------------------------------------------------------------------------------------------
def test_contig_ids_from_contig_id():
    result = contig_ids("16", release=100, species="homo_sapiens")
    assert isinstance(result, list)
    assert "16" in result


def test_contig_ids_from_contig_with_chr():
    ret_by_id = contig_ids("chr4", release=100, species="homo_sapiens")
    assert "4" in ret_by_id


def test_contig_ids_from_exon_id():
    result = contig_ids("ENSE00003826864", release=100, species="homo_sapiens")
    assert isinstance(result, list)
    assert "16" in result


def test_contig_ids_from_gene():
    ret_by_id = contig_ids("ENSG00000149925", release=100, species="homo_sapiens")
    ret_by_name = contig_ids("ALDOA")
    assert "16" in ret_by_id
    assert "16" in ret_by_name
    assert ret_by_id == ret_by_name


def test_contig_ids_from_protein_id():
    result = contig_ids("ENSP00000494188", release=100, species="homo_sapiens")
    assert isinstance(result, list)
    assert "16" in result


def test_contig_ids_from_transcript():
    ret_by_id = contig_ids("ENST00000643777", release=100, species="homo_sapiens")
    ret_by_name = contig_ids("ALDOA-219")
    assert "16" in ret_by_id
    assert "16" in ret_by_name
    assert ret_by_id == ret_by_name


def test_exon_ids_from_contig_id():
    result = exon_ids("16", release=100, species="homo_sapiens")
    assert isinstance(result, list)
    assert "ENSE00003826864" in result


def test_exon_ids_from_exon_id():
    result = exon_ids("ENSE00003826864", release=100, species="homo_sapiens")
    assert isinstance(result, list)
    assert "ENSE00003826864" in result


def test_exon_ids_from_gene():
    ret_by_id = exon_ids("ENSG00000149925", release=100, species="homo_sapiens")
    ret_by_name = exon_ids("ALDOA", release=100, species="homo_sapiens")
    assert "ENSE00003826864" in ret_by_id
    assert "ENSE00003826864" in ret_by_name
    assert ret_by_id == ret_by_name


def test_exon_ids_from_protein_id():
    result = exon_ids("ENSP00000494188", release=100, species="homo_sapiens")
    assert isinstance(result, list)
    assert "ENSE00003826864" in result


def test_exon_ids_from_transcript():
    ret_by_id = exon_ids("ENST00000643777", release=100, species="homo_sapiens")
    ret_by_name = exon_ids("ALDOA-219", release=100, species="homo_sapiens")
    assert "ENSE00003826864" in ret_by_id
    assert "ENSE00003826864" in ret_by_name
    assert ret_by_id == ret_by_name


def test_gene_ids_different_gene_ids():
    # Gene ID changed bewteen 69 and 100
    # ENSG00000166748 AGBL1
    # ENSG00000273540 AGBL1
    assert gene_ids("AGBL1", release=69, species="homo_sapiens") == ["ENSG00000166748"]
    assert gene_ids("AGBL1", release=100, species="homo_sapiens") == ["ENSG00000273540"]


def test_gene_ids_different_gene_names():
    # Gene name changed bewteen 69 and 100
    # ENSG00000118058 MLL
    # ENSG00000118058 KMT2A
    assert gene_ids("MLL", release=69, species="homo_sapiens") == ["ENSG00000118058"]
    assert gene_ids("KMT2A", release=100, species="homo_sapiens") == ["ENSG00000118058"]


def test_gene_ids_from_contig_id():
    result = gene_ids("16", release=100, species="homo_sapiens")
    assert isinstance(result, list)
    assert "ENSG00000149925" in result


def test_gene_ids_from_exon_id():
    result = gene_ids("ENSE00003826864", release=100, species="homo_sapiens")
    assert isinstance(result, list)
    assert "ENSG00000149925" in result


def test_gene_ids_from_gene():
    ret_by_id = gene_ids("ENSG00000149925", release=100, species="homo_sapiens")
    ret_by_name = gene_ids("ALDOA", release=100, species="homo_sapiens")
    assert "ENSG00000149925" in ret_by_id
    assert "ENSG00000149925" in ret_by_name
    assert ret_by_id == ret_by_name


def test_gene_ids_from_protein_id():
    result = gene_ids("ENSP00000494188", release=100, species="homo_sapiens")
    assert isinstance(result, list)
    assert "ENSG00000149925" in result


def test_gene_ids_from_transcript():
    ret_by_id = gene_ids("ENST00000643777", release=100, species="homo_sapiens")
    ret_by_name = gene_ids("ALDOA-219", release=100, species="homo_sapiens")
    assert "ENSG00000149925" in ret_by_id
    assert "ENSG00000149925" in ret_by_name
    assert ret_by_id == ret_by_name


def test_gene_names_different_gene_ids():
    # Gene ID changed bewteen 69 and 100
    # ENSG00000204645 SSX4
    # ENSG00000268009 SSX4
    assert gene_names("ENSG00000204645", release=69, species="homo_sapiens") == ["SSX4"]
    assert gene_names("ENSG00000268009", release=100, species="homo_sapiens") == ["SSX4"]


def test_gene_names_different_gene_names():
    # Gene name changed bewteen 69 and 100
    # ENSG00000118058 MLL
    # ENSG00000118058 KMT2A
    assert gene_names("ENSG00000118058", release=69, species="homo_sapiens") == ["MLL"]
    assert gene_names("ENSG00000118058", release=100, species="homo_sapiens") == ["KMT2A"]


def test_gene_names_from_contig_id():
    result = gene_names("16", release=100, species="homo_sapiens")
    assert isinstance(result, list)
    assert "ALDOA" in result


def test_gene_names_from_exon_id():
    result = gene_names("ENSE00003826864", release=100, species="homo_sapiens")
    assert isinstance(result, list)
    assert "ALDOA" in result


def test_gene_names_from_gene():
    ret_by_id = gene_names("ENSG00000149925", release=100, species="homo_sapiens")
    ret_by_name = gene_names("ALDOA", release=100, species="homo_sapiens")
    assert "ALDOA" in ret_by_id
    assert "ALDOA" in ret_by_name
    assert ret_by_id == ret_by_name


def test_gene_names_from_protein_id():
    result = gene_names("ENSP00000494188", release=100, species="homo_sapiens")
    assert isinstance(result, list)
    assert "ALDOA" in result


def test_gene_names_from_refseq_transcript():
    ret_by_id = gene_names("NM_000244", release=100, species="homo_sapiens")
    assert "MEN1" in ret_by_id
    ret_by_id = gene_names("NM_007194", release=100, species="homo_sapiens")
    assert "CHEK2" in ret_by_id
    ret_by_id = gene_names("NM_001128849", release=100, species="homo_sapiens")
    assert "SMARCA4" in ret_by_id
    ret_by_id = gene_names("NM_000314", release=100, species="homo_sapiens")
    assert "PTEN" in ret_by_id


def test_gene_names_from_transcript():
    ret_by_id = gene_names("ENST00000643777", release=100, species="homo_sapiens")
    ret_by_name = gene_names("ALDOA-219", release=100, species="homo_sapiens")
    assert "ALDOA" in ret_by_id
    assert "ALDOA" in ret_by_name
    assert ret_by_id == ret_by_name


def test_protein_ids_from_contig_id():
    result = protein_ids("16", release=100, species="homo_sapiens")
    assert isinstance(result, list)
    assert "ENSP00000494188" in result


def test_protein_ids_from_exon_id():
    result = protein_ids("ENSE00003826864", release=100, species="homo_sapiens")
    assert isinstance(result, list)
    assert "ENSP00000494188" in result


def test_protein_ids_from_gene():
    ret_by_id = protein_ids("ENSG00000149925", release=100, species="homo_sapiens")
    ret_by_name = protein_ids("ALDOA", release=100, species="homo_sapiens")
    assert "ENSP00000494188" in ret_by_id
    assert "ENSP00000494188" in ret_by_name
    assert ret_by_id == ret_by_name


def test_protein_ids_from_protein_id():
    result = protein_ids("ENSP00000494188", release=100, species="homo_sapiens")
    assert isinstance(result, list)
    assert "ENSP00000494188" in result


def test_protein_ids_from_transcript():
    ret_by_id = protein_ids("ENST00000643777", release=100, species="homo_sapiens")
    ret_by_name = protein_ids("ALDOA-219", release=100, species="homo_sapiens")
    assert "ENSP00000494188" in ret_by_id
    assert "ENSP00000494188" in ret_by_name
    assert ret_by_id == ret_by_name


def test_transcript_ids_from_contig_id():
    result = transcript_ids("16", release=100, species="homo_sapiens")
    assert isinstance(result, list)
    assert "ENST00000643777" in result


def test_transcript_ids_from_exon_id():
    result = transcript_ids("ENSE00003826864", release=100, species="homo_sapiens")
    assert isinstance(result, list)
    assert "ENST00000643777" in result


def test_transcript_ids_from_gene():
    ret_by_id = transcript_ids("ENSG00000149925", release=100, species="homo_sapiens")
    ret_by_name = transcript_ids("ALDOA", release=100, species="homo_sapiens")
    assert "ENST00000643777" in ret_by_id
    assert "ENST00000643777" in ret_by_name
    assert ret_by_id == ret_by_name


def test_transcript_ids_from_protein_id():
    result = transcript_ids("ENSP00000494188", release=100, species="homo_sapiens")
    assert isinstance(result, list)
    assert "ENST00000643777" in result


def test_transcript_ids_from_transcript():
    ret_by_id = transcript_ids("ENST00000643777", release=100, species="homo_sapiens")
    ret_by_name = transcript_ids("ALDOA-219", release=100, species="homo_sapiens")
    assert "ENST00000643777" in ret_by_id
    assert "ENST00000643777" in ret_by_name
    assert ret_by_id == ret_by_name


def test_transcript_names_from_contig_id():
    result = transcript_names("16", release=100, species="homo_sapiens")
    assert isinstance(result, list)
    assert "ALDOA-219" in result


def test_transcript_names_from_exon_id():
    result = transcript_names("ENSE00003826864", release=100, species="homo_sapiens")
    assert isinstance(result, list)
    assert "ALDOA-219" in result


def test_transcript_names_from_gene():
    ret_by_id = transcript_names("ENSG00000149925", release=100, species="homo_sapiens")
    ret_by_name = transcript_names("ALDOA", release=100, species="homo_sapiens")
    assert "ALDOA-219" in ret_by_id
    assert "ALDOA-219" in ret_by_name
    assert ret_by_id == ret_by_name


def test_transcript_names_from_protein_id():
    result = transcript_names("ENSP00000494188", release=100, species="homo_sapiens")
    assert isinstance(result, list)
    assert "ALDOA-219" in result


def test_transcript_names_from_transcript():
    ret_by_id = transcript_names("ENST00000643777", release=100, species="homo_sapiens")
    ret_by_name = transcript_names("ALDOA-219", release=100, species="homo_sapiens")
    assert "ALDOA-219" in ret_by_id
    assert "ALDOA-219" in ret_by_name
    assert ret_by_id == ret_by_name


# -------------------------------------------------------------------------------------------------
# test is_<feature>
# -------------------------------------------------------------------------------------------------
def test_is_contig_false():
    assert is_contig("spam", release=100, species="homo_sapiens") is False


def test_is_contig_true():
    assert is_contig("4", release=100, species="homo_sapiens") is True


def test_is_exon_false():
    assert is_exon("spam", release=100, species="homo_sapiens") is False


def test_is_exon_true():
    assert is_exon("ENSE00003826864", release=100, species="homo_sapiens") is True


def test_is_gene_false():
    assert is_gene("spam", release=100, species="homo_sapiens") is False


def test_is_gene_true_from_id():
    assert is_gene("ENSG00000149925", release=100, species="homo_sapiens") is True


def test_is_gene_true_from_name():
    assert is_gene("ALDOA", release=100, species="homo_sapiens") is True


def test_is_protein_false():
    assert is_protein("spam", release=100, species="homo_sapiens") is False


def test_is_protein_true():
    assert is_protein("ENSP00000494188", release=100, species="homo_sapiens") is True


def test_is_transcript_false():
    assert is_transcript("spam", release=100, species="homo_sapiens") is False


def test_is_transcript_true_from_id():
    assert is_transcript("ENST00000643777", release=100, species="homo_sapiens") is True


def test_is_transcript_true_from_name():
    assert is_transcript("ALDOA-219", release=100, species="homo_sapiens") is True


# -------------------------------------------------------------------------------------------------
# test <feature>_sequence
# -------------------------------------------------------------------------------------------------
def test_cds_sequence():
    assert cds_sequence("ENST00000288135", 7, 9, release=100, species="homo_sapiens") == "GGC"


def test_dna_sequence():
    assert dna_sequence("4", 54695512, 54695514, "+", release=100, species="homo_sapiens") == "GCT"
    assert dna_sequence("4", 54695512, 54695514, "-", release=100, species="homo_sapiens") == "AGC"


def test_peptide_sequence():
    assert peptide_sequence("ENSP00000288135", 3, 4, release=100, species="homo_sapiens") == "GA"


def test_rna_sequence():
    assert rna_sequence("ENST00000288135", 126, 128, release=100, species="homo_sapiens") == "GCT"


# -------------------------------------------------------------------------------------------------
# test normalize_feature(feature, feature_type)
# -------------------------------------------------------------------------------------------------
def test_normalize_feature_contig_id():
    assert normalize_feature("4", release=100, species="homo_sapiens") == [("4", CONTIG_ID)]


def test_normalize_feature_exon_id():
    assert normalize_feature("ENSE00003826864", release=100, species="homo_sapiens") == [
        ("ENSE00003826864", EXON_ID)
    ]


def test_normalize_feature_gene_id():
    assert normalize_feature("ENSG00000149925", release=100, species="homo_sapiens") == [
        ("ENSG00000149925", GENE_ID)
    ]


def test_normalize_feature_gene_name():
    assert normalize_feature("ALDOA", release=100, species="homo_sapiens") == [("ALDOA", GENE_NAME)]


def test_normalize_feature_protein_id():
    assert normalize_feature("ENSP00000494188", release=100, species="homo_sapiens") == [
        ("ENSP00000494188", PROTEIN_ID)
    ]


def test_normalize_feature_refseq_transcript_id():
    assert normalize_feature("NM_000314.4", release=100, species="homo_sapiens") == [
        ("ENST00000371953", TRANSCRIPT_ID)
    ]


def test_normalize_feature_transcript_id():
    assert normalize_feature("ENST00000643777", release=100, species="homo_sapiens") == [
        ("ENST00000643777", TRANSCRIPT_ID)
    ]


def test_normalize_feature_transcript_name():
    assert normalize_feature("ALDOA-219", release=100, species="homo_sapiens") == [
        ("ALDOA-219", TRANSCRIPT_NAME)
    ]


# -------------------------------------------------------------------------------------------------
# test get_<feature>
# -------------------------------------------------------------------------------------------------
def test_get_cdna(ensembl100):
    assert get_cdna("ENST00000310581", release=100, species="homo_sapiens") == [
        CdnaMappablePosition(
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
    ]


def test_get_dna(ensembl100):
    assert get_dna("TERT", release=100, species="homo_sapiens") == [
        DnaMappablePosition(
            _data=ensembl100,
            contig_id="5",
            start=1,
            start_offset=0,
            end=181538259,
            end_offset=0,
            strand="-",
        )
    ]


def test_get_exon(ensembl100):
    assert get_exons("ENSE00003896691", release=100, species="homo_sapiens") == [
        ExonMappablePosition(
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
    ]


def test_get_gene(ensembl100):
    assert get_genes("ENSG00000164362", release=100, species="homo_sapiens") == [
        DnaMappablePosition(
            _data=ensembl100,
            contig_id="5",
            start=1253147,
            start_offset=0,
            end=1295068,
            end_offset=0,
            strand="-",
        )
    ]


def test_get_transcript(ensembl100):
    assert get_transcripts("ENST00000310581", release=100, species="homo_sapiens") == [
        RnaMappablePosition(
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
    ]


# -------------------------------------------------------------------------------------------------
# test <feature>_to_<feature>
# -------------------------------------------------------------------------------------------------
def test_cdna_to_cdna(ensembl100):
    position = CdnaMappablePosition(
        _data=ensembl100,
        contig_id="5",
        start=124,
        start_offset=0,
        end=126,
        end_offset=0,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
        protein_id="ENSP00000309572",
    )
    assert cdna_to_cdna("ENST00000310581", 124, end=126) == [position]
