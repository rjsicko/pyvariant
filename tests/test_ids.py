from coordinate_mapper.constants import CDS, CONTIG, EXON, GENE, PROTEIN, TRANSCRIPT
from coordinate_mapper.ids import (
    contig_ids,
    exon_ids,
    gene_ids,
    gene_names,
    get_exons,
    get_genes,
    get_transcripts,
    is_contig,
    is_exon,
    is_gene,
    is_protein,
    is_transcript,
    normalize_feature,
    protein_ids,
    transcript_ids,
    transcript_names,
)


def test_contig_ids_from_cds():
    ret_by_id = contig_ids("ENST00000288135", CDS)
    ret_by_name = contig_ids("KIT-201", TRANSCRIPT)
    assert "4" in ret_by_id
    assert "4" in ret_by_name
    assert ret_by_id, ret_by_name


def test_contig_ids_from_contig():
    ret_by_id = contig_ids("12", CONTIG)
    assert "12" in ret_by_id


def test_contig_ids_from_exon():
    ret_by_id = contig_ids("ENSE00001074448", EXON)
    assert "4" in ret_by_id


def test_contig_ids_from_gene():
    ret_by_id = contig_ids("ENSG00000157404", GENE)
    ret_by_name = contig_ids("KIT", GENE)
    assert "4" in ret_by_id
    assert "4" in ret_by_name
    assert ret_by_id, ret_by_name


def test_contig_ids_from_protein():
    ret_by_id = contig_ids("ENSP00000288135", PROTEIN)
    assert "4" in ret_by_id


def test_contig_ids_from_transcript():
    ret_by_id = contig_ids("ENST00000288135", TRANSCRIPT)
    ret_by_name = contig_ids("KIT-201", TRANSCRIPT)
    assert "4" in ret_by_id
    assert "4" in ret_by_name
    assert ret_by_id, ret_by_name


def test_contig_ids_from_contig_with_chr():
    ret_by_id = contig_ids("chr12", CONTIG)
    assert "12" in ret_by_id


def test_exon_ids_from_cds():
    ret_by_id = exon_ids("ENST00000288135", CDS)
    ret_by_name = exon_ids("KIT-201", CDS)
    assert "ENSE00001074448" in ret_by_id
    assert "ENSE00001074448" in ret_by_name
    assert ret_by_id, ret_by_name


def test_exon_ids_from_exon():
    ret_by_id = exon_ids("ENSE00001074448", EXON)
    assert "ENSE00001074448" in ret_by_id


def test_exon_ids_from_gene():
    ret_by_id = exon_ids("ENSG00000157404", GENE)
    ret_by_name = exon_ids("KIT", GENE)
    assert "ENSE00001074448" in ret_by_id
    assert "ENSE00001074448" in ret_by_name
    assert ret_by_id, ret_by_name


def test_exon_ids_from_protein():
    ret_by_id = exon_ids("ENSP00000288135", PROTEIN)
    assert "ENSE00001074448" in ret_by_id


def test_exon_ids_from_transcript():
    ret_by_id = exon_ids("ENST00000288135", TRANSCRIPT)
    ret_by_name = exon_ids("KIT-201", TRANSCRIPT)
    assert "ENSE00001074448" in ret_by_id
    assert "ENSE00001074448" in ret_by_name
    assert ret_by_id, ret_by_name


def test_exon_from_cds():
    ret_by_id = get_exons("ENST00000288135", CDS)
    ret_by_name = get_exons("KIT-201", CDS)
    assert "ENSE00001074448" in [i.exon_id for i in ret_by_id]
    assert "ENSE00001074448" in [i.exon_id for i in ret_by_name]
    assert ret_by_id, ret_by_name


def test_exon_from_exon():
    ret_by_id = get_exons("ENSE00001074448", EXON)
    assert "ENSE00001074448" in [i.exon_id for i in ret_by_id]


def test_exon_from_gene():
    ret_by_id = get_exons("ENSG00000157404", GENE)
    ret_by_name = get_exons("KIT", GENE)
    assert "ENSE00001074448" in [i.exon_id for i in ret_by_id]
    assert "ENSE00001074448" in [i.exon_id for i in ret_by_name]
    assert ret_by_id, ret_by_name


def test_exon_from_protein():
    ret_by_id = get_exons("ENSP00000288135", PROTEIN)
    assert "ENSE00001074448" in [i.exon_id for i in ret_by_id]


def test_exon_from_transcript():
    ret_by_id = get_exons("ENST00000288135", TRANSCRIPT)
    ret_by_name = get_exons("KIT-201", TRANSCRIPT)
    assert "ENSE00001074448" in [i.exon_id for i in ret_by_id]
    assert "ENSE00001074448" in [i.exon_id for i in ret_by_name]
    assert ret_by_id, ret_by_name


def test_gene_ids_from_cds():
    ret_by_id = gene_ids("ENST00000288135", CDS)
    ret_by_name = gene_ids("KIT-201", CDS)
    assert "ENSG00000157404" in ret_by_id
    assert "ENSG00000157404" in ret_by_name
    assert ret_by_id, ret_by_name


def test_gene_ids_from_exon():
    ret_by_id = gene_ids("ENSE00001074448", EXON)
    assert "ENSG00000157404" in ret_by_id


def test_gene_ids_from_gene():
    ret_by_id = gene_ids("ENSG00000157404", GENE)
    ret_by_name = gene_ids("KIT", GENE)
    assert "ENSG00000157404" in ret_by_id
    assert "ENSG00000157404" in ret_by_name
    assert ret_by_id, ret_by_name


def test_gene_ids_from_protein():
    ret_by_id = gene_ids("ENSP00000288135", PROTEIN)
    assert "ENSG00000157404" in ret_by_id


def test_gene_ids_from_transcript():
    ret_by_id = gene_ids("ENST00000288135", TRANSCRIPT)
    ret_by_name = gene_ids("KIT-201", TRANSCRIPT)
    assert "ENSG00000157404" in ret_by_id
    assert "ENSG00000157404" in ret_by_name
    assert ret_by_id, ret_by_name


def test_gene_names_from_cds():
    ret_by_id = gene_names("ENST00000288135", CDS)
    ret_by_name = gene_names("KIT-201", CDS)
    assert "KIT" in ret_by_id
    assert "KIT" in ret_by_name
    assert ret_by_id, ret_by_name


def test_gene_names_from_exon():
    ret_by_id = gene_names("ENSE00001074448", EXON)
    assert "KIT" in ret_by_id


def test_gene_names_from_gene():
    ret_by_id = gene_names("ENSG00000157404", GENE)
    ret_by_name = gene_names("KIT", GENE)
    assert "KIT" in ret_by_id
    assert "KIT" in ret_by_name
    assert ret_by_id, ret_by_name


def test_gene_names_from_protein():
    ret_by_id = gene_names("ENSP00000288135", PROTEIN)
    assert "KIT" in ret_by_id


def test_gene_names_from_transcript():
    ret_by_id = gene_names("ENST00000288135", TRANSCRIPT)
    ret_by_name = gene_names("KIT-201", TRANSCRIPT)
    assert "KIT" in ret_by_id
    assert "KIT" in ret_by_name
    assert ret_by_id, ret_by_name


def test_gene_names_from_refseq_transcript():
    # from SDEV-2909:
    ret_by_id = gene_names("NM_000244.3", TRANSCRIPT)
    assert "MEN1" in ret_by_id
    ret_by_id = gene_names("NM_007194.3", TRANSCRIPT)
    assert "CHEK2" in ret_by_id
    ret_by_id = gene_names("NM_001128849.1", TRANSCRIPT)
    assert "SMARCA4" in ret_by_id
    ret_by_id = gene_names("NM_000314.4", TRANSCRIPT)
    assert "PTEN" in ret_by_id


def test_gene_from_cds():
    ret_by_id = get_genes("ENST00000288135", CDS)
    ret_by_name = get_genes("KIT-201", CDS)
    assert "ENSG00000157404" in [i.gene_id for i in ret_by_id]
    assert "ENSG00000157404" in [i.gene_id for i in ret_by_name]
    assert ret_by_id, ret_by_name


def test_gene_from_exon():
    ret_by_id = get_genes("ENSE00001074448", EXON)
    assert "ENSG00000157404" in [i.gene_id for i in ret_by_id]


def test_gene_from_gene():
    ret_by_id = get_genes("ENSG00000157404", GENE)
    ret_by_name = get_genes("KIT", GENE)
    assert "ENSG00000157404" in [i.gene_id for i in ret_by_id]
    assert "ENSG00000157404" in [i.gene_id for i in ret_by_name]
    assert ret_by_id, ret_by_name


def test_gene_from_gene_old_name():
    # FAM175A was renamed to ABRAXAS1 in HG38
    ret_by_id = get_genes("ENSG00000163322", GENE)
    ret_by_name = get_genes("FAM175A", GENE)
    assert "ABRAXAS1" in [i.gene_name for i in ret_by_id]
    assert "ABRAXAS1" in [i.gene_name for i in ret_by_name]
    assert ret_by_id, ret_by_name


def test_gene_from_protein():
    ret_by_id = get_genes("ENSP00000288135", PROTEIN)
    assert "ENSG00000157404" in [i.gene_id for i in ret_by_id]


def test_gene_from_transcript():
    ret_by_id = get_genes("ENST00000288135", TRANSCRIPT)
    ret_by_name = get_genes("KIT-201", TRANSCRIPT)
    assert "ENSG00000157404" in [i.gene_id for i in ret_by_id]
    assert "ENSG00000157404" in [i.gene_id for i in ret_by_name]
    assert ret_by_id, ret_by_name


def test_is_contig_false():
    assert is_contig("spam") is False


def test_is_contig_true():
    assert is_contig("2") is True


def test_is_contig_true_alias():
    assert is_contig("chr2") is True


def test_is_exon_false():
    assert is_exon("spam") is False


def test_is_exon_true():
    assert is_exon("ENSE00003659301") is True


def test_is_gene_false():
    assert is_gene("spam") is False


def test_is_gene_true_from_id():
    assert is_gene("ENSG00000157404") is True


def test_is_gene_true_from_name():
    assert is_gene("KIT") is True


def test_is_protein_false():
    assert is_protein("spam") is False


def test_is_protein_true():
    assert is_protein("ENSP00000288135") is True


def test_is_transcript_false():
    assert is_transcript("spam") is False


def test_is_transcript_true_from_id():
    assert is_transcript("ENST00000288135") is True


def test_is_transcript_true_from_name():
    assert is_transcript("KIT-201") is True


def test_normalize_feature_1():
    assert normalize_feature("17") == [("17", CONTIG)]


def test_normalize_feature_2():
    assert normalize_feature("ENSE00003659301") == [("ENSE00003659301", EXON)]


def test_normalize_feature_3():
    assert normalize_feature("ENSG00000157404") == [("ENSG00000157404", GENE)]


def test_normalize_feature_4():
    assert normalize_feature("KIT") == [("KIT", GENE)]


def test_normalize_feature_5():
    assert normalize_feature("ENSP00000288135") == [("ENSP00000288135", PROTEIN)]


def test_normalize_feature_6():
    assert normalize_feature("ENST00000288135") == [("ENST00000288135", TRANSCRIPT)]


def test_normalize_feature_7():
    assert normalize_feature("KIT-201") == [("KIT-201", TRANSCRIPT)]


def test_normalize_feature_8():
    assert normalize_feature("NM_000314.4") == [
        ("ENST00000371953", TRANSCRIPT),
        ("ENST00000645317", TRANSCRIPT),
    ]


def test_protein_ids_from_cds():
    ret_by_id = protein_ids("ENST00000288135", CDS)
    ret_by_name = protein_ids("KIT-201", TRANSCRIPT)
    assert "ENSP00000288135" in ret_by_id
    assert "ENSP00000288135" in ret_by_name
    assert ret_by_id, ret_by_name


def test_protein_ids_from_exon():
    ret_by_id = protein_ids("ENSE00001074448", EXON)
    assert "ENSP00000288135" in ret_by_id


def test_protein_ids_from_gene():
    ret_by_id = protein_ids("ENSG00000157404", GENE)
    ret_by_name = protein_ids("KIT", GENE)
    assert "ENSP00000288135" in ret_by_id
    assert "ENSP00000288135" in ret_by_name
    assert ret_by_id, ret_by_name


def test_protein_ids_from_protein():
    ret_by_id = protein_ids("ENSP00000288135", PROTEIN)
    assert "ENSP00000288135" in ret_by_id


def test_protein_ids_from_transcript():
    ret_by_id = protein_ids("ENST00000288135", TRANSCRIPT)
    ret_by_name = protein_ids("KIT-201", TRANSCRIPT)
    assert "ENSP00000288135" in ret_by_id
    assert "ENSP00000288135" in ret_by_name
    assert ret_by_id, ret_by_name


def test_transcript_ids_from_cds():
    ret_by_id = transcript_ids("ENST00000288135", CDS)
    ret_by_name = transcript_ids("KIT-201", CDS)
    assert "ENST00000288135" in ret_by_id
    assert "ENST00000288135" in ret_by_name
    assert ret_by_id, ret_by_name


def test_transcript_ids_from_exon():
    ret_by_id = transcript_ids("ENSE00001074448", EXON)
    assert "ENST00000288135" in ret_by_id


def test_transcript_ids_from_gene():
    ret_by_id = transcript_ids("ENSG00000157404", GENE)
    ret_by_name = transcript_ids("KIT", GENE)
    assert "ENST00000288135" in ret_by_id
    assert "ENST00000288135" in ret_by_name
    assert ret_by_id, ret_by_name


def test_transcript_ids_from_protein():
    ret_by_id = transcript_ids("ENSP00000288135", PROTEIN)
    assert "ENST00000288135" in ret_by_id


def test_transcript_ids_from_transcript():
    ret_by_id = transcript_ids("ENST00000288135", TRANSCRIPT)
    ret_by_name = transcript_ids("KIT-201", TRANSCRIPT)
    assert "ENST00000288135" in ret_by_id
    assert "ENST00000288135" in ret_by_name
    assert ret_by_id, ret_by_name


def test_transcript_ids_from_gene_best_only():
    # assumes ENST00000288135 is the best transcript for ENSG00000157404 (KIT)
    ret_all = transcript_ids("ENSG00000157404", GENE, best_only=False)
    ret_best = transcript_ids("ENSG00000157404", GENE, best_only=True)
    assert "ENST00000412167" in ret_all
    assert "ENST00000412167" not in ret_best


def test_transcript_names_from_cds():
    ret_by_id = transcript_names("ENST00000288135", CDS)
    ret_by_name = transcript_names("KIT-201", CDS)
    assert "KIT-201" in ret_by_id
    assert "KIT-201" in ret_by_name
    assert ret_by_id, ret_by_name


def test_transcript_names_from_exon():
    ret_by_id = transcript_names("ENSE00001074448", EXON)
    assert "KIT-201" in ret_by_id


def test_transcript_names_from_gene():
    ret_by_id = transcript_names("ENSG00000157404", GENE)
    ret_by_name = transcript_names("KIT", GENE)
    assert "KIT-201" in ret_by_id
    assert "KIT-201" in ret_by_name
    assert ret_by_id, ret_by_name


def test_transcript_names_from_protein():
    ret_by_id = transcript_names("ENSP00000288135", PROTEIN)
    assert "KIT-201" in ret_by_id


def test_transcript_names_from_transcript():
    ret_by_id = transcript_names("ENST00000288135", TRANSCRIPT)
    ret_by_name = transcript_names("KIT-201", TRANSCRIPT)
    assert "KIT-201" in ret_by_id
    assert "KIT-201" in ret_by_name
    assert ret_by_id, ret_by_name


def test_transcript_from_cds():
    ret_by_id = get_transcripts("ENST00000288135", CDS)
    ret_by_name = get_transcripts("KIT-201", CDS)
    assert "ENST00000288135" in [i.transcript_id for i in ret_by_id]
    assert "ENST00000288135" in [i.transcript_id for i in ret_by_name]
    assert ret_by_id, ret_by_name


def test_transcript_from_exon():
    ret_by_id = get_transcripts("ENSE00003659301", EXON)
    assert "ENST00000380152" in [i.transcript_id for i in ret_by_id]


def test_transcript_from_gene():
    ret_by_id = get_transcripts("ENSG00000157404", GENE)
    ret_by_name = get_transcripts("KIT", GENE)
    assert "ENST00000288135" in [i.transcript_id for i in ret_by_id]
    assert "ENST00000288135" in [i.transcript_id for i in ret_by_name]
    assert ret_by_id, ret_by_name


def test_transcript_from_protein():
    ret_by_id = get_transcripts("ENSP00000288135", PROTEIN)
    assert "ENST00000288135" in [i.transcript_id for i in ret_by_id]


def test_transcript_from_transcript():
    ret_by_id = get_transcripts("ENST00000288135", TRANSCRIPT)
    ret_by_name = get_transcripts("KIT-201", TRANSCRIPT)
    assert "ENST00000288135" in [i.transcript_id for i in ret_by_id]
    assert "ENST00000288135" in [i.transcript_id for i in ret_by_name]
    assert ret_by_id, ret_by_name
