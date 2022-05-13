import pytest

from coordinate_mapper.constants import CDS, CONTIG, EXON, GENE, PROTEIN, TRANSCRIPT
from coordinate_mapper.map import (
    _transcript_ids_with_exon,
    cds_to_cds,
    cds_to_contig,
    cds_to_exon,
    cds_to_gene,
    cds_to_protein,
    cds_to_transcript,
    contig_to_cds,
    contig_to_contig,
    contig_to_exon,
    contig_to_gene,
    contig_to_protein,
    contig_to_transcript,
    exon_to_cds,
    exon_to_contig,
    exon_to_exon,
    exon_to_gene,
    exon_to_protein,
    exon_to_transcript,
    gene_to_cds,
    gene_to_contig,
    gene_to_exon,
    gene_to_gene,
    gene_to_protein,
    gene_to_transcript,
    get_map_func,
    protein_to_cds,
    protein_to_contig,
    protein_to_exon,
    protein_to_gene,
    protein_to_protein,
    protein_to_transcript,
    transcript_to_cds,
    transcript_to_contig,
    transcript_to_exon,
    transcript_to_gene,
    transcript_to_protein,
    transcript_to_transcript,
)


def test_cds_to_cds_pos_strand_by_transcript_id():
    pos = cds_to_cds("ENST00000380152", 1)
    assert isinstance(pos, list)
    assert pos[0].transcript_id, "ENST00000380152"
    assert pos[0].start == 1
    assert pos[0].end == 1
    assert pos[0].strand == "+"


def test_cds_to_cds_pos_strand_by_transcript_id_2():
    pos = cds_to_cds("ENST00000380152", 100)
    assert isinstance(pos, list)
    assert pos[0].transcript_id, "ENST00000380152"
    assert pos[0].start == 100
    assert pos[0].end == 100
    assert pos[0].strand == "+"


def test_cds_to_cds_pos_strand_by_transcript_id_3():
    pos = cds_to_cds("ENST00000380152", 100, 101)
    assert isinstance(pos, list)
    assert pos[0].transcript_id, "ENST00000380152"
    assert pos[0].start == 100
    assert pos[0].end == 101
    assert pos[0].strand == "+"


def test_cds_to_cds_pos_strand_by_transcript_name():
    pos = cds_to_cds("BRCA2-201", 1)
    assert isinstance(pos, list)
    assert pos[0].transcript_id, "ENST00000380152"
    assert pos[0].start == 1
    assert pos[0].end == 1
    assert pos[0].strand == "+"


def test_cds_to_cds_pos_strand_by_transcript_name_2():
    pos = cds_to_cds("BRCA2-201", 360)
    assert isinstance(pos, list)
    assert pos[0].transcript_id, "ENST00000380152"
    assert pos[0].start == 360
    assert pos[0].end == 360
    assert pos[0].strand == "+"


def test_cds_to_cds_pos_strand_by_transcript_name_3():
    pos = cds_to_cds("BRCA2-201", 360, 363)
    assert isinstance(pos, list)
    assert pos[0].transcript_id, "ENST00000380152"
    assert pos[0].start == 360
    assert pos[0].end == 363
    assert pos[0].strand == "+"


def test_cds_to_cds_by_contig_id():
    # NOTE: I don't see this query being useful because it will return coordinates for every transcript on chr12.
    # result = sorted([i.transcript_id for i in cds_to_cds("12", 1)])
    pass


def test_cds_to_cds_by_exon_id():
    result = sorted([i.transcript_id for i in cds_to_cds("ENSE00000936617", 1)])
    assert result == ["ENST00000256078", "ENST00000311936", "ENST00000556131", "ENST00000557334"]


def test_cds_to_cds_by_gene_id():
    result = sorted([i.transcript_id for i in cds_to_cds("ENSG00000133703", 1)])
    assert result == ["ENST00000256078", "ENST00000311936", "ENST00000556131", "ENST00000557334"]


def test_cds_to_cds_by_gene_name():
    result = sorted([i.transcript_id for i in cds_to_cds("KRAS", 1)])
    assert result == ["ENST00000256078", "ENST00000311936", "ENST00000556131", "ENST00000557334"]


def test_cds_to_cds_by_protein_id():
    result = sorted([i.transcript_id for i in cds_to_cds("ENSP00000256078", 1)])
    assert result == ["ENST00000256078"]


def test_cds_to_contig_pos_strand_by_transcript_id():
    pos = cds_to_contig("ENST00000380152", 1)
    assert isinstance(pos, list)
    assert pos[0].contig == "13"
    assert pos[0].start == 32316461
    assert pos[0].end == 32316461
    assert pos[0].strand == "+"


def test_cds_to_contig_pos_strand_by_transcript_id_2():
    pos = cds_to_contig("ENST00000380152", 100)
    assert isinstance(pos, list)
    assert pos[0].contig == "13"
    assert pos[0].start == 32319109
    assert pos[0].end == 32319109
    assert pos[0].strand == "+"


def test_cds_to_contig_pos_strand_by_transcript_id_3():
    pos = cds_to_contig("ENST00000380152", 100, 103)
    assert isinstance(pos, list)
    assert pos[0].contig == "13"
    assert pos[0].start == 32319109
    assert pos[0].end == 32319112
    assert pos[0].strand == "+"


def test_cds_to_contig_pos_strand_by_transcript_name():
    pos = cds_to_contig("BRCA2-201", 1)
    assert isinstance(pos, list)
    assert pos[0].contig == "13"
    assert pos[0].start == 32316461
    assert pos[0].end == 32316461
    assert pos[0].strand == "+"


def test_cds_to_contig_pos_strand_by_transcript_name_2():
    pos = cds_to_contig("BRCA2-201", 100)
    assert isinstance(pos, list)
    assert pos[0].contig == "13"
    assert pos[0].start == 32319109
    assert pos[0].end == 32319109
    assert pos[0].strand == "+"


def test_cds_to_contig_pos_strand_by_transcript_name_3():
    pos = cds_to_contig("BRCA2-201", 100, 103)
    assert isinstance(pos, list)
    assert pos[0].contig == "13"
    assert pos[0].start == 32319109
    assert pos[0].end == 32319112
    assert pos[0].strand == "+"


def test_cds_to_contig_neg_strand_by_transcript_id():
    pos = cds_to_contig("ENST00000256078", 1)
    assert isinstance(pos, list)
    assert pos[0].contig == "12"
    assert pos[0].start == 25245384
    assert pos[0].end == 25245384
    assert pos[0].strand == "+"


def test_cds_to_contig_neg_strand_by_transcript_id_2():
    pos = cds_to_contig("ENST00000256078", 100)
    assert isinstance(pos, list)
    assert pos[0].contig == "12"
    assert pos[0].start == 25245285
    assert pos[0].end == 25245285
    assert pos[0].strand == "+"


def test_cds_to_contig_neg_strand_by_transcript_id_3():
    pos = cds_to_contig("ENST00000256078", 100, 103)
    assert isinstance(pos, list)
    assert pos[0].contig == "12"
    assert pos[0].start == 25245282
    assert pos[0].end == 25245285
    assert pos[0].strand == "+"


def test_cds_to_contig_neg_strand_by_transcript_name():
    pos = cds_to_contig("KRAS-201", 1)
    assert isinstance(pos, list)
    assert pos[0].contig == "12"
    assert pos[0].start == 25245384
    assert pos[0].end == 25245384
    assert pos[0].strand == "+"


def test_cds_to_contig_neg_strand_by_transcript_name_2():
    pos = cds_to_contig("KRAS-201", 100)
    assert isinstance(pos, list)
    assert pos[0].contig == "12"
    assert pos[0].start == 25245285
    assert pos[0].end == 25245285
    assert pos[0].strand == "+"


def test_cds_to_contig_neg_strand_by_transcript_name_3():
    pos = cds_to_contig("KRAS-201", 100, 103)
    assert isinstance(pos, list)
    assert pos[0].contig == "12"
    assert pos[0].start == 25245282
    assert pos[0].end == 25245285
    assert pos[0].strand == "+"


def test_cds_to_contig_by_contig_id():
    result = sorted([i.contig for i in cds_to_contig("ENSE00000936617", 1)])
    assert result == ["12"]


def test_cds_to_contig_by_exon_id():
    result = sorted([i.contig for i in cds_to_contig("ENSE00000936617", 1)])
    assert result == ["12"]


def test_cds_to_contig_by_gene_id():
    result = sorted([i.contig for i in cds_to_contig("ENSG00000133703", 1)])
    assert result == ["12"]


def test_cds_to_contig_by_gene_name():
    result = sorted([i.contig for i in cds_to_contig("KRAS", 1)])
    assert result == ["12"]


def test_cds_to_contig_by_protein_id():
    result = sorted([i.contig for i in cds_to_contig("ENSP00000256078", 1)])
    assert result == ["12"]


def test_cds_to_exon_pos_strand_by_transcript_id():
    pos = cds_to_exon("ENST00000380152", 1)
    assert isinstance(pos, list)
    assert pos[0].transcript_id, "ENST00000380152"
    assert pos[0].exon_id, "ENSE00001484009"
    assert pos[0].start == 32316422
    assert pos[0].end == 32316527
    assert pos[0].strand == "+"


def test_cds_to_exon_pos_strand_by_transcript_id_2():
    pos = cds_to_exon("ENST00000380152", 360)
    assert isinstance(pos, list)
    assert pos[0].transcript_id, "ENST00000380152"
    assert pos[0].exon_id, "ENSE00003659301"
    assert pos[0].start == 32325076
    assert pos[0].end == 32325184
    assert pos[0].strand == "+"


def test_cds_to_exon_pos_strand_by_transcript_id_3():
    pos = cds_to_exon("ENST00000380152", 360, 380)
    assert isinstance(pos, list)
    assert pos[0].transcript_id, "ENST00000380152"
    assert pos[0].exon_id, "ENSE00003659301"
    assert pos[0].start == 32325076
    assert pos[0].end == 32325184
    assert pos[0].strand == "+"


def test_cds_to_exon_pos_strand_by_transcript_name():
    pos = cds_to_exon("BRCA2-201", 1)
    assert isinstance(pos, list)
    assert pos[0].transcript_id, "ENST00000380152"
    assert pos[0].exon_id, "ENSE00001484009"
    assert pos[0].start == 32316422
    assert pos[0].end == 32316527
    assert pos[0].strand == "+"


def test_cds_to_exon_pos_strand_by_transcript_name_2():
    pos = cds_to_exon("BRCA2-201", 360)
    assert isinstance(pos, list)
    assert pos[0].transcript_id, "ENST00000380152"
    assert pos[0].exon_id, "ENSE00003659301"
    assert pos[0].start == 32325076
    assert pos[0].end == 32325184
    assert pos[0].strand == "+"


def test_cds_to_exon_pos_strand_by_transcript_name_3():
    pos = cds_to_exon("BRCA2-201", 360, 380)
    assert isinstance(pos, list)
    assert pos[0].transcript_id, "ENST00000380152"
    assert pos[0].exon_id, "ENSE00003659301"
    assert pos[0].start == 32325076
    assert pos[0].end == 32325184
    assert pos[0].strand == "+"


def test_cds_to_exon_neg_strand_by_transcript_id():
    pos = cds_to_exon("ENST00000256078", 1)
    assert isinstance(pos, list)
    assert pos[0].transcript_id, "ENST00000256078"
    assert pos[0].exon_id, "ENSE00000936617"
    assert pos[0].start == 25245274
    assert pos[0].end == 25245395
    assert pos[0].strand == "-"


def test_cds_to_exon_neg_strand_by_transcript_id_2():
    pos = cds_to_exon("ENST00000256078", 144)
    assert isinstance(pos, list)
    assert pos[0].transcript_id, "ENST00000256078"
    assert pos[0].exon_id, "ENSE00001719809"
    assert pos[0].start == 25227234
    assert pos[0].end == 25227412
    assert pos[0].strand == "-"


def test_cds_to_exon_neg_strand_by_transcript_id_3():
    pos = cds_to_exon("ENST00000256078", 144, 156)
    assert isinstance(pos, list)
    assert pos[0].transcript_id, "ENST00000256078"
    assert pos[0].exon_id, "ENSE00001719809"
    assert pos[0].start == 25227234
    assert pos[0].end == 25227412
    assert pos[0].strand == "-"


def test_cds_to_exon_neg_strand_by_transcript_name():
    pos = cds_to_exon("KRAS-201", 1)
    assert isinstance(pos, list)
    assert pos[0].transcript_id, "ENST00000256078"
    assert pos[0].exon_id, "ENSE00000936617"
    assert pos[0].start == 25245274
    assert pos[0].end == 25245395
    assert pos[0].strand == "-"


def test_cds_to_exon_neg_strand_by_transcript_name_2():
    pos = cds_to_exon("KRAS-201", 144)
    assert isinstance(pos, list)
    assert pos[0].transcript_id, "ENST00000256078"
    assert pos[0].exon_id, "ENSE00001719809"
    assert pos[0].start == 25227234
    assert pos[0].end == 25227412
    assert pos[0].strand == "-"


def test_cds_to_exon_neg_strand_by_transcript_name_3():
    pos = cds_to_exon("KRAS-201", 144, 156)
    assert isinstance(pos, list)
    assert pos[0].transcript_id, "ENST00000256078"
    assert pos[0].exon_id, "ENSE00001719809"
    assert pos[0].start == 25227234
    assert pos[0].end == 25227412
    assert pos[0].strand == "-"


def test_cds_to_exon_different_exons():
    # positions are on different exons
    with pytest.raises(ValueError):
        _ = cds_to_exon("ENST00000380152", 300, 1000)


def test_cds_to_exon_by_contig_id():
    # NOTE: I don't see this query being useful because it will return coordinates for every transcript on chr12.
    # result = sorted([i.exon_id for i in cds_to_exon("12", 1)])
    pass


def test_cds_to_exon_by_exon_id():
    result = sorted([i.exon_id for i in cds_to_exon("ENSE00000936617", 1)])
    # the exon id appears 4 times because the exon belongs to 4 different transcripts
    assert result == ["ENSE00000936617", "ENSE00000936617", "ENSE00000936617", "ENSE00000936617"]


def test_cds_to_exon_by_gene_id():
    result = sorted([i.exon_id for i in cds_to_exon("ENSG00000133703", 1)])
    # the exon id appears 4 times because the exon belongs to 4 different transcripts
    assert result == ["ENSE00000936617", "ENSE00000936617", "ENSE00000936617", "ENSE00000936617"]


def test_cds_to_exon_by_gene_name():
    result = sorted([i.exon_id for i in cds_to_exon("KRAS", 1)])
    # the exon id appears 4 times because the exon belongs to 4 different transcripts
    assert result == ["ENSE00000936617", "ENSE00000936617", "ENSE00000936617", "ENSE00000936617"]


def test_cds_to_exon_by_protein_id():
    result = sorted([i.exon_id for i in cds_to_exon("ENSP00000256078", 1)])
    assert result == ["ENSE00000936617"]


def test_cds_to_gene_pos_strand_by_transcript_id():
    pos = cds_to_gene("ENST00000288135", 1)
    assert isinstance(pos, list)
    assert pos[0].gene_id, "ENSG00000157404"
    assert pos[0].start == 54658015
    assert pos[0].end == 54658015
    assert pos[0].strand == "+"


def test_cds_to_gene_pos_strand_by_transcript_id_2():
    pos = cds_to_gene("ENST00000288135", 213)
    assert isinstance(pos, list)
    assert pos[0].gene_id, "ENSG00000157404"
    assert pos[0].start == 54695657
    assert pos[0].end == 54695657
    assert pos[0].strand == "+"


def test_cds_to_gene_pos_strand_by_transcript_id_3():
    pos = cds_to_gene("ENST00000288135", 213, 221)
    assert isinstance(pos, list)
    assert pos[0].gene_id, "ENSG00000157404"
    assert pos[0].start == 54695657
    assert pos[0].end == 54695665
    assert pos[0].strand == "+"


def test_cds_to_gene_pos_strand_by_transcript_name():
    pos = cds_to_gene("KIT-201", 1)
    assert isinstance(pos, list)
    assert pos[0].gene_id, "ENSG00000157404"
    assert pos[0].start == 54658015
    assert pos[0].end == 54658015
    assert pos[0].strand == "+"


def test_cds_to_gene_pos_strand_by_transcript_name_2():
    pos = cds_to_gene("KIT-201", 213)
    assert isinstance(pos, list)
    assert pos[0].gene_id, "ENSG00000157404"
    assert pos[0].start == 54695657
    assert pos[0].end == 54695657
    assert pos[0].strand == "+"


def test_cds_to_gene_pos_strand_by_transcript_name_3():
    pos = cds_to_gene("KIT-201", 213, 221)
    assert isinstance(pos, list)
    assert pos[0].gene_id, "ENSG00000157404"
    assert pos[0].start == 54695657
    assert pos[0].end == 54695665
    assert pos[0].strand == "+"


def test_cds_to_gene_neg_strand_by_transcript_id():
    pos = cds_to_gene("ENST00000256078", 1)
    assert isinstance(pos, list)
    assert pos[0].gene_id, "ENSG00000133703"
    assert pos[0].start == 25245384
    assert pos[0].end == 25245384
    assert pos[0].strand == "-"


def test_cds_to_gene_neg_strand_by_transcript_id_2():
    pos = cds_to_gene("ENST00000256078", 201)
    assert isinstance(pos, list)
    assert pos[0].gene_id, "ENSG00000133703"
    assert pos[0].start == 25227323
    assert pos[0].end == 25227323
    assert pos[0].strand == "-"


def test_cds_to_gene_neg_strand_by_transcript_id_3():
    pos = cds_to_gene("ENST00000256078", 201, 301)
    assert isinstance(pos, list)
    assert pos[0].gene_id, "ENSG00000133703"
    assert pos[0].start == 25225763
    assert pos[0].end == 25227323
    assert pos[0].strand == "-"


def test_cds_to_gene_neg_strand_by_transcript_name():
    pos = cds_to_gene("KRAS-201", 1)
    assert isinstance(pos, list)
    assert pos[0].gene_id, "ENSG00000133703"
    assert pos[0].start == 25245384
    assert pos[0].end == 25245384
    assert pos[0].strand == "-"


def test_cds_to_gene_neg_strand_by_transcript_name_2():
    pos = cds_to_gene("KRAS-201", 201)
    assert isinstance(pos, list)
    assert pos[0].gene_id, "ENSG00000133703"
    assert pos[0].start == 25227323
    assert pos[0].end == 25227323
    assert pos[0].strand == "-"


def test_cds_to_gene_neg_strand_by_transcript_name_3():
    pos = cds_to_gene("KRAS-201", 201, 301)
    assert isinstance(pos, list)
    assert pos[0].gene_id, "ENSG00000133703"
    assert pos[0].start == 25225763
    assert pos[0].end == 25227323
    assert pos[0].strand == "-"


def test_cds_to_gene_by_contig_id():
    # NOTE: I don't see this query being useful because it will return coordinates for every transcript on chr12.
    # result = sorted([i.gene_id for i in cds_to_gene("12", 1)])
    pass


def test_cds_to_gene_by_exon_id():
    result = sorted([i.gene_id for i in cds_to_gene("ENSE00000936617", 1)])
    assert result == ["ENSG00000133703"]


def test_cds_to_gene_by_gene_id():
    result = sorted([i.gene_id for i in cds_to_gene("ENSG00000133703", 1)])
    assert result == ["ENSG00000133703"]


def test_cds_to_gene_by_gene_name():
    result = sorted([i.gene_id for i in cds_to_gene("KRAS", 1)])
    assert result == ["ENSG00000133703"]


def test_cds_to_gene_by_protein_id():
    result = sorted([i.gene_id for i in cds_to_gene("ENSP00000256078", 1)])
    assert result == ["ENSG00000133703"]


def test_cds_to_protein_pos_strand_by_transcript_id():
    pos = cds_to_protein("ENST00000288135", 1)
    assert isinstance(pos, list)
    assert pos[0].protein_id, "ENSP00000288135"
    assert pos[0].start == 1
    assert pos[0].end == 1
    assert pos[0].strand is None


def test_cds_to_protein_pos_strand_by_transcript_id_2():
    pos = cds_to_protein("ENST00000288135", 213)
    assert isinstance(pos, list)
    assert pos[0].protein_id, "ENSP00000288135"
    assert pos[0].start == 71
    assert pos[0].end == 71
    assert pos[0].strand is None


def test_cds_to_protein_pos_strand_by_transcript_id_3():
    pos = cds_to_protein("ENST00000288135", 213, 221)
    assert isinstance(pos, list)
    assert pos[0].protein_id, "ENSP00000288135"
    assert pos[0].start == 71
    assert pos[0].end == 74
    assert pos[0].strand is None


def test_cds_to_protein_pos_strand_by_transcript_name():
    pos = cds_to_protein("KIT-201", 1)
    assert isinstance(pos, list)
    assert pos[0].protein_id, "ENSP00000288135"
    assert pos[0].start == 1
    assert pos[0].end == 1
    assert pos[0].strand is None


def test_cds_to_protein_pos_strand_by_transcript_name_2():
    pos = cds_to_protein("KIT-201", 213)
    assert isinstance(pos, list)
    assert pos[0].protein_id, "ENSP00000288135"
    assert pos[0].start == 71
    assert pos[0].end == 71
    assert pos[0].strand is None


def test_cds_to_protein_pos_strand_by_transcript_name_3():
    pos = cds_to_protein("KIT-201", 213, 221)
    assert isinstance(pos, list)
    assert pos[0].protein_id, "ENSP00000288135"
    assert pos[0].start == 71
    assert pos[0].end == 74
    assert pos[0].strand is None


def test_cds_to_protein_neg_strand_by_transcript_id():
    pos = cds_to_protein("ENST00000256078", 1)
    assert isinstance(pos, list)
    assert pos[0].protein_id, "ENSP00000256078"
    assert pos[0].start == 1
    assert pos[0].end == 1
    assert pos[0].strand is None


def test_cds_to_protein_neg_strand_by_transcript_id_2():
    pos = cds_to_protein("ENST00000256078", 26)
    assert isinstance(pos, list)
    assert pos[0].protein_id, "ENSP00000256078"
    assert pos[0].start == 9
    assert pos[0].end == 9
    assert pos[0].strand is None


def test_cds_to_protein_neg_strand_by_transcript_id_3():
    pos = cds_to_protein("ENST00000256078", 26, 56)
    assert isinstance(pos, list)
    assert pos[0].protein_id, "ENSP00000256078"
    assert pos[0].start == 9
    assert pos[0].end == 19
    assert pos[0].strand is None


def test_cds_to_protein_neg_strand_by_transcript_name():
    pos = cds_to_protein("KRAS-201", 1)
    assert isinstance(pos, list)
    assert pos[0].protein_id, "ENSP00000256078"
    assert pos[0].start == 1
    assert pos[0].end == 1
    assert pos[0].strand is None


def test_cds_to_protein_neg_strand_by_transcript_name_2():
    pos = cds_to_protein("KRAS-201", 26)
    assert isinstance(pos, list)
    assert pos[0].protein_id, "ENSP00000256078"
    assert pos[0].start == 9
    assert pos[0].end == 9
    assert pos[0].strand is None


def test_cds_to_protein_neg_strand_by_transcript_name_3():
    pos = cds_to_protein("KRAS-201", 26, 56)
    assert isinstance(pos, list)
    assert pos[0].protein_id, "ENSP00000256078"
    assert pos[0].start == 9
    assert pos[0].end == 19
    assert pos[0].strand is None


def test_cds_to_protein_by_contig_id():
    # NOTE: I don't see this query being useful because it will return coordinates for every transcript on chr12.
    # result = sorted([i.protein_id for i in cds_to_protein("12", 1)])
    pass


def test_cds_to_protein_by_exon_id():
    result = sorted([i.protein_id for i in cds_to_protein("ENSE00000936617", 1)])
    assert result == ["ENSP00000256078", "ENSP00000308495", "ENSP00000451856", "ENSP00000452512"]


def test_cds_to_protein_by_gene_id():
    result = sorted([i.protein_id for i in cds_to_protein("ENSG00000133703", 1)])
    assert result == ["ENSP00000256078", "ENSP00000308495", "ENSP00000451856", "ENSP00000452512"]


def test_cds_to_protein_by_gene_name():
    result = sorted([i.protein_id for i in cds_to_protein("KRAS", 1)])
    assert result == ["ENSP00000256078", "ENSP00000308495", "ENSP00000451856", "ENSP00000452512"]


def test_cds_to_protein_by_protein_id():
    result = sorted([i.protein_id for i in cds_to_protein("ENSP00000256078", 1)])
    assert result == ["ENSP00000256078"]


def test_cds_to_transcript_pos_strand_by_transcript_id():
    pos = cds_to_transcript("ENST00000288135", 1)
    assert isinstance(pos, list)
    assert pos[0].transcript_id, "ENST00000288135"
    assert pos[0].start == 59
    assert pos[0].end == 59
    assert pos[0].strand == "+"


def test_cds_to_transcript_pos_strand_by_transcript_id_2():
    pos = cds_to_transcript("ENST00000288135", 213)
    assert isinstance(pos, list)
    assert pos[0].transcript_id, "ENST00000288135"
    assert pos[0].start == 271
    assert pos[0].end == 271
    assert pos[0].strand == "+"


def test_cds_to_transcript_pos_strand_by_transcript_id_3():
    pos = cds_to_transcript("ENST00000288135", 213, 221)
    assert isinstance(pos, list)
    assert pos[0].transcript_id, "ENST00000288135"
    assert pos[0].start == 271
    assert pos[0].end == 279
    assert pos[0].strand == "+"


def test_cds_to_transcript_pos_strand_by_transcript_name():
    pos = cds_to_transcript("KIT-201", 1)
    assert isinstance(pos, list)
    assert pos[0].transcript_id, "ENST00000288135"
    assert pos[0].start == 59
    assert pos[0].end == 59
    assert pos[0].strand == "+"


def test_cds_to_transcript_pos_strand_by_transcript_name_2():
    pos = cds_to_transcript("KIT-201", 213)
    assert isinstance(pos, list)
    assert pos[0].transcript_id, "ENST00000288135"
    assert pos[0].start == 271
    assert pos[0].end == 271
    assert pos[0].strand == "+"


def test_cds_to_transcript_pos_strand_by_transcript_name_3():
    pos = cds_to_transcript("KIT-201", 213, 221)
    assert isinstance(pos, list)
    assert pos[0].transcript_id, "ENST00000288135"
    assert pos[0].start == 271
    assert pos[0].end == 279
    assert pos[0].strand == "+"


def test_cds_to_transcript_neg_strand_by_transcript_id():
    pos = cds_to_transcript("ENST00000256078", 1)
    assert isinstance(pos, list)
    assert pos[0].transcript_id, "ENST00000256078"
    assert pos[0].start == 191
    assert pos[0].end == 191
    assert pos[0].strand == "-"


def test_cds_to_transcript_neg_strand_by_transcript_id_2():
    pos = cds_to_transcript("ENST00000256078", 87)
    assert isinstance(pos, list)
    assert pos[0].transcript_id, "ENST00000256078"
    assert pos[0].start == 277
    assert pos[0].end == 277
    assert pos[0].strand == "-"


def test_cds_to_transcript_neg_strand_by_transcript_id_3():
    pos = cds_to_transcript("ENST00000256078", 87, 92)
    assert isinstance(pos, list)
    assert pos[0].transcript_id, "ENST00000256078"
    assert pos[0].start == 277
    assert pos[0].end == 282
    assert pos[0].strand == "-"


def test_cds_to_transcript_neg_strand_by_transcript_name():
    pos = cds_to_transcript("KRAS-201", 1)
    assert isinstance(pos, list)
    assert pos[0].transcript_id, "ENST00000256078"
    assert pos[0].start == 191
    assert pos[0].end == 191
    assert pos[0].strand == "-"


def test_cds_to_transcript_neg_strand_by_transcript_name_2():
    pos = cds_to_transcript("KRAS-201", 87)
    assert isinstance(pos, list)
    assert pos[0].transcript_id, "ENST00000256078"
    assert pos[0].start == 277
    assert pos[0].end == 277
    assert pos[0].strand == "-"


def test_cds_to_transcript_neg_strand_by_transcript_name_3():
    pos = cds_to_transcript("KRAS-201", 87, 92)
    assert isinstance(pos, list)
    assert pos[0].transcript_id, "ENST00000256078"
    assert pos[0].start == 277
    assert pos[0].end == 282
    assert pos[0].strand == "-"


def test_cds_to_transcript_by_contig_id():
    # NOTE: I don't see this query being useful because it will return coordinates for every transcript on chr12.
    # result = sorted([i.transcript_id for i in cds_to_transcript("12", 1)])
    pass


def test_cds_to_transcript_by_exon_id():
    result = sorted([i.transcript_id for i in cds_to_transcript("ENSE00000936617", 1)])
    assert result == ["ENST00000256078", "ENST00000311936", "ENST00000556131", "ENST00000557334"]


def test_cds_to_transcript_by_gene_id():
    result = sorted([i.transcript_id for i in cds_to_transcript("ENSG00000133703", 1)])
    assert result == ["ENST00000256078", "ENST00000311936", "ENST00000556131", "ENST00000557334"]


def test_cds_to_transcript_by_gene_name():
    result = sorted([i.transcript_id for i in cds_to_transcript("KRAS", 1)])
    assert result == ["ENST00000256078", "ENST00000311936", "ENST00000556131", "ENST00000557334"]


def test_cds_to_transcript_by_protein_id():
    result = sorted([i.transcript_id for i in cds_to_transcript("ENSP00000256078", 1)])
    assert result == ["ENST00000256078"]


def test_contig_to_cds():
    pos = contig_to_cds("5", 1294501)
    assert isinstance(pos, list)
    assert pos[0].transcript_id, "ENST00000310581"
    assert pos[0].start == 385
    assert pos[0].end == 385
    assert pos[0].strand == "-"


def test_contig_to_cds_2():
    pos = contig_to_cds("5", 1294497, 1294501)
    assert isinstance(pos, list)
    assert pos[0].transcript_id, "ENST00000310581"
    assert pos[0].start == 385
    assert pos[0].end == 389
    assert pos[0].strand == "-"


def test_contig_to_cds_by_exon_id():
    result = sorted([i.transcript_id for i in contig_to_cds("ENSE00000936617", 25245275)])
    assert result == ["ENST00000256078", "ENST00000311936", "ENST00000556131", "ENST00000557334"]


def test_contig_to_cds_by_gene_id():
    result = sorted([i.transcript_id for i in contig_to_cds("ENSG00000133703", 25245275)])
    assert result == ["ENST00000256078", "ENST00000311936", "ENST00000556131", "ENST00000557334"]


def test_contig_to_cds_by_gene_name():
    result = sorted([i.transcript_id for i in contig_to_cds("KRAS", 25245275)])
    assert result == ["ENST00000256078", "ENST00000311936", "ENST00000556131", "ENST00000557334"]


def test_contig_to_cds_by_protein_id():
    result = sorted([i.transcript_id for i in contig_to_cds("ENSP00000256078", 25245275)])
    assert result == ["ENST00000256078", "ENST00000311936", "ENST00000556131", "ENST00000557334"]


def test_contig_to_cds_by_transcript_id():
    result = sorted([i.transcript_id for i in contig_to_cds("ENST00000256078", 25245275)])
    assert result == ["ENST00000256078", "ENST00000311936", "ENST00000556131", "ENST00000557334"]


def test_contig_to_cds_by_transcript_name():
    result = sorted([i.transcript_id for i in contig_to_cds("KRAS-201", 25245275)])
    assert result == ["ENST00000256078", "ENST00000311936", "ENST00000556131", "ENST00000557334"]


def test_contig_to_contig():
    pos = contig_to_contig("5", 1253147, 1295068)
    assert isinstance(pos, list)
    assert pos[0].contig == "5"
    assert pos[0].start == 1253147
    assert pos[0].end == 1295068
    assert pos[0].strand == "+"


def test_contig_to_contig_2():
    pos = contig_to_contig("5", 1254000)
    assert isinstance(pos, list)
    assert pos[0].contig == "5"
    assert pos[0].start == 1254000
    assert pos[0].end == 1254000
    assert pos[0].strand == "+"


def test_contig_to_contig_by_exon_id():
    result = sorted([i.contig for i in contig_to_contig("ENSE00000936617", 25245275)])
    assert result == ["12"]


def test_contig_to_contig_by_gene_id():
    result = sorted([i.contig for i in contig_to_contig("ENSG00000133703", 25245275)])
    assert result == ["12"]


def test_contig_to_contig_by_gene_name():
    result = sorted([i.contig for i in contig_to_contig("KRAS", 25245275)])
    assert result == ["12"]


def test_contig_to_contig_by_protein_id():
    result = sorted([i.contig for i in contig_to_contig("ENSP00000256078", 25245275)])
    assert result == ["12"]


def test_contig_to_contig_by_transcript_id():
    result = sorted([i.contig for i in contig_to_contig("ENST00000256078", 25245275)])
    assert result == ["12"]


def test_contig_to_contig_by_transcript_name():
    result = sorted([i.contig for i in contig_to_contig("KRAS-201", 25245275)])
    assert result == ["12"]


def test_contig_to_contig_with_chr():
    pos = contig_to_contig("chr5", 1253147, 1295068, feature_type=CONTIG)
    assert isinstance(pos, list)
    assert pos[0].contig == "5"
    assert pos[0].start == 1253147
    assert pos[0].end == 1295068
    assert pos[0].strand == "+"


def test_contig_to_gene():
    pos = contig_to_gene("5", 1253147, 1295068)
    assert isinstance(pos, list)
    assert pos[0].gene_id, "ENSG00000164362"
    assert pos[0].start == 1253147
    assert pos[0].end == 1295068
    assert pos[0].strand == "-"


def test_contig_to_gene_2():
    pos = contig_to_gene("5", 1254000)
    assert isinstance(pos, list)
    assert pos[0].gene_id, "ENSG00000164362"
    assert pos[0].start == 1254000
    assert pos[0].end == 1254000
    assert pos[0].strand == "-"


def test_contig_to_gene_by_exon_id():
    result = sorted([i.gene_id for i in contig_to_gene("ENSE00000936617", 25245275)])
    assert result == ["ENSG00000133703"]


def test_contig_to_gene_by_gene_id():
    result = sorted([i.gene_id for i in contig_to_gene("ENSG00000133703", 25245275)])
    assert result == ["ENSG00000133703"]


def test_contig_to_gene_by_gene_name():
    result = sorted([i.gene_id for i in contig_to_gene("KRAS", 25245275)])
    assert result == ["ENSG00000133703"]


def test_contig_to_gene_by_protein_id():
    result = sorted([i.gene_id for i in contig_to_gene("ENSP00000256078", 25245275)])
    assert result == ["ENSG00000133703"]


def test_contig_to_gene_by_transcript_id():
    result = sorted([i.gene_id for i in contig_to_gene("ENST00000256078", 25245275)])
    assert result == ["ENSG00000133703"]


def test_contig_to_gene_by_transcript_name():
    result = sorted([i.gene_id for i in contig_to_gene("KRAS-201", 25245275)])
    assert result == ["ENSG00000133703"]


def test_contig_to_protein():
    pos = contig_to_protein("5", 1294501)
    assert isinstance(pos, list)
    assert pos[0].protein_id, "ENSP00000309572"
    assert pos[0].start == 129
    assert pos[0].end == 129
    assert pos[0].strand is None


def test_contig_to_protein_2():
    pos = contig_to_protein("5", 1294497, 1294501)
    assert isinstance(pos, list)
    assert pos[0].protein_id, "ENSP00000309572"
    assert pos[0].start == 129
    assert pos[0].end == 130
    assert pos[0].strand is None


def test_contig_to_protein_by_exon_id():
    result = sorted([i.protein_id for i in contig_to_protein("ENSE00000936617", 25245275)])
    assert result == ["ENSP00000256078", "ENSP00000308495", "ENSP00000451856", "ENSP00000452512"]


def test_contig_to_protein_by_gene_id():
    result = sorted([i.protein_id for i in contig_to_protein("ENSG00000133703", 25245275)])
    assert result == ["ENSP00000256078", "ENSP00000308495", "ENSP00000451856", "ENSP00000452512"]


def test_contig_to_protein_by_gene_name():
    result = sorted([i.protein_id for i in contig_to_protein("KRAS", 25245275)])
    assert result == ["ENSP00000256078", "ENSP00000308495", "ENSP00000451856", "ENSP00000452512"]


def test_contig_to_protein_by_protein_id():
    result = sorted([i.protein_id for i in contig_to_protein("ENSP00000256078", 25245275)])
    assert result == ["ENSP00000256078", "ENSP00000308495", "ENSP00000451856", "ENSP00000452512"]


def test_contig_to_protein_by_transcript_id():
    result = sorted([i.protein_id for i in contig_to_protein("ENST00000256078", 25245275)])
    assert result == ["ENSP00000256078", "ENSP00000308495", "ENSP00000451856", "ENSP00000452512"]


def test_contig_to_protein_by_transcript_name():
    result = sorted([i.protein_id for i in contig_to_protein("KRAS-201", 25245275)])
    assert result == ["ENSP00000256078", "ENSP00000308495", "ENSP00000451856", "ENSP00000452512"]


def test_contig_to_transcript():
    pos = contig_to_transcript("5", 1294501)
    assert isinstance(pos, list)
    assert pos[0].transcript_id, "ENST00000310581"
    assert pos[0].start == 464
    assert pos[0].end == 464
    assert pos[0].strand == "-"


def test_contig_to_transcript_2():
    pos = contig_to_transcript("5", 1294497, 1294501)
    assert isinstance(pos, list)
    assert pos[0].transcript_id, "ENST00000310581"
    assert pos[0].start == 464
    assert pos[0].end == 468
    assert pos[0].strand == "-"


def test_contig_to_transcript_by_exon_id():
    result = sorted([i.transcript_id for i in contig_to_transcript("ENSE00000936617", 25245275)])

    assert result == ["ENST00000256078", "ENST00000311936", "ENST00000556131", "ENST00000557334"]


def test_contig_to_transcript_by_gene_id():
    result = sorted([i.transcript_id for i in contig_to_transcript("ENSG00000133703", 25245275)])

    assert result == ["ENST00000256078", "ENST00000311936", "ENST00000556131", "ENST00000557334"]


def test_contig_to_transcript_by_gene_name():
    result = sorted([i.transcript_id for i in contig_to_transcript("KRAS", 25245275)])
    assert result == ["ENST00000256078", "ENST00000311936", "ENST00000556131", "ENST00000557334"]


def test_contig_to_transcript_by_protein_id():
    result = sorted([i.transcript_id for i in contig_to_transcript("ENSP00000256078", 25245275)])

    assert result == ["ENST00000256078", "ENST00000311936", "ENST00000556131", "ENST00000557334"]


def test_contig_to_transcript_by_transcript_id():
    result = sorted([i.transcript_id for i in contig_to_transcript("ENST00000256078", 25245275)])

    assert result == ["ENST00000256078", "ENST00000311936", "ENST00000556131", "ENST00000557334"]


def test_contig_to_transcript_by_transcript_name():
    result = sorted([i.transcript_id for i in contig_to_transcript("KRAS-201", 25245275)])
    assert result == ["ENST00000256078", "ENST00000311936", "ENST00000556131", "ENST00000557334"]


def test_exon_to_cds_pos_strand_by_exon_id():
    pos = exon_to_cds("ENSE00003659301")
    assert isinstance(pos, list)
    assert pos[0].transcript_id, "ENST00000380152"
    assert pos[0].start == 317
    assert pos[0].end == 425
    assert pos[0].strand == "+"


def test_exon_to_cds_neg_strand_by_exon_id():
    pos = exon_to_cds("ENSE00001719809")
    assert isinstance(pos, list)
    assert pos[0].transcript_id, "ENST00000256078"
    assert pos[0].start == 112
    assert pos[0].end == 290
    assert pos[0].strand == "-"


def test_exon_to_cds_neg_strand_by_exon_id_2():
    pos = exon_to_cds("ENSE00000936617")
    assert isinstance(pos, list)
    assert pos[0].transcript_id, "ENST00000256078"
    assert pos[0].start == 1
    assert pos[0].end == 111
    assert pos[0].strand == "-"


def test_exon_to_cds_neg_strand_by_exon_id_3():
    # 'ENSE00001189804' is not part of the CDS
    with pytest.raises(ValueError):
        _ = exon_to_cds("ENSE00001189804")


def test_exon_to_cds_by_contig_id():
    # NOTE: I don't see this query being useful because it will return coordinates for every exon on chr12.
    # result = sorted([i.transcript_id for i in exon_to_cds("12")])
    pass


def test_exon_to_cds_by_gene_id():
    result = sorted([i.transcript_id for i in exon_to_cds("ENSG00000133703")])
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


def test_exon_to_cds_by_gene_name():
    result = sorted([i.transcript_id for i in exon_to_cds("KRAS")])
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


def test_exon_to_cds_by_protein_id():
    result = sorted([i.transcript_id for i in exon_to_cds("ENSP00000256078")])
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


def test_exon_to_cds_by_transcript_id():
    result = sorted([i.transcript_id for i in exon_to_cds("ENST00000256078")])
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


def test_exon_to_cds_by_transcript_name():
    result = sorted([i.transcript_id for i in exon_to_cds("KRAS-201")])
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


def test_exon_to_contig_pos_strand_by_exon_id():
    pos = exon_to_contig("ENSE00003659301")
    assert isinstance(pos, list)
    assert pos[0].contig == "13"
    assert pos[0].start == 32325076
    assert pos[0].end == 32325184
    assert pos[0].strand == "+"


def test_exon_to_contig_neg_strand_by_exon_id():
    pos = exon_to_contig("ENSE00001189807")
    assert isinstance(pos, list)
    assert pos[0].contig == "12"
    assert pos[0].start == 25215437
    assert pos[0].end == 25215560
    assert pos[0].strand == "+"


def test_exon_to_contig_by_contig_id():
    # NOTE: I don't see this query being useful because it will return coordinates for every exon on chr12.
    # result = sorted([i.contig for i in exon_to_contig("12")])
    pass


def test_exon_to_contig_by_gene_id():
    result = sorted([i.contig for i in exon_to_contig("ENSG00000133703")])
    # multiple results because they are for different positions on chr12
    assert result == ["12", "12", "12", "12", "12", "12", "12", "12", "12", "12"]


def test_exon_to_contig_by_gene_name():
    result = sorted([i.contig for i in exon_to_contig("KRAS")])
    # multiple results because they are for different positions on chr12
    assert result == ["12", "12", "12", "12", "12", "12", "12", "12", "12", "12"]


def test_exon_to_contig_by_protein_id():
    result = sorted([i.contig for i in exon_to_contig("ENSP00000256078")])
    # multiple results because they are for different positions on chr12
    assert result == ["12", "12", "12", "12", "12", "12"]


def test_exon_to_contig_by_transcript_id():
    result = sorted([i.contig for i in exon_to_contig("ENST00000256078")])
    # multiple results because they are for different positions on chr12
    assert result == ["12", "12", "12", "12", "12", "12"]


def test_exon_to_contig_by_transcript_name():
    result = sorted([i.contig for i in exon_to_contig("KRAS-201")])
    # multiple results because they are for different positions on chr12
    assert result == ["12", "12", "12", "12", "12", "12"]


def test_exon_to_exon_pos_strand_by_exon_id():
    pos = exon_to_exon("ENSE00003659301")
    assert isinstance(pos, list)
    assert pos[0].transcript_id, "ENST00000380152"
    assert pos[0].exon_id, "ENSE00003659301"
    assert pos[0].start == 32325076
    assert pos[0].end == 32325184
    assert pos[0].strand == "+"


def test_exon_to_exon_neg_strand_by_exon_id():
    pos = exon_to_exon("ENSE00001189807")
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


def test_exon_to_gene_pos_strand_by_exon_id():
    pos = exon_to_gene("ENSE00003659301")
    assert isinstance(pos, list)
    assert pos[0].gene_id, "ENSG00000139618"
    assert pos[0].start == 32325076
    assert pos[0].end == 32325184
    assert pos[0].strand == "+"


def test_exon_to_gene_neg_strand_by_exon_id():
    pos = exon_to_gene("ENSE00001189807")
    assert isinstance(pos, list)
    assert pos[0].gene_id, "ENSG00000133703"
    assert pos[0].start == 25215437
    assert pos[0].end == 25215560
    assert pos[0].strand == "-"


def test_exon_to_gene_by_contig_id():
    # NOTE: I don't see this query being useful because it will return coordinates for every exon on chr12.
    # result = sorted([i.gene_id for i in exon_to_gene("12")])
    pass


def test_exon_to_gene_by_gene_id():
    result = sorted([i.gene_id for i in exon_to_gene("ENSG00000133703")])
    assert result == [
        "ENSG00000133703",
        "ENSG00000133703",
        "ENSG00000133703",
        "ENSG00000133703",
        "ENSG00000133703",
        "ENSG00000133703",
        "ENSG00000133703",
        "ENSG00000133703",
        "ENSG00000133703",
        "ENSG00000133703",
    ]


def test_exon_to_gene_by_gene_name():
    result = sorted([i.gene_id for i in exon_to_gene("KRAS")])
    assert result == [
        "ENSG00000133703",
        "ENSG00000133703",
        "ENSG00000133703",
        "ENSG00000133703",
        "ENSG00000133703",
        "ENSG00000133703",
        "ENSG00000133703",
        "ENSG00000133703",
        "ENSG00000133703",
        "ENSG00000133703",
    ]


def test_exon_to_gene_by_protein_id():
    result = sorted([i.gene_id for i in exon_to_gene("ENSP00000256078")])
    assert result == [
        "ENSG00000133703",
        "ENSG00000133703",
        "ENSG00000133703",
        "ENSG00000133703",
        "ENSG00000133703",
        "ENSG00000133703",
    ]


def test_exon_to_gene_by_transcript_id():
    result = sorted([i.gene_id for i in exon_to_gene("ENST00000256078")])
    assert result == [
        "ENSG00000133703",
        "ENSG00000133703",
        "ENSG00000133703",
        "ENSG00000133703",
        "ENSG00000133703",
        "ENSG00000133703",
    ]


def test_exon_to_gene_by_transcript_name():
    result = sorted([i.gene_id for i in exon_to_gene("KRAS-201")])
    assert result == [
        "ENSG00000133703",
        "ENSG00000133703",
        "ENSG00000133703",
        "ENSG00000133703",
        "ENSG00000133703",
        "ENSG00000133703",
    ]


def test_exon_to_protein_pos_strand_by_exon_id():
    pos = exon_to_protein("ENSE00003659301")
    assert isinstance(pos, list)
    assert pos[0].protein_id, "ENSP00000369497"
    assert pos[0].start == 106
    assert pos[0].end == 142
    assert pos[0].strand is None


def test_exon_to_protein_neg_strand_by_exon_id():
    pos = exon_to_protein("ENSE00001719809")
    assert isinstance(pos, list)
    assert pos[0].protein_id, "ENSP00000256078"
    assert pos[0].start == 38
    assert pos[0].end == 97
    assert pos[0].strand is None


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


def test_exon_to_transcript_pos_strand_by_exon_id():
    pos = exon_to_transcript("ENSE00003659301")
    assert isinstance(pos, list)
    assert pos[0].transcript_id, "ENST00000380152"
    assert pos[0].start == 550
    assert pos[0].end == 658
    assert pos[0].strand == "+"


def test_exon_to_transcript_neg_strand_by_exon_id():
    pos = exon_to_transcript("ENSE00001719809")
    assert isinstance(pos, list)
    assert pos[0].transcript_id, "ENST00000256078"
    assert pos[0].start == 302
    assert pos[0].end == 480
    assert pos[0].strand == "-"


def test_exon_to_transcript_by_contig_id():
    # NOTE: I don't see this query being useful because it will return coordinates for every exon on chr12.
    # result = sorted([i.transcript_id for i in exon_to_transcript("12")])
    pass


def test_exon_to_transcript_by_gene_id():
    result = sorted([i.transcript_id for i in exon_to_transcript("ENSG00000133703")])
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


def test_exon_to_transcript_by_gene_name():
    result = sorted([i.transcript_id for i in exon_to_transcript("KRAS")])
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


def test_exon_to_transcript_by_protein_id():
    result = sorted([i.transcript_id for i in exon_to_transcript("ENSP00000256078")])
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


def test_exon_to_transcript_by_transcript_id():
    result = sorted([i.transcript_id for i in exon_to_transcript("ENST00000256078")])
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


def test_exon_to_transcript_by_transcript_name():
    result = sorted([i.transcript_id for i in exon_to_transcript("KRAS-201")])
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


def test_gene_to_cds_pos_strand_by_gene_id():
    pos = gene_to_cds("ENSG00000157404", 54658015)
    assert isinstance(pos, list)
    assert pos[0].transcript_id, "ENST00000288135"
    assert pos[0].start == 1
    assert pos[0].end == 1
    assert pos[0].strand == "+"


def test_gene_to_cds_pos_strand_by_gene_id_2():
    pos = gene_to_cds("ENSG00000157404", 54695656)
    assert isinstance(pos, list)
    assert pos[0].transcript_id, "ENST00000288135"
    assert pos[0].start == 212
    assert pos[0].end == 212
    assert pos[0].strand == "+"


def test_gene_to_cds_pos_strand_by_gene_id_3():
    pos = gene_to_cds("ENSG00000157404", 54695656, 54695664)
    assert isinstance(pos, list)
    assert pos[0].transcript_id, "ENST00000288135"
    assert pos[0].start == 212
    assert pos[0].end == 220
    assert pos[0].strand == "+"


def test_gene_to_cds_pos_strand_by_gene_name():
    pos = gene_to_cds("KIT", 54658015)
    assert isinstance(pos, list)
    assert pos[0].transcript_id, "ENST00000288135"
    assert pos[0].start == 1
    assert pos[0].end == 1
    assert pos[0].strand == "+"


def test_gene_to_cds_pos_strand_by_gene_name_2():
    pos = gene_to_cds("KIT", 54695656)
    assert isinstance(pos, list)
    assert pos[0].transcript_id, "ENST00000288135"
    assert pos[0].start == 212
    assert pos[0].end == 212
    assert pos[0].strand == "+"


def test_gene_to_cds_pos_strand_by_gene_name_3():
    pos = gene_to_cds("KIT", 54695656, 54695664)
    assert isinstance(pos, list)
    assert pos[0].transcript_id, "ENST00000288135"
    assert pos[0].start == 212
    assert pos[0].end == 220
    assert pos[0].strand == "+"


def test_gene_to_cds_neg_strand_by_gene_id():
    pos = gene_to_cds("ENSG00000133703", 25245384)
    assert isinstance(pos, list)
    assert pos[0].transcript_id, "ENST00000256078"
    assert pos[0].start == 1
    assert pos[0].end == 1
    assert pos[0].strand == "-"


def test_gene_to_cds_neg_strand_by_gene_id_2():
    pos = gene_to_cds("ENSG00000133703", 25225740)
    assert isinstance(pos, list)
    assert pos[0].transcript_id, "ENST00000256078"
    assert pos[0].start == 324
    assert pos[0].end == 324
    assert pos[0].strand == "-"


def test_gene_to_cds_neg_strand_by_gene_id_3():
    pos = gene_to_cds("ENSG00000133703", 25225740, 25225750)
    assert isinstance(pos, list)
    assert pos[0].transcript_id, "ENST00000256078"
    assert pos[0].start == 314
    assert pos[0].end == 324
    assert pos[0].strand == "-"


def test_gene_to_cds_neg_strand_by_gene_name():
    pos = gene_to_cds("KRAS", 25245384)
    assert isinstance(pos, list)
    assert pos[0].transcript_id, "ENST00000256078"
    assert pos[0].start == 1
    assert pos[0].end == 1
    assert pos[0].strand == "-"


def test_gene_to_cds_neg_strand_by_gene_name_2():
    pos = gene_to_cds("KRAS", 25225740)
    assert isinstance(pos, list)
    assert pos[0].transcript_id, "ENST00000256078"
    assert pos[0].start == 324
    assert pos[0].end == 324
    assert pos[0].strand == "-"


def test_gene_to_cds_neg_strand_by_gene_name_3():
    pos = gene_to_cds("KRAS", 25225740, 25225750)
    assert isinstance(pos, list)
    assert pos[0].transcript_id, "ENST00000256078"
    assert pos[0].start == 314
    assert pos[0].end == 324
    assert pos[0].strand == "-"


def test_gene_to_cds_by_contig_id():
    # TODO: This query can be useful but right now it is very slow because it iterates through every gene on chr12.
    # result = sorted([i.transcript_id for i in gene_to_cds("12", 25245275)])
    # assert result == ["ENST00000256078", "ENST00000311936", "ENST00000556131", "ENST00000557334"]
    pass


def test_gene_to_cds_by_exon_id():
    result = sorted([i.transcript_id for i in gene_to_cds("ENSE00000936617", 25245275)])
    assert result == ["ENST00000256078", "ENST00000311936", "ENST00000556131", "ENST00000557334"]


def test_gene_to_cds_by_protein_id():
    result = sorted([i.transcript_id for i in gene_to_cds("ENSP00000256078", 25245275)])
    assert result == ["ENST00000256078", "ENST00000311936", "ENST00000556131", "ENST00000557334"]


def test_gene_to_cds_by_transcript_id():
    result = sorted([i.transcript_id for i in gene_to_cds("ENST00000256078", 25245275)])
    assert result == ["ENST00000256078", "ENST00000311936", "ENST00000556131", "ENST00000557334"]


def test_gene_to_cds_by_transcript_name():
    result = sorted([i.transcript_id for i in gene_to_cds("KRAS-201", 25245275)])
    assert result == ["ENST00000256078", "ENST00000311936", "ENST00000556131", "ENST00000557334"]


def test_gene_to_contig_pos_strand_by_gene_id():
    pos = gene_to_contig("ENSG00000157404", 54658015)
    assert isinstance(pos, list)
    assert pos[0].contig == "4"
    assert pos[0].start == 54658015
    assert pos[0].end == 54658015
    assert pos[0].strand == "+"


def test_gene_to_contig_pos_strand_by_gene_id_2():
    pos = gene_to_contig("ENSG00000157404", 54695656)
    assert isinstance(pos, list)
    assert pos[0].contig == "4"
    assert pos[0].start == 54695656
    assert pos[0].end == 54695656
    assert pos[0].strand == "+"


def test_gene_to_contig_pos_strand_by_gene_id_3():
    pos = gene_to_contig("ENSG00000157404", 54695656, 54695664)
    assert isinstance(pos, list)
    assert pos[0].contig == "4"
    assert pos[0].start == 54695656
    assert pos[0].end == 54695664
    assert pos[0].strand == "+"


def test_gene_to_contig_pos_strand_by_gene_name():
    pos = gene_to_contig("KIT", 54658015)
    assert isinstance(pos, list)
    assert pos[0].contig == "4"
    assert pos[0].start == 54658015
    assert pos[0].end == 54658015
    assert pos[0].strand == "+"


def test_gene_to_contig_pos_strand_by_gene_name_2():
    pos = gene_to_contig("KIT", 54695656)
    assert isinstance(pos, list)
    assert pos[0].contig == "4"
    assert pos[0].start == 54695656
    assert pos[0].end == 54695656
    assert pos[0].strand == "+"


def test_gene_to_contig_pos_strand_by_gene_name_3():
    pos = gene_to_contig("KIT", 54695656, 54695664)
    assert isinstance(pos, list)
    assert pos[0].contig == "4"
    assert pos[0].start == 54695656
    assert pos[0].end == 54695664
    assert pos[0].strand == "+"


def test_gene_to_contig_neg_strand_by_gene_id():
    pos = gene_to_contig("ENSG00000133703", 25245384)
    assert isinstance(pos, list)
    assert pos[0].contig == "12"
    assert pos[0].start == 25245384
    assert pos[0].end == 25245384
    assert pos[0].strand == "+"


def test_gene_to_contig_neg_strand_by_gene_id_2():
    pos = gene_to_contig("ENSG00000133703", 25225740)
    assert isinstance(pos, list)
    assert pos[0].contig == "12"
    assert pos[0].start == 25225740
    assert pos[0].end == 25225740
    assert pos[0].strand == "+"


def test_gene_to_contig_neg_strand_by_gene_id_3():
    pos = gene_to_contig("ENSG00000133703", 25225740, 25225750)
    assert isinstance(pos, list)
    assert pos[0].contig == "12"
    assert pos[0].start == 25225740
    assert pos[0].end == 25225750
    assert pos[0].strand == "+"


def test_gene_to_contig_neg_strand_by_gene_name():
    pos = gene_to_contig("KRAS", 25245384)
    assert isinstance(pos, list)
    assert pos[0].contig == "12"
    assert pos[0].start == 25245384
    assert pos[0].end == 25245384
    assert pos[0].strand == "+"


def test_gene_to_contig_neg_strand_by_gene_name_2():
    pos = gene_to_contig("KRAS", 25225740)
    assert isinstance(pos, list)
    assert pos[0].contig == "12"
    assert pos[0].start == 25225740
    assert pos[0].end == 25225740
    assert pos[0].strand == "+"


def test_gene_to_contig_neg_strand_by_gene_name_3():
    pos = gene_to_contig("KRAS", 25225740, 25225750)
    assert isinstance(pos, list)
    assert pos[0].contig == "12"
    assert pos[0].start == 25225740
    assert pos[0].end == 25225750
    assert pos[0].strand == "+"


def test_gene_to_contig_by_contig_id():
    # TODO: This query can be useful but right now it is very slow because it iterates through every gene on chr12.
    # result = sorted([i.contig for i in gene_to_contig("12", 25245275)])
    # assert result == ["12"]
    pass


def test_gene_to_contig_by_exon_id():
    result = sorted([i.contig for i in gene_to_contig("ENSE00000936617", 25245275)])
    assert result == ["12"]


def test_gene_to_contig_by_protein_id():
    result = sorted([i.contig for i in gene_to_contig("ENSP00000256078", 25245275)])
    assert result == ["12"]


def test_gene_to_contig_by_transcript_id():
    result = sorted([i.contig for i in gene_to_contig("ENST00000256078", 25245275)])
    assert result == ["12"]


def test_gene_to_contig_by_transcript_name():
    result = sorted([i.contig for i in gene_to_contig("KRAS-201", 25245275)])
    assert result == ["12"]


def test_gene_to_exon_pos_strand_by_gene_id():
    pos = gene_to_exon("ENSG00000157404", 54657918)
    assert isinstance(pos, list)
    assert pos[0].transcript_id, "ENST00000412167"
    assert pos[0].exon_id, "ENSE00001905199"
    assert pos[0].start == 54657918
    assert pos[0].end == 54658081
    assert pos[0].strand == "+"


def test_gene_to_exon_pos_strand_by_gene_id_2():
    pos = gene_to_exon("ENSG00000157404", 54698301)
    assert isinstance(pos, list)
    assert pos[0].transcript_id, "ENST00000288135"
    assert pos[0].exon_id, "ENSE00001074448"
    assert pos[0].start == 54698284
    assert pos[0].end == 54698565
    assert pos[0].strand == "+"


def test_gene_to_exon_pos_strand_by_gene_id_3():
    pos = gene_to_exon("ENSG00000157404", 54698301, 54698330)
    assert isinstance(pos, list)
    assert pos[0].transcript_id, "ENST00000288135"
    assert pos[0].exon_id, "ENSE00001074448"
    assert pos[0].start == 54698284
    assert pos[0].end == 54698565
    assert pos[0].strand == "+"


def test_gene_to_exon_pos_strand_by_gene_name():
    pos = gene_to_exon("KIT", 54657918)
    assert isinstance(pos, list)
    assert pos[0].transcript_id, "ENST00000412167"
    assert pos[0].exon_id, "ENSE00001905199"
    assert pos[0].start == 54657918
    assert pos[0].end == 54658081
    assert pos[0].strand == "+"


def test_gene_to_exon_pos_strand_by_gene_name_2():
    pos = gene_to_exon("KIT", 54698301)
    assert isinstance(pos, list)
    assert pos[0].transcript_id, "ENST00000288135"
    assert pos[0].exon_id, "ENSE00001074448"
    assert pos[0].start == 54698284
    assert pos[0].end == 54698565
    assert pos[0].strand == "+"


def test_gene_to_exon_pos_strand_by_gene_name_3():
    pos = gene_to_exon("KIT", 54698301, 54698330)
    assert isinstance(pos, list)
    assert pos[0].transcript_id, "ENST00000288135"
    assert pos[0].exon_id, "ENSE00001074448"
    assert pos[0].start == 54698284
    assert pos[0].end == 54698565
    assert pos[0].strand == "+"


def test_gene_to_exon_neg_strand_by_gene_id():
    pos = gene_to_exon("ENSG00000133703", 25250929)
    assert isinstance(pos, list)
    assert pos[0].transcript_id, "ENST00000256078"
    assert pos[0].exon_id, "ENSE00000000028"
    assert pos[0].start == 25250751
    assert pos[0].end == 25250929
    assert pos[0].strand == "-"


def test_gene_to_exon_neg_strand_by_gene_id_2():
    pos = gene_to_exon("ENSG00000133703", 25227400)
    assert isinstance(pos, list)
    assert pos[0].transcript_id, "ENST00000256078"
    assert pos[0].exon_id, "ENSE00001719809"
    assert pos[0].start == 25227234
    assert pos[0].end == 25227412
    assert pos[0].strand == "-"


def test_gene_to_exon_neg_strand_by_gene_id_3():
    pos = gene_to_exon("ENSG00000133703", 25227380, 25227400)
    assert isinstance(pos, list)
    assert pos[0].transcript_id, "ENST00000256078"
    assert pos[0].exon_id, "ENSE00001719809"
    assert pos[0].start == 25227234
    assert pos[0].end == 25227412
    assert pos[0].strand == "-"


def test_gene_to_exon_neg_strand_by_gene_name():
    pos = gene_to_exon("KRAS", 25250929)
    assert isinstance(pos, list)
    assert pos[0].transcript_id, "ENST00000256078"
    assert pos[0].exon_id, "ENSE00000000028"
    assert pos[0].start == 25250751
    assert pos[0].end == 25250929
    assert pos[0].strand == "-"


def test_gene_to_exon_neg_strand_by_gene_name_2():
    pos = gene_to_exon("KRAS", 25227400)
    assert isinstance(pos, list)
    assert pos[0].transcript_id, "ENST00000256078"
    assert pos[0].exon_id, "ENSE00001719809"
    assert pos[0].start == 25227234
    assert pos[0].end == 25227412
    assert pos[0].strand == "-"


def test_gene_to_exon_neg_strand_by_gene_name_3():
    pos = gene_to_exon("KRAS", 25227380, 25227400)
    assert isinstance(pos, list)
    assert pos[0].transcript_id, "ENST00000256078"
    assert pos[0].exon_id, "ENSE00001719809"
    assert pos[0].start == 25227234
    assert pos[0].end == 25227412
    assert pos[0].strand == "-"


def test_gene_to_exon_by_contig_id():
    # TODO: This query can be useful but right now it is very slow because it iterates through every gene on chr12.
    # result = sorted([i.exon_id for i in gene_to_exon("12", 25245275)])
    # assert result == ['ENSE00000936617', 'ENSE00000936617', 'ENSE00000936617', 'ENSE00000936617']
    pass


def test_gene_to_exon_by_exon_id():
    result = sorted([i.exon_id for i in gene_to_exon("ENSE00000936617", 25245275)])
    assert result == ["ENSE00000936617", "ENSE00000936617", "ENSE00000936617", "ENSE00000936617"]


def test_gene_to_exon_by_protein_id():
    result = sorted([i.exon_id for i in gene_to_exon("ENSP00000256078", 25245275)])
    assert result == ["ENSE00000936617", "ENSE00000936617", "ENSE00000936617", "ENSE00000936617"]


def test_gene_to_exon_by_transcript_id():
    result = sorted([i.exon_id for i in gene_to_exon("ENST00000256078", 25245275)])
    assert result == ["ENSE00000936617", "ENSE00000936617", "ENSE00000936617", "ENSE00000936617"]


def test_gene_to_exon_by_transcript_name():
    result = sorted([i.exon_id for i in gene_to_exon("KRAS-201", 25245275)])
    assert result == ["ENSE00000936617", "ENSE00000936617", "ENSE00000936617", "ENSE00000936617"]


def test_gene_to_gene_pos_strand_by_gene_id():
    pos = gene_to_gene("ENSG00000157404", 54657918)
    assert isinstance(pos, list)
    assert pos[0].gene_id, "ENSG00000157404"
    assert pos[0].start == 54657918
    assert pos[0].end == 54657918
    assert pos[0].strand == "+"


def test_gene_to_gene_pos_strand_by_gene_id_2():
    pos = gene_to_gene("ENSG00000157404", 54658000)
    assert isinstance(pos, list)
    assert pos[0].gene_id, "ENSG00000157404"
    assert pos[0].start == 54658000
    assert pos[0].end == 54658000
    assert pos[0].strand == "+"


def test_gene_to_gene_pos_strand_by_gene_id_3():
    pos = gene_to_gene("ENSG00000157404", 54658000, 54658001)
    assert isinstance(pos, list)
    assert pos[0].gene_id, "ENSG00000157404"
    assert pos[0].start == 54658000
    assert pos[0].end == 54658001
    assert pos[0].strand == "+"


def test_gene_to_gene_pos_strand_by_gene_name():
    pos = gene_to_gene("KIT", 54657918)
    assert isinstance(pos, list)
    assert pos[0].gene_id, "ENSG00000157404"
    assert pos[0].start == 54657918
    assert pos[0].end == 54657918
    assert pos[0].strand == "+"


def test_gene_to_gene_pos_strand_by_gene_name_2():
    pos = gene_to_gene("KIT", 54658000)
    assert isinstance(pos, list)
    assert pos[0].gene_id, "ENSG00000157404"
    assert pos[0].start == 54658000
    assert pos[0].end == 54658000
    assert pos[0].strand == "+"


def test_gene_to_gene_pos_strand_by_gene_name_3():
    pos = gene_to_gene("KIT", 54658000, 54658001)
    assert isinstance(pos, list)
    assert pos[0].gene_id, "ENSG00000157404"
    assert pos[0].start == 54658000
    assert pos[0].end == 54658001
    assert pos[0].strand == "+"


def test_gene_to_gene_neg_strand_by_gene_id():
    pos = gene_to_gene("ENSG00000133703", 25205246)
    assert isinstance(pos, list)
    assert pos[0].gene_id, "ENSG00000133703"
    assert pos[0].start == 25205246
    assert pos[0].end == 25205246
    assert pos[0].strand == "-"


def test_gene_to_gene_neg_strand_by_gene_id_2():
    pos = gene_to_gene("ENSG00000133703", 25205300)
    assert isinstance(pos, list)
    assert pos[0].gene_id, "ENSG00000133703"
    assert pos[0].start == 25205300
    assert pos[0].end == 25205300
    assert pos[0].strand == "-"


def test_gene_to_gene_neg_strand_by_gene_id_3():
    pos = gene_to_gene("ENSG00000133703", 25205300, 25205301)
    assert isinstance(pos, list)
    assert pos[0].gene_id, "ENSG00000133703"
    assert pos[0].start == 25205300
    assert pos[0].end == 25205301
    assert pos[0].strand == "-"


def test_gene_to_gene_neg_strand_by_gene_name():
    pos = gene_to_gene("KRAS", 25205246)
    assert isinstance(pos, list)
    assert pos[0].gene_id, "ENSG00000133703"
    assert pos[0].start == 25205246
    assert pos[0].end == 25205246
    assert pos[0].strand == "-"


def test_gene_to_gene_neg_strand_by_gene_name_2():
    pos = gene_to_gene("KRAS", 25205300)
    assert isinstance(pos, list)
    assert pos[0].gene_id, "ENSG00000133703"
    assert pos[0].start == 25205300
    assert pos[0].end == 25205300
    assert pos[0].strand == "-"


def test_gene_to_gene_neg_strand_by_gene_name_3():
    pos = gene_to_gene("KRAS", 25205300, 25205301)
    assert isinstance(pos, list)
    assert pos[0].gene_id, "ENSG00000133703"
    assert pos[0].start == 25205300
    assert pos[0].end == 25205301
    assert pos[0].strand == "-"


def test_gene_to_gene_by_contig_id():
    # TODO: This query can be useful but right now it is very slow because it iterates through every gene on chr12.
    # result = sorted([i.gene_id for i in gene_to_gene("12", 25245275)])
    # assert result == ['ENSE00000936617', 'ENSE00000936617', 'ENSE00000936617', 'ENSE00000936617']
    pass


def test_gene_to_gene_by_exon_id():
    result = sorted([i.gene_id for i in gene_to_gene("ENSE00000936617", 25245275)])
    assert result == ["ENSG00000133703"]


def test_gene_to_gene_by_protein_id():
    result = sorted([i.gene_id for i in gene_to_gene("ENSP00000256078", 25245275)])
    assert result == ["ENSG00000133703"]


def test_gene_to_gene_by_transcript_id():
    result = sorted([i.gene_id for i in gene_to_gene("ENST00000256078", 25245275)])
    assert result == ["ENSG00000133703"]


def test_gene_to_gene_by_transcript_name():
    result = sorted([i.gene_id for i in gene_to_gene("KRAS-201", 25245275)])
    assert result == ["ENSG00000133703"]


def test_gene_to_protein_pos_strand_by_gene_id():
    pos = gene_to_protein("ENSG00000157404", 54658015)
    assert isinstance(pos, list)
    assert pos[0].protein_id, "ENSP00000288135"
    assert pos[0].start == 1
    assert pos[0].end == 1
    assert pos[0].strand is None


def test_gene_to_protein_pos_strand_by_gene_id_2():
    pos = gene_to_protein("ENSG00000157404", 54695656)
    assert isinstance(pos, list)
    assert pos[0].protein_id, "ENSP00000288135"
    assert pos[0].start == 71
    assert pos[0].end == 71
    assert pos[0].strand is None


def test_gene_to_protein_pos_strand_by_gene_id_3():
    pos = gene_to_protein("ENSG00000157404", 54695656, 54695664)
    assert isinstance(pos, list)
    assert pos[0].protein_id, "ENSP00000288135"
    assert pos[0].start == 71
    assert pos[0].end == 74
    assert pos[0].strand is None


def test_gene_to_protein_pos_strand_by_gene_name():
    pos = gene_to_protein("KIT", 54658015)
    assert isinstance(pos, list)
    assert pos[0].protein_id, "ENSP00000288135"
    assert pos[0].start == 1
    assert pos[0].end == 1
    assert pos[0].strand is None


def test_gene_to_protein_pos_strand_by_gene_name_2():
    pos = gene_to_protein("KIT", 54695656)
    assert isinstance(pos, list)
    assert pos[0].protein_id, "ENSP00000288135"
    assert pos[0].start == 71
    assert pos[0].end == 71
    assert pos[0].strand is None


def test_gene_to_protein_pos_strand_by_gene_name_3():
    pos = gene_to_protein("KIT", 54695656, 54695664)
    assert isinstance(pos, list)
    assert pos[0].protein_id, "ENSP00000288135"
    assert pos[0].start == 71
    assert pos[0].end == 74
    assert pos[0].strand is None


def test_gene_to_protein_neg_strand_by_gene_id():
    pos = gene_to_protein("ENSG00000133703", 25245384)
    assert isinstance(pos, list)
    assert pos[0].protein_id, "ENSP00000256078"
    assert pos[0].start == 1
    assert pos[0].end == 1
    assert pos[0].strand is None


def test_gene_to_protein_neg_strand_by_gene_id_2():
    pos = gene_to_protein("ENSG00000133703", 25227402)
    assert isinstance(pos, list)
    assert pos[0].protein_id, "ENSP00000256078"
    assert pos[0].start == 41
    assert pos[0].end == 41
    assert pos[0].strand is None


def test_gene_to_protein_neg_strand_by_gene_id_3():
    pos = gene_to_protein("ENSG00000133703", 25227382, 25227402)
    assert isinstance(pos, list)
    assert pos[0].protein_id, "ENSP00000256078"
    assert pos[0].start == 41
    assert pos[0].end == 48
    assert pos[0].strand is None


def test_gene_to_protein_neg_strand_by_gene_name():
    pos = gene_to_protein("KRAS", 25245384)
    assert isinstance(pos, list)
    assert pos[0].protein_id, "ENSP00000256078"
    assert pos[0].start == 1
    assert pos[0].end == 1
    assert pos[0].strand is None


def test_gene_to_protein_neg_strand_by_gene_name_2():
    pos = gene_to_protein("KRAS", 25227402)
    assert isinstance(pos, list)
    assert pos[0].protein_id, "ENSP00000256078"
    assert pos[0].start == 41
    assert pos[0].end == 41
    assert pos[0].strand is None


def test_gene_to_protein_neg_strand_by_gene_name_3():
    pos = gene_to_protein("KRAS", 25227382, 25227402)
    assert isinstance(pos, list)
    assert pos[0].protein_id, "ENSP00000256078"
    assert pos[0].start == 41
    assert pos[0].end == 48
    assert pos[0].strand is None


def test_gene_to_protein_by_contig_id():
    # TODO: This query can be useful but right now it is very slow because it iterates through every gene on chr12.
    # result = sorted([i.protein_id for i in gene_to_protein("12", 25245275)])
    # assert result == ['ENSP00000256078', 'ENSP00000308495', 'ENSP00000451856', 'ENSP00000452512']
    pass


def test_gene_to_protein_by_exon_id():
    result = sorted([i.protein_id for i in gene_to_protein("ENSE00000936617", 25245275)])
    assert result == ["ENSP00000256078", "ENSP00000308495", "ENSP00000451856", "ENSP00000452512"]


def test_gene_to_protein_by_protein_id():
    result = sorted([i.protein_id for i in gene_to_protein("ENSP00000256078", 25245275)])
    assert result == ["ENSP00000256078", "ENSP00000308495", "ENSP00000451856", "ENSP00000452512"]


def test_gene_to_protein_by_transcript_id():
    result = sorted([i.protein_id for i in gene_to_protein("ENST00000256078", 25245275)])
    assert result == ["ENSP00000256078", "ENSP00000308495", "ENSP00000451856", "ENSP00000452512"]


def test_gene_to_protein_by_transcript_name():
    result = sorted([i.protein_id for i in gene_to_protein("KRAS-201", 25245275)])
    assert result == ["ENSP00000256078", "ENSP00000308495", "ENSP00000451856", "ENSP00000452512"]


def test_gene_to_transcript_pos_strand_by_gene_id():
    pos = gene_to_transcript("ENSG00000188554", 43170481)
    assert isinstance(pos, list)
    assert pos[0].transcript_id, "ENST00000341165"
    assert pos[0].start == 1
    assert pos[0].end == 1
    assert pos[0].strand == "+"


def test_gene_to_transcript_pos_strand_by_gene_id_2():
    pos = gene_to_transcript("ENSG00000188554", 43189698)
    assert isinstance(pos, list)
    assert pos[0].transcript_id, "ENST00000341165"
    assert pos[0].start == 731
    assert pos[0].end == 731
    assert pos[0].strand == "+"


def test_gene_to_transcript_pos_strand_by_gene_id_3():
    pos = gene_to_transcript("ENSG00000188554", 43189698, 43189784)
    assert isinstance(pos, list)
    assert pos[0].transcript_id, "ENST00000341165"
    assert pos[0].start == 731
    assert pos[0].end == 817
    assert pos[0].strand == "+"


def test_gene_to_transcript_pos_strand_by_gene_name():
    pos = gene_to_transcript("NBR1", 43170481)
    assert isinstance(pos, list)
    assert pos[0].transcript_id, "ENST00000341165"
    assert pos[0].start == 1
    assert pos[0].end == 1
    assert pos[0].strand == "+"


def test_gene_to_transcript_pos_strand_by_gene_name_2():
    pos = gene_to_transcript("NBR1", 43189698)
    assert isinstance(pos, list)
    assert pos[0].transcript_id, "ENST00000341165"
    assert pos[0].start == 731
    assert pos[0].end == 731
    assert pos[0].strand == "+"


def test_gene_to_transcript_pos_strand_by_gene_name_3():
    pos = gene_to_transcript("NBR1", 43189698, 43189784)
    assert isinstance(pos, list)
    assert pos[0].transcript_id, "ENST00000341165"
    assert pos[0].start == 731
    assert pos[0].end == 817
    assert pos[0].strand == "+"


def test_gene_to_transcript_neg_strand_by_gene_id():
    pos = gene_to_transcript("ENSG00000133703", 25250929)
    assert isinstance(pos, list)
    assert pos[0].transcript_id, "ENST00000256078"
    assert pos[0].start == 1
    assert pos[0].end == 1
    assert pos[0].strand == "-"


def test_gene_to_transcript_neg_strand_by_gene_id_2():
    pos = gene_to_transcript("ENSG00000133703", 25225634)
    assert isinstance(pos, list)
    assert pos[0].transcript_id, "ENST00000256078"
    assert pos[0].start == 620
    assert pos[0].end == 620
    assert pos[0].strand == "-"


def test_gene_to_transcript_neg_strand_by_gene_id_3():
    pos = gene_to_transcript("ENSG00000133703", 25225634, 25225651)
    assert isinstance(pos, list)
    assert pos[0].transcript_id, "ENST00000256078"
    assert pos[0].start == 603
    assert pos[0].end == 620
    assert pos[0].strand == "-"


def test_gene_to_transcript_neg_strand_by_gene_name():
    pos = gene_to_transcript("KRAS", 25225634, 25225651)
    assert isinstance(pos, list)
    assert pos[0].transcript_id, "ENST00000256078"
    assert pos[0].start == 603
    assert pos[0].end == 620
    assert pos[0].strand == "-"


def test_gene_to_transcript_neg_strand_by_gene_name_2():
    pos = gene_to_transcript("KRAS", 25225634, 25225651)
    assert isinstance(pos, list)
    assert pos[0].transcript_id, "ENST00000256078"
    assert pos[0].start == 603
    assert pos[0].end == 620
    assert pos[0].strand == "-"


def test_gene_to_transcript_neg_strand_by_gene_name_3():
    pos = gene_to_transcript("KRAS", 25225634, 25225651)
    assert isinstance(pos, list)
    assert pos[0].transcript_id, "ENST00000256078"
    assert pos[0].start == 603
    assert pos[0].end == 620
    assert pos[0].strand == "-"


def test_gene_to_transcript_by_contig_id():
    # TODO: This query can be useful but right now it is very slow because it iterates through every gene on chr12.
    # result = sorted([i.transcript_id for i in gene_to_transcript("12", 25245275)])
    # assert result == ['ENST00000256078', 'ENST00000311936', 'ENST00000556131', 'ENST00000557334']
    pass


def test_gene_to_transcript_by_exon_id():
    result = sorted([i.transcript_id for i in gene_to_transcript("ENSE00000936617", 25245275)])
    assert result == ["ENST00000256078", "ENST00000311936", "ENST00000556131", "ENST00000557334"]


def test_gene_to_transcript_by_protein_id():
    result = sorted([i.transcript_id for i in gene_to_transcript("ENSP00000256078", 25245275)])
    assert result == ["ENST00000256078", "ENST00000311936", "ENST00000556131", "ENST00000557334"]


def test_gene_to_transcript_by_transcript_id():
    result = sorted([i.transcript_id for i in gene_to_transcript("ENST00000256078", 25245275)])
    assert result == ["ENST00000256078", "ENST00000311936", "ENST00000556131", "ENST00000557334"]


def test_gene_to_transcript_by_transcript_name():
    result = sorted([i.transcript_id for i in gene_to_transcript("KRAS-201", 25245275)])
    assert result == ["ENST00000256078", "ENST00000311936", "ENST00000556131", "ENST00000557334"]


def test_protein_to_cds_pos_strand_by_protein_id():
    pos = protein_to_cds("ENSP00000288135", 1)
    assert isinstance(pos, list)
    assert pos[0].transcript_id, "ENST00000288135"
    assert pos[0].start == 1
    assert pos[0].end == 3
    assert pos[0].strand == "+"


def test_protein_to_cds_pos_strand_by_protein_id_2():
    pos = protein_to_cds("ENSP00000288135", 12)
    assert isinstance(pos, list)
    assert pos[0].transcript_id, "ENST00000288135"
    assert pos[0].start == 34
    assert pos[0].end == 36
    assert pos[0].strand == "+"


def test_protein_to_cds_pos_strand_by_protein_id_3():
    pos = protein_to_cds("ENSP00000288135", 12, 15)
    assert isinstance(pos, list)
    assert pos[0].transcript_id, "ENST00000288135"
    assert pos[0].start == 34
    assert pos[0].end == 45
    assert pos[0].strand == "+"


def test_protein_to_cds_neg_strand_by_protein_id():
    pos = protein_to_cds("ENSP00000308495", 1)
    assert isinstance(pos, list)
    assert pos[0].transcript_id, "ENST00000311936"
    assert pos[0].start == 1
    assert pos[0].end == 3
    assert pos[0].strand == "-"


def test_protein_to_cds_neg_strand_by_protein_id_2():
    pos = protein_to_cds("ENSP00000308495", 123)
    assert isinstance(pos, list)
    assert pos[0].transcript_id, "ENST00000311936"
    assert pos[0].start == 367
    assert pos[0].end == 369
    assert pos[0].strand == "-"


def test_protein_to_cds_neg_strand_by_protein_id_3():
    pos = protein_to_cds("ENSP00000308495", 123, 124)
    assert isinstance(pos, list)
    assert pos[0].transcript_id, "ENST00000311936"
    assert pos[0].start == 367
    assert pos[0].end == 372
    assert pos[0].strand == "-"


def test_protein_to_cds_by_contig_id():
    # NOTE: I don't see this query being useful because it will return coordinates for every protein on chr12.
    # result = sorted([i.transcript_id for i in protein_to_cds("12", 1)])
    pass


def test_protein_to_cds_by_exon_id():
    result = sorted([i.transcript_id for i in protein_to_cds("ENSE00000936617", 1)])
    assert result == ["ENST00000256078", "ENST00000311936", "ENST00000556131", "ENST00000557334"]


def test_protein_to_cds_by_gene_id():
    result = sorted([i.transcript_id for i in protein_to_cds("ENSG00000133703", 1)])
    assert result == ["ENST00000256078", "ENST00000311936", "ENST00000556131", "ENST00000557334"]


def test_protein_to_cds_by_gene_name():
    result = sorted([i.transcript_id for i in protein_to_cds("KRAS", 1)])
    assert result == ["ENST00000256078", "ENST00000311936", "ENST00000556131", "ENST00000557334"]


def test_protein_to_cds_by_transcript_id():
    result = sorted([i.transcript_id for i in protein_to_cds("ENST00000256078", 1)])
    assert result == ["ENST00000256078"]


def test_protein_to_cds_by_transcript_name():
    result = sorted([i.transcript_id for i in protein_to_cds("KRAS-201", 1)])
    assert result == ["ENST00000256078"]


def test_protein_to_contig_pos_strand_by_protein_id():
    pos = protein_to_contig("ENSP00000288135", 1)
    assert isinstance(pos, list)
    assert pos[0].contig == "4"
    assert pos[0].start == 54658015
    assert pos[0].end == 54658017
    assert pos[0].strand == "+"


def test_protein_to_contig_pos_strand_by_protein_id_2():
    pos = protein_to_contig("ENSP00000288135", 71)
    assert isinstance(pos, list)
    assert pos[0].contig == "4"
    assert pos[0].start == 54695655
    assert pos[0].end == 54695657
    assert pos[0].strand == "+"


def test_protein_to_contig_pos_strand_by_protein_id_3():
    pos = protein_to_contig("ENSP00000288135", 71, 74)
    assert isinstance(pos, list)
    assert pos[0].contig == "4"
    assert pos[0].start == 54695655
    assert pos[0].end == 54695666
    assert pos[0].strand == "+"


def test_protein_to_contig_neg_strand_by_protein_id():
    pos = protein_to_contig("ENSP00000256078", 1)
    assert isinstance(pos, list)
    assert pos[0].contig == "12"
    assert pos[0].start == 25245382
    assert pos[0].end == 25245384
    assert pos[0].strand == "+"


def test_protein_to_contig_neg_strand_by_protein_id_2():
    pos = protein_to_contig("ENSP00000256078", 6)
    assert isinstance(pos, list)
    assert pos[0].contig == "12"
    assert pos[0].start == 25245367
    assert pos[0].end == 25245369
    assert pos[0].strand == "+"


def test_protein_to_contig_neg_strand_by_protein_id_3():
    pos = protein_to_contig("ENSP00000256078", 6, 9)
    assert isinstance(pos, list)
    assert pos[0].contig == "12"
    assert pos[0].start == 25245358
    assert pos[0].end == 25245369
    assert pos[0].strand == "+"


def test_protein_to_contig_by_contig_id():
    # NOTE: I don't see this query being useful because it will return coordinates for every protein on chr12.
    # result = sorted([i.contig for i in protein_to_contig("12", 1)])
    pass


def test_protein_to_contig_by_exon_id():
    result = sorted([i.contig for i in protein_to_contig("ENSE00000936617", 1)])
    assert result == ["12"]


def test_protein_to_contig_by_gene_id():
    result = sorted([i.contig for i in protein_to_contig("ENSG00000133703", 1)])
    assert result == ["12"]


def test_protein_to_contig_by_gene_name():
    result = sorted([i.contig for i in protein_to_contig("KRAS", 1)])
    assert result == ["12"]


def test_protein_to_contig_by_transcript_id():
    result = sorted([i.contig for i in protein_to_contig("ENST00000256078", 1)])
    assert result == ["12"]


def test_protein_to_contig_by_transcript_name():
    result = sorted([i.contig for i in protein_to_contig("KRAS-201", 1)])
    assert result == ["12"]


def test_protein_to_exon_pos_strand_by_protein_id():
    pos = protein_to_exon("ENSP00000288135", 1)
    assert isinstance(pos, list)
    assert pos[0].transcript_id, "ENST00000288135"
    assert pos[0].exon_id, "ENSE00000000233"
    assert pos[0].start == 54657957
    assert pos[0].end == 54658081
    assert pos[0].strand == "+"


def test_protein_to_exon_pos_strand_by_protein_id_2():
    pos = protein_to_exon("ENSP00000288135", 4)
    assert isinstance(pos, list)
    assert pos[0].transcript_id, "ENST00000288135"
    assert pos[0].exon_id, "ENSE00000000233"
    assert pos[0].start == 54657957
    assert pos[0].end == 54658081
    assert pos[0].strand == "+"


def test_protein_to_exon_pos_strand_by_protein_id_3():
    pos = protein_to_exon("ENSP00000288135", 4, 6)
    assert isinstance(pos, list)
    assert pos[0].transcript_id, "ENST00000288135"
    assert pos[0].exon_id, "ENSE00000000233"
    assert pos[0].start == 54657957
    assert pos[0].end == 54658081
    assert pos[0].strand == "+"


def test_protein_to_exon_neg_strand_by_protein_id():
    pos = protein_to_exon("ENSP00000308495", 1)
    assert isinstance(pos, list)
    assert pos[0].transcript_id, "ENST00000311936"
    assert pos[0].exon_id, "ENSE00000936617"
    assert pos[0].start == 25245274
    assert pos[0].end == 25245395
    assert pos[0].strand == "-"


def test_protein_to_exon_neg_strand_by_protein_id_2():
    pos = protein_to_exon("ENSP00000308495", 4)
    assert isinstance(pos, list)
    assert pos[0].transcript_id, "ENST00000311936"
    assert pos[0].exon_id, "ENSE00000936617"
    assert pos[0].start == 25245274
    assert pos[0].end == 25245395
    assert pos[0].strand == "-"


def test_protein_to_exon_neg_strand_by_protein_id_3():
    pos = protein_to_exon("ENSP00000308495", 4, 6)
    assert isinstance(pos, list)
    assert pos[0].transcript_id, "ENST00000311936"
    assert pos[0].exon_id, "ENSE00000936617"
    assert pos[0].start == 25245274
    assert pos[0].end == 25245395
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


def test_protein_to_gene_pos_strand_by_protein_id():
    pos = protein_to_gene("ENSP00000288135", 1)
    assert isinstance(pos, list)
    assert pos[0].gene_id, "ENSG00000157404"
    assert pos[0].start == 54658015
    assert pos[0].end == 54658017
    assert pos[0].strand == "+"


def test_protein_to_gene_pos_strand_by_protein_id_2():
    pos = protein_to_gene("ENSP00000288135", 71)
    assert isinstance(pos, list)
    assert pos[0].gene_id, "ENSG00000157404"
    assert pos[0].start == 54695655
    assert pos[0].end == 54695657
    assert pos[0].strand == "+"


def test_protein_to_gene_pos_strand_by_protein_id_3():
    pos = protein_to_gene("ENSP00000288135", 71, 74)
    assert isinstance(pos, list)
    assert pos[0].gene_id, "ENSG00000157404"
    assert pos[0].start == 54695655
    assert pos[0].end == 54695666
    assert pos[0].strand == "+"


def test_protein_to_gene_neg_strand_by_protein_id():
    pos = protein_to_gene("ENSP00000256078", 1)
    assert isinstance(pos, list)
    assert pos[0].gene_id, "ENSG00000133703"
    assert pos[0].start == 25245382
    assert pos[0].end == 25245384
    assert pos[0].strand == "-"


def test_protein_to_gene_neg_strand_by_protein_id_2():
    pos = protein_to_gene("ENSP00000256078", 6)
    assert isinstance(pos, list)
    assert pos[0].gene_id, "ENSG00000133703"
    assert pos[0].start == 25245367
    assert pos[0].end == 25245369
    assert pos[0].strand == "-"


def test_protein_to_gene_neg_strand_by_protein_id_3():
    pos = protein_to_gene("ENSP00000256078", 6, 9)
    assert isinstance(pos, list)
    assert pos[0].gene_id, "ENSG00000133703"
    assert pos[0].start == 25245358
    assert pos[0].end == 25245369
    assert pos[0].strand == "-"


def test_protein_to_gene_by_contig_id():
    # NOTE: I don't see this query being useful because it will return coordinates for every protein on chr12.
    # result = sorted([i.gene_id for i in protein_to_gene("12", 1)])
    pass


def test_protein_to_gene_by_exon_id():
    result = sorted([i.gene_id for i in protein_to_gene("ENSE00000936617", 1)])
    assert result == ["ENSG00000133703"]


def test_protein_to_gene_by_gene_id():
    result = sorted([i.gene_id for i in protein_to_gene("ENSG00000133703", 1)])
    assert result == ["ENSG00000133703"]


def test_protein_to_gene_by_gene_name():
    result = sorted([i.gene_id for i in protein_to_gene("KRAS", 1)])
    assert result == ["ENSG00000133703"]


def test_protein_to_gene_by_transcript_id():
    result = sorted([i.gene_id for i in protein_to_gene("ENST00000256078", 1)])
    assert result == ["ENSG00000133703"]


def test_protein_to_gene_by_transcript_name():
    result = sorted([i.gene_id for i in protein_to_gene("KRAS-201", 1)])
    assert result == ["ENSG00000133703"]


def test_protein_to_protein_pos_strand_by_protein_id():
    pos = protein_to_protein("ENSP00000288135", 1)
    assert isinstance(pos, list)
    assert pos[0].protein_id, "ENSP00000288135"
    assert pos[0].start == 1
    assert pos[0].end == 1
    assert pos[0].strand is None


def test_protein_to_protein_pos_strand_by_protein_id_2():
    pos = protein_to_protein("ENSP00000288135", 12)
    assert isinstance(pos, list)
    assert pos[0].protein_id, "ENSP00000288135"
    assert pos[0].start == 12
    assert pos[0].end == 12
    assert pos[0].strand is None


def test_protein_to_protein_pos_strand_by_protein_id_3():
    pos = protein_to_protein("ENSP00000288135", 12, 21)
    assert isinstance(pos, list)
    assert pos[0].protein_id, "ENSP00000288135"
    assert pos[0].start == 12
    assert pos[0].end == 21
    assert pos[0].strand is None


def test_protein_to_protein_neg_strand_by_protein_id():
    pos = protein_to_protein("ENSP00000260947", 1)
    assert isinstance(pos, list)
    assert pos[0].protein_id, "ENSP00000260947"
    assert pos[0].start == 1
    assert pos[0].end == 1
    assert pos[0].strand is None


def test_protein_to_protein_neg_strand_by_protein_id_2():
    pos = protein_to_protein("ENSP00000260947", 71)
    assert isinstance(pos, list)
    assert pos[0].protein_id, "ENSP00000260947"
    assert pos[0].start == 71
    assert pos[0].end == 71
    assert pos[0].strand is None


def test_protein_to_protein_neg_strand_by_protein_id_3():
    pos = protein_to_protein("ENSP00000260947", 71, 74)
    assert isinstance(pos, list)
    assert pos[0].protein_id, "ENSP00000260947"
    assert pos[0].start == 71
    assert pos[0].end == 74
    assert pos[0].strand is None


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


def test_protein_to_transcript_pos_strand_by_protein_id():
    pos = protein_to_transcript("ENSP00000288135", 1)
    assert isinstance(pos, list)
    assert pos[0].transcript_id, "ENST00000288135"
    assert pos[0].start == 59
    assert pos[0].end == 61
    assert pos[0].strand == "+"


def test_protein_to_transcript_pos_strand_by_protein_id_2():
    pos = protein_to_transcript("ENSP00000288135", 12)
    assert isinstance(pos, list)
    assert pos[0].transcript_id, "ENST00000288135"
    assert pos[0].start == 92
    assert pos[0].end == 94
    assert pos[0].strand == "+"


def test_protein_to_transcript_pos_strand_by_protein_id_3():
    pos = protein_to_transcript("ENSP00000288135", 12, 21)
    assert isinstance(pos, list)
    assert pos[0].transcript_id, "ENST00000288135"
    assert pos[0].start == 92
    assert pos[0].end == 121
    assert pos[0].strand == "+"


def test_protein_to_transcript_neg_strand_by_protein_id():
    pos = protein_to_transcript("ENSP00000260947", 1)
    assert isinstance(pos, list)
    assert pos[0].transcript_id, "ENST00000260947"
    assert pos[0].start == 115
    assert pos[0].end == 117
    assert pos[0].strand == "-"


def test_protein_to_transcript_neg_strand_by_protein_id_2():
    pos = protein_to_transcript("ENSP00000260947", 71)
    assert isinstance(pos, list)
    assert pos[0].transcript_id, "ENST00000260947"
    assert pos[0].start == 325
    assert pos[0].end == 327
    assert pos[0].strand == "-"


def test_protein_to_transcript_neg_strand_by_protein_id_3():
    pos = protein_to_transcript("ENSP00000260947", 71, 74)
    assert isinstance(pos, list)
    assert pos[0].transcript_id, "ENST00000260947"
    assert pos[0].start == 325
    assert pos[0].end == 336
    assert pos[0].strand == "-"


def test_protein_to_transcript_by_contig_id():
    # NOTE: I don't see this query being useful because it will return coordinates for every protein on chr12.
    # result = sorted([i.transcript_id for i in protein_to_transcript("12", 1)])
    pass


def test_protein_to_transcript_by_exon_id():
    result = sorted([i.transcript_id for i in protein_to_transcript("ENSE00000936617", 1)])
    assert result == ["ENST00000256078", "ENST00000311936", "ENST00000556131", "ENST00000557334"]


def test_protein_to_transcript_by_gene_id():
    result = sorted([i.transcript_id for i in protein_to_transcript("ENSG00000133703", 1)])
    assert result == ["ENST00000256078", "ENST00000311936", "ENST00000556131", "ENST00000557334"]


def test_protein_to_transcript_by_gene_name():
    result = sorted([i.transcript_id for i in protein_to_transcript("KRAS", 1)])
    assert result == ["ENST00000256078", "ENST00000311936", "ENST00000556131", "ENST00000557334"]


def test_protein_to_transcript_by_transcript_id():
    result = sorted([i.transcript_id for i in protein_to_transcript("ENST00000256078", 1)])
    assert result == ["ENST00000256078"]


def test_protein_to_transcript_by_transcript_name():
    result = sorted([i.transcript_id for i in protein_to_transcript("KRAS-201", 1)])
    assert result == ["ENST00000256078"]


def test_transcript_to_cds_pos_strand_by_transcript_id():
    pos = transcript_to_cds("ENST00000288135", 98)
    assert isinstance(pos, list)
    assert pos[0].transcript_id, "ENST00000288135"
    assert pos[0].start == 40
    assert pos[0].end == 40
    assert pos[0].strand == "+"


def test_transcript_to_cds_pos_strand_by_transcript_id_2():
    pos = transcript_to_cds("ENST00000288135", 100)
    assert isinstance(pos, list)
    assert pos[0].transcript_id, "ENST00000288135"
    assert pos[0].start == 42
    assert pos[0].end == 42
    assert pos[0].strand == "+"


def test_transcript_to_cds_pos_strand_by_transcript_id_3():
    pos = transcript_to_cds("ENST00000288135", 98, 100)
    assert isinstance(pos, list)
    assert pos[0].transcript_id, "ENST00000288135"
    assert pos[0].start == 40
    assert pos[0].end == 42
    assert pos[0].strand == "+"


def test_transcript_to_cds_pos_strand_by_transcript_name():
    pos = transcript_to_cds("KIT-201", 98)
    assert isinstance(pos, list)
    assert pos[0].transcript_id, "ENST00000288135"
    assert pos[0].start == 40
    assert pos[0].end == 40
    assert pos[0].strand == "+"


def test_transcript_to_cds_pos_strand_by_transcript_name_2():
    pos = transcript_to_cds("KIT-201", 100)
    assert isinstance(pos, list)
    assert pos[0].transcript_id, "ENST00000288135"
    assert pos[0].start == 42
    assert pos[0].end == 42
    assert pos[0].strand == "+"


def test_transcript_to_cds_pos_strand_by_transcript_name_3():
    pos = transcript_to_cds("KIT-201", 98, 100)
    assert isinstance(pos, list)
    assert pos[0].transcript_id, "ENST00000288135"
    assert pos[0].start == 40
    assert pos[0].end == 42
    assert pos[0].strand == "+"


def test_transcript_to_cds_neg_strand_by_transcript_id():
    pos = transcript_to_cds("ENST00000256078", 191)
    assert isinstance(pos, list)
    assert pos[0].transcript_id, "ENST00000256078"
    assert pos[0].start == 1
    assert pos[0].end == 1
    assert pos[0].strand == "-"


def test_transcript_to_cds_neg_strand_by_transcript_id_2():
    pos = transcript_to_cds("ENST00000256078", 301)
    assert isinstance(pos, list)
    assert pos[0].transcript_id, "ENST00000256078"
    assert pos[0].start == 111
    assert pos[0].end == 111
    assert pos[0].strand == "-"


def test_transcript_to_cds_neg_strand_by_transcript_id_3():
    pos = transcript_to_cds("ENST00000256078", 301, 323)
    assert isinstance(pos, list)
    assert pos[0].transcript_id, "ENST00000256078"
    assert pos[0].start == 111
    assert pos[0].end == 133
    assert pos[0].strand == "-"


def test_transcript_to_cds_neg_strand_by_transcript_name():
    pos = transcript_to_cds("KRAS-201", 191)
    assert isinstance(pos, list)
    assert pos[0].transcript_id, "ENST00000256078"
    assert pos[0].start == 1
    assert pos[0].end == 1
    assert pos[0].strand == "-"


def test_transcript_to_cds_neg_strand_by_transcript_name_2():
    pos = transcript_to_cds("KRAS-201", 301)
    assert isinstance(pos, list)
    assert pos[0].transcript_id, "ENST00000256078"
    assert pos[0].start == 111
    assert pos[0].end == 111
    assert pos[0].strand == "-"


def test_transcript_to_cds_neg_strand_by_transcript_name_3():
    pos = transcript_to_cds("KRAS-201", 301, 323)
    assert isinstance(pos, list)
    assert pos[0].transcript_id, "ENST00000256078"
    assert pos[0].start == 111
    assert pos[0].end == 133
    assert pos[0].strand == "-"


def test_transcript_to_cds_by_contig_id():
    # NOTE: I don't see this query being useful because it will return coordinates for every transcript on chr12.
    # result = sorted([i.transcript_id for i in transcript_to_cds("12", 301)])
    pass


def test_transcript_to_cds_by_exon_id():
    result = sorted([i.transcript_id for i in transcript_to_cds("ENSE00000936617", 301)])
    assert result == ["ENST00000256078", "ENST00000311936", "ENST00000556131", "ENST00000557334"]


def test_transcript_to_cds_by_gene_id():
    result = sorted([i.transcript_id for i in transcript_to_cds("ENSG00000133703", 301)])
    assert result == ["ENST00000256078", "ENST00000311936", "ENST00000556131", "ENST00000557334"]


def test_transcript_to_cds_by_gene_name():
    result = sorted([i.transcript_id for i in transcript_to_cds("KRAS", 301)])
    assert result == ["ENST00000256078", "ENST00000311936", "ENST00000556131", "ENST00000557334"]


def test_transcript_to_cds_by_protein_id():
    result = sorted([i.transcript_id for i in transcript_to_cds("ENSP00000256078", 301)])
    assert result == ["ENST00000256078"]


def test_transcript_to_contig_pos_strand_by_transcript_id():
    pos = transcript_to_contig("ENST00000341165", 1)
    assert isinstance(pos, list)
    assert pos[0].contig == "17"
    assert pos[0].start == 43170481
    assert pos[0].end == 43170481
    assert pos[0].strand == "+"


def test_transcript_to_contig_pos_strand_by_transcript_id_2():
    pos = transcript_to_contig("ENST00000341165", 731)
    assert isinstance(pos, list)
    assert pos[0].contig == "17"
    assert pos[0].start == 43189698
    assert pos[0].end == 43189698
    assert pos[0].strand == "+"


def test_transcript_to_contig_pos_strand_by_transcript_id_3():
    pos = transcript_to_contig("ENST00000341165", 731, 817)
    assert isinstance(pos, list)
    assert pos[0].contig == "17"
    assert pos[0].start == 43189698
    assert pos[0].end == 43189784
    assert pos[0].strand == "+"


def test_transcript_to_contig_pos_strand_by_transcript_name():
    pos = transcript_to_contig("NBR1-201", 1)
    assert isinstance(pos, list)
    assert pos[0].contig == "17"
    assert pos[0].start == 43170481
    assert pos[0].end == 43170481
    assert pos[0].strand == "+"


def test_transcript_to_contig_pos_strand_by_transcript_name_2():
    pos = transcript_to_contig("NBR1-201", 731)
    assert isinstance(pos, list)
    assert pos[0].contig == "17"
    assert pos[0].start == 43189698
    assert pos[0].end == 43189698
    assert pos[0].strand == "+"


def test_transcript_to_contig_pos_strand_by_transcript_name_3():
    pos = transcript_to_contig("NBR1-201", 731, 817)
    assert isinstance(pos, list)
    assert pos[0].contig == "17"
    assert pos[0].start == 43189698
    assert pos[0].end == 43189784
    assert pos[0].strand == "+"


def test_transcript_to_contig_neg_strand_by_transcript_id():
    pos = transcript_to_contig("ENST00000256078", 1)
    assert isinstance(pos, list)
    assert pos[0].contig == "12"
    assert pos[0].start == 25250929
    assert pos[0].end == 25250929
    assert pos[0].strand == "+"


def test_transcript_to_contig_neg_strand_by_transcript_id_2():
    pos = transcript_to_contig("ENST00000256078", 466)
    assert isinstance(pos, list)
    assert pos[0].contig == "12"
    assert pos[0].start == 25227248
    assert pos[0].end == 25227248
    assert pos[0].strand == "+"


def test_transcript_to_contig_neg_strand_by_transcript_id_3():
    pos = transcript_to_contig("ENST00000256078", 466, 501)
    assert isinstance(pos, list)
    assert pos[0].contig == "12"
    assert pos[0].start == 25225753
    assert pos[0].end == 25227248
    assert pos[0].strand == "+"


def test_transcript_to_contig_neg_strand_by_transcript_name():
    pos = transcript_to_contig("KRAS-201", 1)
    assert isinstance(pos, list)
    assert pos[0].contig == "12"
    assert pos[0].start == 25250929
    assert pos[0].end == 25250929
    assert pos[0].strand == "+"


def test_transcript_to_contig_neg_strand_by_transcript_name_2():
    pos = transcript_to_contig("KRAS-201", 466)
    assert isinstance(pos, list)
    assert pos[0].contig == "12"
    assert pos[0].start == 25227248
    assert pos[0].end == 25227248
    assert pos[0].strand == "+"


def test_transcript_to_contig_neg_strand_by_transcript_name_3():
    pos = transcript_to_contig("KRAS-201", 466, 501)
    assert isinstance(pos, list)
    assert pos[0].contig == "12"
    assert pos[0].start == 25225753
    assert pos[0].end == 25227248
    assert pos[0].strand == "+"


def test_transcript_to_contig_by_refseq_id():
    pos = transcript_to_contig("NM_001004491", 1)
    assert isinstance(pos, list)
    assert pos[0].contig == "1"
    assert pos[0].start == 247965233
    assert pos[0].end == 247965233
    assert pos[0].strand == "+"


def test_transcript_to_contig_by_contig_id():
    # NOTE: I don't see this query being useful because it will return coordinates for every transcript on chr12.
    # result = sorted([i.contig for i in transcript_to_contig("12", 301)])
    pass


def test_transcript_to_contig_by_exon_id():
    result = sorted([i.contig for i in transcript_to_contig("ENSE00000936617", 301)])
    assert result == ["12", "12", "12"]


def test_transcript_to_contig_by_gene_id():
    result = sorted([i.contig for i in transcript_to_contig("ENSG00000133703", 301)])
    assert result == ["12", "12", "12"]


def test_transcript_to_contig_by_gene_name():
    result = sorted([i.contig for i in transcript_to_contig("KRAS", 301)])
    assert result == ["12", "12", "12"]


def test_transcript_to_contig_by_protein_id():
    result = sorted([i.contig for i in transcript_to_contig("ENSP00000256078", 301)])
    assert result == ["12"]


def test_transcript_to_exon_pos_strand_by_transcript_id():
    pos = transcript_to_exon("ENST00000380152", 1)
    assert isinstance(pos, list)
    assert pos[0].transcript_id, "ENST00000380152"
    assert pos[0].exon_id, "ENSE00001184784"
    assert pos[0].start == 32315474
    assert pos[0].end == 32315667
    assert pos[0].strand == "+"


def test_transcript_to_exon_pos_strand_by_transcript_id_2():
    pos = transcript_to_exon("ENST00000380152", 500)
    assert isinstance(pos, list)
    assert pos[0].transcript_id, "ENST00000380152"
    assert pos[0].exon_id, "ENSE00003666217"
    assert pos[0].start == 32319077
    assert pos[0].end == 32319325
    assert pos[0].strand == "+"


def test_transcript_to_exon_pos_strand_by_transcript_id_3():
    pos = transcript_to_exon("ENST00000380152", 500, 510)
    assert isinstance(pos, list)
    assert pos[0].transcript_id, "ENST00000380152"
    assert pos[0].exon_id, "ENSE00003666217"
    assert pos[0].start == 32319077
    assert pos[0].end == 32319325
    assert pos[0].strand == "+"


def test_transcript_to_exon_pos_strand_by_transcript_name():
    pos = transcript_to_exon("BRCA2-201", 1)
    assert isinstance(pos, list)
    assert pos[0].transcript_id, "ENST00000380152"
    assert pos[0].exon_id, "ENSE00001184784"
    assert pos[0].start == 32315474
    assert pos[0].end == 32315667
    assert pos[0].strand == "+"


def test_transcript_to_exon_pos_strand_by_transcript_name_2():
    pos = transcript_to_exon("BRCA2-201", 500)
    assert isinstance(pos, list)
    assert pos[0].transcript_id, "ENST00000380152"
    assert pos[0].exon_id, "ENSE00003666217"
    assert pos[0].start == 32319077
    assert pos[0].end == 32319325
    assert pos[0].strand == "+"


def test_transcript_to_exon_pos_strand_by_transcript_name_3():
    pos = transcript_to_exon("BRCA2-201", 500, 510)
    assert isinstance(pos, list)
    assert pos[0].transcript_id, "ENST00000380152"
    assert pos[0].exon_id, "ENSE00003666217"
    assert pos[0].start == 32319077
    assert pos[0].end == 32319325
    assert pos[0].strand == "+"


def test_transcript_to_exon_neg_strand_by_transcript_id():
    pos = transcript_to_exon("ENST00000256078", 1)
    assert isinstance(pos, list)
    assert pos[0].transcript_id, "ENST00000256078"
    assert pos[0].exon_id, "ENSE00000000028"
    assert pos[0].start == 25250751
    assert pos[0].end == 25250929
    assert pos[0].strand == "-"


def test_transcript_to_exon_neg_strand_by_transcript_id_2():
    pos = transcript_to_exon("ENST00000256078", 400)
    assert isinstance(pos, list)
    assert pos[0].transcript_id, "ENST00000256078"
    assert pos[0].exon_id, "ENSE00001719809"
    assert pos[0].start == 25227234
    assert pos[0].end == 25227412
    assert pos[0].strand == "-"


def test_transcript_to_exon_neg_strand_by_transcript_id_3():
    pos = transcript_to_exon("ENST00000256078", 400, 420)
    assert isinstance(pos, list)
    assert pos[0].transcript_id, "ENST00000256078"
    assert pos[0].exon_id, "ENSE00001719809"
    assert pos[0].start == 25227234
    assert pos[0].end == 25227412
    assert pos[0].strand == "-"


def test_transcript_to_exon_neg_strand_by_transcript_name():
    pos = transcript_to_exon("KRAS-201", 1)
    assert isinstance(pos, list)
    assert pos[0].transcript_id, "ENST00000256078"
    assert pos[0].exon_id, "ENSE00000000028"
    assert pos[0].start == 25250751
    assert pos[0].end == 25250929
    assert pos[0].strand == "-"


def test_transcript_to_exon_neg_strand_by_transcript_name_2():
    pos = transcript_to_exon("KRAS-201", 400)
    assert isinstance(pos, list)
    assert pos[0].transcript_id, "ENST00000256078"
    assert pos[0].exon_id, "ENSE00001719809"
    assert pos[0].start == 25227234
    assert pos[0].end == 25227412
    assert pos[0].strand == "-"


def test_transcript_to_exon_neg_strand_by_transcript_name_3():
    pos = transcript_to_exon("KRAS-201", 400, 420)
    assert isinstance(pos, list)
    assert pos[0].transcript_id, "ENST00000256078"
    assert pos[0].exon_id, "ENSE00001719809"
    assert pos[0].start == 25227234
    assert pos[0].end == 25227412
    assert pos[0].strand == "-"


def test_transcript_to_exon_by_contig_id():
    # NOTE: I don't see this query being useful because it will return coordinates for every transcript on chr12.
    # result = sorted([i.exon_id for i in transcript_to_exon("12", 1)])
    pass


def test_transcript_to_exon_by_exon_id():
    result = sorted([i.exon_id for i in transcript_to_exon("ENSE00000936617", 1)])
    assert result == ["ENSE00000000028", "ENSE00001189804", "ENSE00002446502", "ENSE00002530521"]


def test_transcript_to_exon_by_gene_id():
    result = sorted([i.exon_id for i in transcript_to_exon("ENSG00000133703", 1)])
    assert result == ["ENSE00000000028", "ENSE00001189804", "ENSE00002446502", "ENSE00002530521"]


def test_transcript_to_exon_by_gene_name():
    result = sorted([i.exon_id for i in transcript_to_exon("KRAS", 1)])
    assert result == ["ENSE00000000028", "ENSE00001189804", "ENSE00002446502", "ENSE00002530521"]


def test_transcript_to_exon_by_protein_id():
    result = sorted([i.exon_id for i in transcript_to_exon("ENSP00000256078", 1)])
    assert result == ["ENSE00000000028"]


def test_transcript_to_gene_pos_strand_by_transcript_id():
    pos = transcript_to_gene("ENST00000341165", 1)
    assert isinstance(pos, list)
    assert pos[0].gene_id, "ENSG00000188554"
    assert pos[0].start == 43170481
    assert pos[0].end == 43170481
    assert pos[0].strand == "+"


def test_transcript_to_gene_pos_strand_by_transcript_id_2():
    pos = transcript_to_gene("ENST00000341165", 731)
    assert isinstance(pos, list)
    assert pos[0].gene_id, "ENSG00000188554"
    assert pos[0].start == 43189698
    assert pos[0].end == 43189698
    assert pos[0].strand == "+"


def test_transcript_to_gene_pos_strand_by_transcript_id_3():
    pos = transcript_to_gene("ENST00000341165", 731, 817)
    assert isinstance(pos, list)
    assert pos[0].gene_id, "ENSG00000188554"
    assert pos[0].start == 43189698
    assert pos[0].end == 43189784
    assert pos[0].strand == "+"


def test_transcript_to_gene_pos_strand_by_transcript_name():
    pos = transcript_to_gene("NBR1-201", 1)
    assert isinstance(pos, list)
    assert pos[0].gene_id, "ENSG00000188554"
    assert pos[0].start == 43170481
    assert pos[0].end == 43170481
    assert pos[0].strand == "+"


def test_transcript_to_gene_pos_strand_by_transcript_name_2():
    pos = transcript_to_gene("NBR1-201", 731)
    assert isinstance(pos, list)
    assert pos[0].gene_id, "ENSG00000188554"
    assert pos[0].start == 43189698
    assert pos[0].end == 43189698
    assert pos[0].strand == "+"


def test_transcript_to_gene_pos_strand_by_transcript_name_3():
    pos = transcript_to_gene("NBR1-201", 731, 817)
    assert isinstance(pos, list)
    assert pos[0].gene_id, "ENSG00000188554"
    assert pos[0].start == 43189698
    assert pos[0].end == 43189784
    assert pos[0].strand == "+"


def test_transcript_to_gene_neg_strand_by_transcript_id():
    pos = transcript_to_gene("ENST00000256078", 1)
    assert isinstance(pos, list)
    assert pos[0].gene_id, "ENSG00000133703"
    assert pos[0].start == 25250929
    assert pos[0].end == 25250929
    assert pos[0].strand == "-"


def test_transcript_to_gene_neg_strand_by_transcript_id_2():
    pos = transcript_to_gene("ENST00000256078", 466)
    assert isinstance(pos, list)
    assert pos[0].gene_id, "ENSG00000133703"
    assert pos[0].start == 25227248
    assert pos[0].end == 25227248
    assert pos[0].strand == "-"


def test_transcript_to_gene_neg_strand_by_transcript_id_3():
    pos = transcript_to_gene("ENST00000256078", 466, 501)
    assert isinstance(pos, list)
    assert pos[0].gene_id, "ENSG00000133703"
    assert pos[0].start == 25225753
    assert pos[0].end == 25227248
    assert pos[0].strand == "-"


def test_transcript_to_gene_neg_strand_by_transcript_name():
    pos = transcript_to_gene("KRAS-201", 1)
    assert isinstance(pos, list)
    assert pos[0].gene_id, "ENSG00000133703"
    assert pos[0].start == 25250929
    assert pos[0].end == 25250929
    assert pos[0].strand == "-"


def test_transcript_to_gene_neg_strand_by_transcript_name_2():
    pos = transcript_to_gene("KRAS-201", 466)
    assert isinstance(pos, list)
    assert pos[0].gene_id, "ENSG00000133703"
    assert pos[0].start == 25227248
    assert pos[0].end == 25227248
    assert pos[0].strand == "-"


def test_transcript_to_gene_neg_strand_by_transcript_name_3():
    pos = transcript_to_gene("KRAS-201", 466, 501)
    assert isinstance(pos, list)
    assert pos[0].gene_id, "ENSG00000133703"
    assert pos[0].start == 25225753
    assert pos[0].end == 25227248
    assert pos[0].strand == "-"


def test_transcript_to_gene_by_refseq_id():
    pos = transcript_to_gene("NM_001004491", 1)
    assert isinstance(pos, list)
    assert pos[0].gene_id, "ENSG00000187080"
    assert pos[0].start == 247965233
    assert pos[0].end == 247965233
    assert pos[0].strand == "+"


def test_transcript_to_gene_by_contig_id():
    # NOTE: I don't see this query being useful because it will return coordinates for every transcript on chr12.
    # result = sorted([i.gene_id for i in transcript_to_gene("12", 1)])
    pass


def test_transcript_to_gene_by_exon_id():
    result = sorted([i.gene_id for i in transcript_to_gene("ENSE00000936617", 1)])
    assert result == ["ENSG00000133703", "ENSG00000133703"]


def test_transcript_to_gene_by_gene_id():
    result = sorted([i.gene_id for i in transcript_to_gene("ENSG00000133703", 1)])
    assert result == ["ENSG00000133703", "ENSG00000133703"]


def test_transcript_to_gene_by_gene_name():
    result = sorted([i.gene_id for i in transcript_to_gene("KRAS", 1)])
    assert result == ["ENSG00000133703", "ENSG00000133703"]


def test_transcript_to_gene_by_protein_id():
    result = sorted([i.gene_id for i in transcript_to_gene("ENSP00000256078", 1)])
    assert result == ["ENSG00000133703"]


def test_transcript_to_protein_pos_strand_by_transcript_id():
    pos = transcript_to_protein("ENST00000288135", 98)
    assert isinstance(pos, list)
    assert pos[0].protein_id, "ENSP00000288135"
    assert pos[0].start == 14
    assert pos[0].end == 14
    assert pos[0].strand is None


def test_transcript_to_protein_pos_strand_by_transcript_id_2():
    pos = transcript_to_protein("ENST00000288135", 128)
    assert isinstance(pos, list)
    assert pos[0].protein_id, "ENSP00000288135"
    assert pos[0].start == 24
    assert pos[0].end == 24
    assert pos[0].strand is None


def test_transcript_to_protein_pos_strand_by_transcript_id_3():
    pos = transcript_to_protein("ENST00000288135", 128, 145)
    assert isinstance(pos, list)
    assert pos[0].protein_id, "ENSP00000288135"
    assert pos[0].start == 24
    assert pos[0].end == 29
    assert pos[0].strand is None


def test_transcript_to_protein_pos_strand_by_transcript_name():
    pos = transcript_to_protein("KIT-201", 98)
    assert isinstance(pos, list)
    assert pos[0].protein_id, "ENSP00000288135"
    assert pos[0].start == 14
    assert pos[0].end == 14
    assert pos[0].strand is None


def test_transcript_to_protein_pos_strand_by_transcript_name_2():
    pos = transcript_to_protein("KIT-201", 128)
    assert isinstance(pos, list)
    assert pos[0].protein_id, "ENSP00000288135"
    assert pos[0].start == 24
    assert pos[0].end == 24
    assert pos[0].strand is None


def test_transcript_to_protein_pos_strand_by_transcript_name_3():
    pos = transcript_to_protein("KIT-201", 128, 145)
    assert isinstance(pos, list)
    assert pos[0].protein_id, "ENSP00000288135"
    assert pos[0].start == 24
    assert pos[0].end == 29
    assert pos[0].strand is None


def test_transcript_to_protein_neg_strand_by_transcript_id():
    pos = transcript_to_protein("ENST00000260947", 115)
    assert isinstance(pos, list)
    assert pos[0].protein_id, "ENSP00000260947"
    assert pos[0].start == 1
    assert pos[0].end == 1
    assert pos[0].strand is None


def test_transcript_to_protein_neg_strand_by_transcript_id_2():
    pos = transcript_to_protein("ENST00000260947", 346)
    assert isinstance(pos, list)
    assert pos[0].protein_id, "ENSP00000260947"
    assert pos[0].start == 78
    assert pos[0].end == 78
    assert pos[0].strand is None


def test_transcript_to_protein_neg_strand_by_transcript_id_3():
    pos = transcript_to_protein("ENST00000260947", 345, 356)
    assert isinstance(pos, list)
    assert pos[0].protein_id, "ENSP00000260947"
    assert pos[0].start == 77
    assert pos[0].end == 81
    assert pos[0].strand is None


def test_transcript_to_protein_neg_strand_by_transcript_name():
    pos = transcript_to_protein("BARD1-201", 115)
    assert isinstance(pos, list)
    assert pos[0].protein_id, "ENSP00000260947"
    assert pos[0].start == 1
    assert pos[0].end == 1
    assert pos[0].strand is None


def test_transcript_to_protein_neg_strand_by_transcript_name_2():
    pos = transcript_to_protein("BARD1-201", 346)
    assert isinstance(pos, list)
    assert pos[0].protein_id, "ENSP00000260947"
    assert pos[0].start == 78
    assert pos[0].end == 78
    assert pos[0].strand is None


def test_transcript_to_protein_neg_strand_by_transcript_name_3():
    pos = transcript_to_protein("BARD1-201", 345, 356)
    assert isinstance(pos, list)
    assert pos[0].protein_id, "ENSP00000260947"
    assert pos[0].start == 77
    assert pos[0].end == 81
    assert pos[0].strand is None


def test_transcript_to_protein_by_contig_id():
    # NOTE: I don't see this query being useful because it will return coordinates for every transcript on chr12.
    # result = sorted([i.protein_id for i in transcript_to_protein("12", 1)])
    pass


def test_transcript_to_protein_by_exon_id():
    result = sorted([i.protein_id for i in transcript_to_protein("ENSE00000936617", 301)])
    assert result == ["ENSP00000256078", "ENSP00000308495", "ENSP00000451856", "ENSP00000452512"]


def test_transcript_to_protein_by_gene_id():
    result = sorted([i.protein_id for i in transcript_to_protein("ENSG00000133703", 301)])
    assert result == ["ENSP00000256078", "ENSP00000308495", "ENSP00000451856", "ENSP00000452512"]


def test_transcript_to_protein_by_gene_name():
    result = sorted([i.protein_id for i in transcript_to_protein("KRAS", 301)])
    assert result == ["ENSP00000256078", "ENSP00000308495", "ENSP00000451856", "ENSP00000452512"]


def test_transcript_to_protein_by_protein_id():
    result = sorted([i.protein_id for i in transcript_to_protein("ENSP00000256078", 301)])
    assert result == ["ENSP00000256078"]


def test_transcript_to_transcript_pos_strand_by_transcript_id():
    pos = transcript_to_transcript("ENST00000288135", 1)
    assert isinstance(pos, list)
    assert pos[0].transcript_id, "ENST00000288135"
    assert pos[0].start == 1
    assert pos[0].end == 1
    assert pos[0].strand == "+"


def test_transcript_to_transcript_pos_strand_by_transcript_id_2():
    pos = transcript_to_transcript("ENST00000288135", 100)
    assert isinstance(pos, list)
    assert pos[0].transcript_id, "ENST00000288135"
    assert pos[0].start == 100
    assert pos[0].end == 100
    assert pos[0].strand == "+"


def test_transcript_to_transcript_pos_strand_by_transcript_id_3():
    pos = transcript_to_transcript("ENST00000288135", 98, 100)
    assert isinstance(pos, list)
    assert pos[0].transcript_id, "ENST00000288135"
    assert pos[0].start == 98
    assert pos[0].end == 100
    assert pos[0].strand == "+"


def test_transcript_to_transcript_pos_strand_by_transcript_name():
    pos = transcript_to_transcript("KIT-201", 1)
    assert isinstance(pos, list)
    assert pos[0].transcript_id, "ENST00000288135"
    assert pos[0].start == 1
    assert pos[0].end == 1
    assert pos[0].strand == "+"


def test_transcript_to_transcript_pos_strand_by_transcript_name_2():
    pos = transcript_to_transcript("KIT-201", 100)
    assert isinstance(pos, list)
    assert pos[0].transcript_id, "ENST00000288135"
    assert pos[0].start == 100
    assert pos[0].end == 100
    assert pos[0].strand == "+"


def test_transcript_to_transcript_pos_strand_by_transcript_name_3():
    pos = transcript_to_transcript("KIT-201", 98, 100)
    assert isinstance(pos, list)
    assert pos[0].transcript_id, "ENST00000288135"
    assert pos[0].start == 98
    assert pos[0].end == 100
    assert pos[0].strand == "+"


def test_transcript_to_transcript_neg_strand_by_transcript_id():
    pos = transcript_to_transcript("ENST00000256078", 1)
    assert isinstance(pos, list)
    assert pos[0].transcript_id, "ENST00000256078"
    assert pos[0].start == 1
    assert pos[0].end == 1
    assert pos[0].strand == "-"


def test_transcript_to_transcript_neg_strand_by_transcript_id_2():
    pos = transcript_to_transcript("ENST00000256078", 111)
    assert isinstance(pos, list)
    assert pos[0].transcript_id, "ENST00000256078"
    assert pos[0].start == 111
    assert pos[0].end == 111
    assert pos[0].strand == "-"


def test_transcript_to_transcript_neg_strand_by_transcript_id_3():
    pos = transcript_to_transcript("ENST00000256078", 111, 133)
    assert isinstance(pos, list)
    assert pos[0].transcript_id, "ENST00000256078"
    assert pos[0].start == 111
    assert pos[0].end == 133
    assert pos[0].strand == "-"


def test_transcript_to_transcript_neg_strand_by_transcript_name():
    pos = transcript_to_transcript("KRAS-201", 1)
    assert isinstance(pos, list)
    assert pos[0].transcript_id, "ENST00000256078"
    assert pos[0].start == 1
    assert pos[0].end == 1
    assert pos[0].strand == "-"


def test_transcript_to_transcript_neg_strand_by_transcript_name_2():
    pos = transcript_to_transcript("KRAS-201", 111)
    assert isinstance(pos, list)
    assert pos[0].transcript_id, "ENST00000256078"
    assert pos[0].start == 111
    assert pos[0].end == 111
    assert pos[0].strand == "-"


def test_transcript_to_transcript_neg_strand_by_transcript_name_3():
    pos = transcript_to_transcript("KRAS-201", 111, 133)
    assert isinstance(pos, list)
    assert pos[0].transcript_id, "ENST00000256078"
    assert pos[0].start == 111
    assert pos[0].end == 133
    assert pos[0].strand == "-"


def test_transcript_to_transcript_by_contig_id():
    # NOTE: I don't see this query being useful because it will return coordinates for essentially every transcript on chr12.
    # result = sorted([i.transcript_id for i in transcript_to_transcript("12", 301)])
    pass


def test_transcript_to_transcript_by_exon_id():
    result = sorted([i.transcript_id for i in transcript_to_transcript("ENSE00000936617", 301)])
    assert result == ["ENST00000256078", "ENST00000311936", "ENST00000556131", "ENST00000557334"]


def test_transcript_to_transcript_by_gene_id():
    result = sorted([i.transcript_id for i in transcript_to_transcript("ENSG00000133703", 301)])
    assert result == ["ENST00000256078", "ENST00000311936", "ENST00000556131", "ENST00000557334"]


def test_transcript_to_transcript_by_gene_name():
    result = sorted([i.transcript_id for i in transcript_to_transcript("KRAS", 301)])
    assert result == ["ENST00000256078", "ENST00000311936", "ENST00000556131", "ENST00000557334"]


def test_transcript_to_transcript_by_protein_id():
    result = sorted([i.transcript_id for i in transcript_to_transcript("ENSP00000256078", 301)])
    assert result == ["ENST00000256078"]


def test_best_transcript_only():
    pos = gene_to_transcript("KIT", 54698359, best_only=False)
    assert isinstance(pos, list)
    assert len(pos), 2

    pos = gene_to_transcript("KIT", 54698359, best_only=True)
    assert isinstance(pos, list)
    assert len(pos), 1


def test_get_map_func():
    assert get_map_func(CDS, CDS) == cds_to_cds
    assert get_map_func(CDS, CONTIG) == cds_to_contig
    assert get_map_func(CDS, EXON) == cds_to_exon
    assert get_map_func(CDS, GENE) == cds_to_gene
    assert get_map_func(CDS, PROTEIN) == cds_to_protein
    assert get_map_func(CDS, TRANSCRIPT) == cds_to_transcript
    assert get_map_func(CONTIG, CDS) == contig_to_cds
    assert get_map_func(CONTIG, CONTIG) == contig_to_contig
    assert get_map_func(CONTIG, EXON) == contig_to_exon
    assert get_map_func(CONTIG, GENE) == contig_to_gene
    assert get_map_func(CONTIG, PROTEIN) == contig_to_protein
    assert get_map_func(CONTIG, TRANSCRIPT) == contig_to_transcript
    assert get_map_func(EXON, CDS) == exon_to_cds
    assert get_map_func(EXON, CONTIG) == exon_to_contig
    assert get_map_func(EXON, EXON) == exon_to_exon
    assert get_map_func(EXON, GENE) == exon_to_gene
    assert get_map_func(EXON, PROTEIN) == exon_to_protein
    assert get_map_func(EXON, TRANSCRIPT) == exon_to_transcript
    assert get_map_func(GENE, CDS) == gene_to_cds
    assert get_map_func(GENE, CONTIG) == gene_to_contig
    assert get_map_func(GENE, EXON) == gene_to_exon
    assert get_map_func(GENE, GENE) == gene_to_gene
    assert get_map_func(GENE, PROTEIN) == gene_to_protein
    assert get_map_func(GENE, TRANSCRIPT) == gene_to_transcript
    assert get_map_func(PROTEIN, CDS) == protein_to_cds
    assert get_map_func(PROTEIN, CONTIG) == protein_to_contig
    assert get_map_func(PROTEIN, EXON) == protein_to_exon
    assert get_map_func(PROTEIN, GENE) == protein_to_gene
    assert get_map_func(PROTEIN, PROTEIN) == protein_to_protein
    assert get_map_func(PROTEIN, TRANSCRIPT) == protein_to_transcript
    assert get_map_func(TRANSCRIPT, CDS) == transcript_to_cds
    assert get_map_func(TRANSCRIPT, CONTIG) == transcript_to_contig
    assert get_map_func(TRANSCRIPT, EXON) == transcript_to_exon
    assert get_map_func(TRANSCRIPT, GENE) == transcript_to_gene
    assert get_map_func(TRANSCRIPT, PROTEIN) == transcript_to_protein
    assert get_map_func(TRANSCRIPT, TRANSCRIPT) == transcript_to_transcript


def test_transcript_ids_with_exon():
    assert _transcript_ids_with_exon("ENSE00001719809") == ["ENST00000256078", "ENST00000311936"]
    assert _transcript_ids_with_exon("ENSE00001719809", 1) == []
    assert _transcript_ids_with_exon("ENSE00001719809", 3) == ["ENST00000256078", "ENST00000311936"]
