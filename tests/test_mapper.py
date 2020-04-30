import unittest

from ensembl_map import set_ensembl_release
from ensembl_map.mapper import (
    cds_to_exon,
    cds_to_gene,
    cds_to_protein,
    cds_to_transcript,
    contig_to_cds,
    contig_to_gene,
    contig_to_protein,
    contig_to_transcript,
    exon_to_cds,
    exon_to_gene,
    exon_to_protein,
    exon_to_transcript,
    gene_to_cds,
    gene_to_exon,
    gene_to_protein,
    gene_to_transcript,
    protein_to_cds,
    protein_to_exon,
    protein_to_gene,
    protein_to_transcript,
    transcript_to_cds,
    transcript_to_exon,
    transcript_to_gene,
    transcript_to_protein,
)


# TODO: Mock Cache object instead of needing to download annotations to run tests.
set_ensembl_release(release=99, species="human", download_if_missing=True)


class TestCdsToExon(unittest.TestCase):
    def test_cds_to_exon_pos_strand_by_transcript_id(self):
        pos = cds_to_exon("ENST00000380152", 1)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0].transcript_id, "ENST00000380152")
        self.assertEqual(pos[0].exon_id, "ENSE00001484009")
        self.assertEqual(pos[0].start, 32316422)
        self.assertEqual(pos[0].end, 32316527)
        self.assertEqual(pos[0].strand, "+")

    def test_cds_to_exon_pos_strand_by_transcript_id_2(self):
        pos = cds_to_exon("ENST00000380152", 360)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0].transcript_id, "ENST00000380152")
        self.assertEqual(pos[0].exon_id, "ENSE00003659301")
        self.assertEqual(pos[0].start, 32325076)
        self.assertEqual(pos[0].end, 32325184)
        self.assertEqual(pos[0].strand, "+")

    def test_cds_to_exon_pos_strand_by_transcript_id_3(self):
        pos = cds_to_exon("ENST00000380152", 360, 380)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0].transcript_id, "ENST00000380152")
        self.assertEqual(pos[0].exon_id, "ENSE00003659301")
        self.assertEqual(pos[0].start, 32325076)
        self.assertEqual(pos[0].end, 32325184)
        self.assertEqual(pos[0].strand, "+")

    def test_cds_to_exon_pos_strand_by_transcript_name(self):
        pos = cds_to_exon("BRCA2-201", 1)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0].transcript_id, "ENST00000380152")
        self.assertEqual(pos[0].exon_id, "ENSE00001484009")
        self.assertEqual(pos[0].start, 32316422)
        self.assertEqual(pos[0].end, 32316527)
        self.assertEqual(pos[0].strand, "+")

    def test_cds_to_exon_pos_strand_by_transcript_name_2(self):
        pos = cds_to_exon("BRCA2-201", 360)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0].transcript_id, "ENST00000380152")
        self.assertEqual(pos[0].exon_id, "ENSE00003659301")
        self.assertEqual(pos[0].start, 32325076)
        self.assertEqual(pos[0].end, 32325184)
        self.assertEqual(pos[0].strand, "+")

    def test_cds_to_exon_pos_strand_by_transcript_name_3(self):
        pos = cds_to_exon("BRCA2-201", 360, 380)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0].transcript_id, "ENST00000380152")
        self.assertEqual(pos[0].exon_id, "ENSE00003659301")
        self.assertEqual(pos[0].start, 32325076)
        self.assertEqual(pos[0].end, 32325184)
        self.assertEqual(pos[0].strand, "+")

    def test_cds_to_exon_neg_strand_by_transcript_id(self):
        pos = cds_to_exon("ENST00000256078", 1)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0].transcript_id, "ENST00000256078")
        self.assertEqual(pos[0].exon_id, "ENSE00000936617")
        self.assertEqual(pos[0].start, 25245274)
        self.assertEqual(pos[0].end, 25245395)
        self.assertEqual(pos[0].strand, "-")

    def test_cds_to_exon_neg_strand_by_transcript_id_2(self):
        pos = cds_to_exon("ENST00000256078", 144)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0].transcript_id, "ENST00000256078")
        self.assertEqual(pos[0].exon_id, "ENSE00001719809")
        self.assertEqual(pos[0].start, 25227234)
        self.assertEqual(pos[0].end, 25227412)
        self.assertEqual(pos[0].strand, "-")

    def test_cds_to_exon_neg_strand_by_transcript_id_3(self):
        pos = cds_to_exon("ENST00000256078", 144, 156)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0].transcript_id, "ENST00000256078")
        self.assertEqual(pos[0].exon_id, "ENSE00001719809")
        self.assertEqual(pos[0].start, 25227234)
        self.assertEqual(pos[0].end, 25227412)
        self.assertEqual(pos[0].strand, "-")

    def test_cds_to_exon_neg_strand_by_transcript_name(self):
        pos = cds_to_exon("KRAS-201", 1)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0].transcript_id, "ENST00000256078")
        self.assertEqual(pos[0].exon_id, "ENSE00000936617")
        self.assertEqual(pos[0].start, 25245274)
        self.assertEqual(pos[0].end, 25245395)
        self.assertEqual(pos[0].strand, "-")

    def test_cds_to_exon_neg_strand_by_transcript_name_2(self):
        pos = cds_to_exon("KRAS-201", 144)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0].transcript_id, "ENST00000256078")
        self.assertEqual(pos[0].exon_id, "ENSE00001719809")
        self.assertEqual(pos[0].start, 25227234)
        self.assertEqual(pos[0].end, 25227412)
        self.assertEqual(pos[0].strand, "-")

    def test_cds_to_exon_neg_strand_by_transcript_name_3(self):
        pos = cds_to_exon("KRAS-201", 144, 156)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0].transcript_id, "ENST00000256078")
        self.assertEqual(pos[0].exon_id, "ENSE00001719809")
        self.assertEqual(pos[0].start, 25227234)
        self.assertEqual(pos[0].end, 25227412)
        self.assertEqual(pos[0].strand, "-")

    def test_cds_to_exon_different_exons(self):
        # positions are on different exons
        self.assertRaises(ValueError, cds_to_exon, "ENST00000380152", 300, 1000)


class TestCdsToGene(unittest.TestCase):
    def test_cds_to_gene_pos_strand_by_transcript_id(self):
        pos = cds_to_gene("ENST00000288135", 1)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0].gene_id, "ENSG00000157404")
        self.assertEqual(pos[0].start, 54658015)
        self.assertEqual(pos[0].end, 54658015)
        self.assertEqual(pos[0].strand, "+")

    def test_cds_to_gene_pos_strand_by_transcript_id_2(self):
        pos = cds_to_gene("ENST00000288135", 213)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0].gene_id, "ENSG00000157404")
        self.assertEqual(pos[0].start, 54695657)
        self.assertEqual(pos[0].end, 54695657)
        self.assertEqual(pos[0].strand, "+")

    def test_cds_to_gene_pos_strand_by_transcript_id_3(self):
        pos = cds_to_gene("ENST00000288135", 213, 221)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0].gene_id, "ENSG00000157404")
        self.assertEqual(pos[0].start, 54695657)
        self.assertEqual(pos[0].end, 54695665)
        self.assertEqual(pos[0].strand, "+")

    def test_cds_to_gene_pos_strand_by_transcript_name(self):
        pos = cds_to_gene("KIT-201", 1)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0].gene_id, "ENSG00000157404")
        self.assertEqual(pos[0].start, 54658015)
        self.assertEqual(pos[0].end, 54658015)
        self.assertEqual(pos[0].strand, "+")

    def test_cds_to_gene_pos_strand_by_transcript_name_2(self):
        pos = cds_to_gene("KIT-201", 213)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0].gene_id, "ENSG00000157404")
        self.assertEqual(pos[0].start, 54695657)
        self.assertEqual(pos[0].end, 54695657)
        self.assertEqual(pos[0].strand, "+")

    def test_cds_to_gene_pos_strand_by_transcript_name_3(self):
        pos = cds_to_gene("KIT-201", 213, 221)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0].gene_id, "ENSG00000157404")
        self.assertEqual(pos[0].start, 54695657)
        self.assertEqual(pos[0].end, 54695665)
        self.assertEqual(pos[0].strand, "+")

    def test_cds_to_gene_neg_strand_by_transcript_id(self):
        pos = cds_to_gene("ENST00000256078", 1)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0].gene_id, "ENSG00000133703")
        self.assertEqual(pos[0].start, 25245384)
        self.assertEqual(pos[0].end, 25245384)
        self.assertEqual(pos[0].strand, "-")

    def test_cds_to_gene_neg_strand_by_transcript_id_2(self):
        pos = cds_to_gene("ENST00000256078", 201)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0].gene_id, "ENSG00000133703")
        self.assertEqual(pos[0].start, 25227323)
        self.assertEqual(pos[0].end, 25227323)
        self.assertEqual(pos[0].strand, "-")

    def test_cds_to_gene_neg_strand_by_transcript_id_3(self):
        pos = cds_to_gene("ENST00000256078", 201, 301)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0].gene_id, "ENSG00000133703")
        self.assertEqual(pos[0].start, 25225763)
        self.assertEqual(pos[0].end, 25227323)
        self.assertEqual(pos[0].strand, "-")

    def test_cds_to_gene_neg_strand_by_transcript_name(self):
        pos = cds_to_gene("KRAS-201", 1)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0].gene_id, "ENSG00000133703")
        self.assertEqual(pos[0].start, 25245384)
        self.assertEqual(pos[0].end, 25245384)
        self.assertEqual(pos[0].strand, "-")

    def test_cds_to_gene_neg_strand_by_transcript_name_2(self):
        pos = cds_to_gene("KRAS-201", 201)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0].gene_id, "ENSG00000133703")
        self.assertEqual(pos[0].start, 25227323)
        self.assertEqual(pos[0].end, 25227323)
        self.assertEqual(pos[0].strand, "-")

    def test_cds_to_gene_neg_strand_by_transcript_name_3(self):
        pos = cds_to_gene("KRAS-201", 201, 301)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0].gene_id, "ENSG00000133703")
        self.assertEqual(pos[0].start, 25225763)
        self.assertEqual(pos[0].end, 25227323)
        self.assertEqual(pos[0].strand, "-")


class TestCdsToProtein(unittest.TestCase):
    def test_cds_to_protein_pos_strand_by_transcript_id(self):
        pos = cds_to_protein("ENST00000288135", 1)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0].protein_id, "ENSP00000288135")
        self.assertEqual(pos[0].start, 1)
        self.assertEqual(pos[0].end, 1)
        self.assertEqual(pos[0].strand, "+")

    def test_cds_to_protein_pos_strand_by_transcript_id_2(self):
        pos = cds_to_protein("ENST00000288135", 213)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0].protein_id, "ENSP00000288135")
        self.assertEqual(pos[0].start, 71)
        self.assertEqual(pos[0].end, 71)
        self.assertEqual(pos[0].strand, "+")

    def test_cds_to_protein_pos_strand_by_transcript_id_3(self):
        pos = cds_to_protein("ENST00000288135", 213, 221)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0].protein_id, "ENSP00000288135")
        self.assertEqual(pos[0].start, 71)
        self.assertEqual(pos[0].end, 74)
        self.assertEqual(pos[0].strand, "+")

    def test_cds_to_protein_pos_strand_by_transcript_name(self):
        pos = cds_to_protein("KIT-201", 1)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0].protein_id, "ENSP00000288135")
        self.assertEqual(pos[0].start, 1)
        self.assertEqual(pos[0].end, 1)
        self.assertEqual(pos[0].strand, "+")

    def test_cds_to_protein_pos_strand_by_transcript_name_2(self):
        pos = cds_to_protein("KIT-201", 213)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0].protein_id, "ENSP00000288135")
        self.assertEqual(pos[0].start, 71)
        self.assertEqual(pos[0].end, 71)
        self.assertEqual(pos[0].strand, "+")

    def test_cds_to_protein_pos_strand_by_transcript_name_3(self):
        pos = cds_to_protein("KIT-201", 213, 221)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0].protein_id, "ENSP00000288135")
        self.assertEqual(pos[0].start, 71)
        self.assertEqual(pos[0].end, 74)
        self.assertEqual(pos[0].strand, "+")

    def test_cds_to_protein_neg_strand_by_transcript_id(self):
        pos = cds_to_protein("ENST00000256078", 1)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0].protein_id, "ENSP00000256078")
        self.assertEqual(pos[0].start, 1)
        self.assertEqual(pos[0].end, 1)
        self.assertEqual(pos[0].strand, "-")

    def test_cds_to_protein_neg_strand_by_transcript_id_2(self):
        pos = cds_to_protein("ENST00000256078", 26)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0].protein_id, "ENSP00000256078")
        self.assertEqual(pos[0].start, 9)
        self.assertEqual(pos[0].end, 9)
        self.assertEqual(pos[0].strand, "-")

    def test_cds_to_protein_neg_strand_by_transcript_id_3(self):
        pos = cds_to_protein("ENST00000256078", 26, 56)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0].protein_id, "ENSP00000256078")
        self.assertEqual(pos[0].start, 9)
        self.assertEqual(pos[0].end, 19)
        self.assertEqual(pos[0].strand, "-")

    def test_cds_to_protein_neg_strand_by_transcript_name(self):
        pos = cds_to_protein("KRAS-201", 1)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0].protein_id, "ENSP00000256078")
        self.assertEqual(pos[0].start, 1)
        self.assertEqual(pos[0].end, 1)
        self.assertEqual(pos[0].strand, "-")

    def test_cds_to_protein_neg_strand_by_transcript_name_2(self):
        pos = cds_to_protein("KRAS-201", 26)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0].protein_id, "ENSP00000256078")
        self.assertEqual(pos[0].start, 9)
        self.assertEqual(pos[0].end, 9)
        self.assertEqual(pos[0].strand, "-")

    def test_cds_to_protein_neg_strand_by_transcript_name_3(self):
        pos = cds_to_protein("KRAS-201", 26, 56)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0].protein_id, "ENSP00000256078")
        self.assertEqual(pos[0].start, 9)
        self.assertEqual(pos[0].end, 19)
        self.assertEqual(pos[0].strand, "-")


class TestCdsToTranscript(unittest.TestCase):
    def test_cds_to_transcript_pos_strand_by_transcript_id(self):
        pos = cds_to_transcript("ENST00000288135", 1)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0].transcript_id, "ENST00000288135")
        self.assertEqual(pos[0].start, 98)
        self.assertEqual(pos[0].end, 98)
        self.assertEqual(pos[0].strand, "+")

    def test_cds_to_transcript_pos_strand_by_transcript_id_2(self):
        pos = cds_to_transcript("ENST00000288135", 213)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0].transcript_id, "ENST00000288135")
        self.assertEqual(pos[0].start, 310)
        self.assertEqual(pos[0].end, 310)
        self.assertEqual(pos[0].strand, "+")

    def test_cds_to_transcript_pos_strand_by_transcript_id_3(self):
        pos = cds_to_transcript("ENST00000288135", 213, 221)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0].transcript_id, "ENST00000288135")
        self.assertEqual(pos[0].start, 310)
        self.assertEqual(pos[0].end, 318)
        self.assertEqual(pos[0].strand, "+")

    def test_cds_to_transcript_pos_strand_by_transcript_name(self):
        pos = cds_to_transcript("KIT-201", 1)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0].transcript_id, "ENST00000288135")
        self.assertEqual(pos[0].start, 98)
        self.assertEqual(pos[0].end, 98)
        self.assertEqual(pos[0].strand, "+")

    def test_cds_to_transcript_pos_strand_by_transcript_name_2(self):
        pos = cds_to_transcript("KIT-201", 213)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0].transcript_id, "ENST00000288135")
        self.assertEqual(pos[0].start, 310)
        self.assertEqual(pos[0].end, 310)
        self.assertEqual(pos[0].strand, "+")

    def test_cds_to_transcript_pos_strand_by_transcript_name_3(self):
        pos = cds_to_transcript("KIT-201", 213, 221)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0].transcript_id, "ENST00000288135")
        self.assertEqual(pos[0].start, 310)
        self.assertEqual(pos[0].end, 318)
        self.assertEqual(pos[0].strand, "+")

    def test_cds_to_transcript_neg_strand_by_transcript_id(self):
        pos = cds_to_transcript("ENST00000256078", 1)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0].transcript_id, "ENST00000256078")
        self.assertEqual(pos[0].start, 191)
        self.assertEqual(pos[0].end, 191)
        self.assertEqual(pos[0].strand, "-")

    def test_cds_to_transcript_neg_strand_by_transcript_id_2(self):
        pos = cds_to_transcript("ENST00000256078", 87)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0].transcript_id, "ENST00000256078")
        self.assertEqual(pos[0].start, 277)
        self.assertEqual(pos[0].end, 277)
        self.assertEqual(pos[0].strand, "-")

    def test_cds_to_transcript_neg_strand_by_transcript_id_3(self):
        pos = cds_to_transcript("ENST00000256078", 87, 92)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0].transcript_id, "ENST00000256078")
        self.assertEqual(pos[0].start, 277)
        self.assertEqual(pos[0].end, 282)
        self.assertEqual(pos[0].strand, "-")

    def test_cds_to_transcript_neg_strand_by_transcript_name(self):
        pos = cds_to_transcript("KRAS-201", 1)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0].transcript_id, "ENST00000256078")
        self.assertEqual(pos[0].start, 191)
        self.assertEqual(pos[0].end, 191)
        self.assertEqual(pos[0].strand, "-")

    def test_cds_to_transcript_neg_strand_by_transcript_name_2(self):
        pos = cds_to_transcript("KRAS-201", 87)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0].transcript_id, "ENST00000256078")
        self.assertEqual(pos[0].start, 277)
        self.assertEqual(pos[0].end, 277)
        self.assertEqual(pos[0].strand, "-")

    def test_cds_to_transcript_neg_strand_by_transcript_name_3(self):
        pos = cds_to_transcript("KRAS-201", 87, 92)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0].transcript_id, "ENST00000256078")
        self.assertEqual(pos[0].start, 277)
        self.assertEqual(pos[0].end, 282)
        self.assertEqual(pos[0].strand, "-")


class TestContigToCds(unittest.TestCase):
    def test_contig_to_cds(self):
        pos = contig_to_cds("5", 1294501)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0].transcript_id, "ENST00000310581")
        self.assertEqual(pos[0].start, 385)
        self.assertEqual(pos[0].end, 385)
        self.assertEqual(pos[0].strand, "-")

    def test_contig_to_cds_2(self):
        pos = contig_to_cds("5", 1294497, 1294501)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0].transcript_id, "ENST00000310581")
        self.assertEqual(pos[0].start, 385)
        self.assertEqual(pos[0].end, 389)
        self.assertEqual(pos[0].strand, "-")


class TestContigToGene(unittest.TestCase):
    def test_contig_to_gene(self):
        pos = contig_to_gene("5", 1253147, 1295068)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0].gene_id, "ENSG00000164362")
        self.assertEqual(pos[0].start, 1253147)
        self.assertEqual(pos[0].end, 1295068)
        self.assertEqual(pos[0].strand, "-")

    def test_contig_to_gene_2(self):
        pos = contig_to_gene("5", 1254000)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0].gene_id, "ENSG00000164362")
        self.assertEqual(pos[0].start, 1254000)
        self.assertEqual(pos[0].end, 1254000)
        self.assertEqual(pos[0].strand, "-")


class TestContigToProtein(unittest.TestCase):
    def test_contig_to_protein(self):
        pos = contig_to_protein("5", 1294501)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0].protein_id, "ENSP00000309572")
        self.assertEqual(pos[0].start, 129)
        self.assertEqual(pos[0].end, 129)
        self.assertEqual(pos[0].strand, "-")

    def test_contig_to_protein_2(self):
        pos = contig_to_protein("5", 1294497, 1294501)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0].protein_id, "ENSP00000309572")
        self.assertEqual(pos[0].start, 129)
        self.assertEqual(pos[0].end, 130)
        self.assertEqual(pos[0].strand, "-")


class TestContigToTranscript(unittest.TestCase):
    def test_contig_to_transcript(self):
        pos = contig_to_transcript("5", 1294501)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0].transcript_id, "ENST00000310581")
        self.assertEqual(pos[0].start, 464)
        self.assertEqual(pos[0].end, 464)
        self.assertEqual(pos[0].strand, "-")

    def test_contig_to_transcript_2(self):
        pos = contig_to_transcript("5", 1294497, 1294501)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0].transcript_id, "ENST00000310581")
        self.assertEqual(pos[0].start, 464)
        self.assertEqual(pos[0].end, 468)
        self.assertEqual(pos[0].strand, "-")


class TestExonToCds(unittest.TestCase):
    def test_exon_to_cds_pos_strand_by_exon_id(self):
        pos = exon_to_cds("ENSE00003659301")
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0].transcript_id, "ENST00000380152")
        self.assertEqual(pos[0].start, 317)
        self.assertEqual(pos[0].end, 425)
        self.assertEqual(pos[0].strand, "+")

    # TODO: https://github.com/mattdoug604/ensembl_map/issues/1
    def test_exon_to_cds_neg_strand_by_exon_id(self):
        pos = exon_to_cds("ENSE00001189807")
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0].transcript_id, "ENST00000256078")
        self.assertEqual(pos[0].start, 451)
        self.assertEqual(pos[0].end, 567)
        self.assertEqual(pos[0].strand, "-")

    def test_exon_to_cds_neg_strand_by_exon_id_2(self):
        pos = exon_to_cds("ENSE00000936617")
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0].transcript_id, "ENST00000256078")
        self.assertEqual(pos[0].start, 1)
        self.assertEqual(pos[0].end, 111)
        self.assertEqual(pos[0].strand, "-")

    def test_exon_to_cds_neg_strand_by_exon_id_3(self):
        # 'ENSE00001189804' is not part of the CDS
        self.assertRaises(ValueError, exon_to_cds, "ENSE00001189804")


class TestExonToGene(unittest.TestCase):
    def test_exon_to_gene_pos_strand_by_exon_id(self):
        pos = exon_to_gene("ENSE00003659301")
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0].gene_id, "ENSG00000139618")
        self.assertEqual(pos[0].start, 32325076)
        self.assertEqual(pos[0].end, 32325184)
        self.assertEqual(pos[0].strand, "+")

    def test_exon_to_gene_neg_strand_by_exon_id(self):
        pos = exon_to_gene("ENSE00001189807")
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0].gene_id, "ENSG00000133703")
        self.assertEqual(pos[0].start, 25215437)
        self.assertEqual(pos[0].end, 25215560)
        self.assertEqual(pos[0].strand, "-")


class TestExonToProtein(unittest.TestCase):
    def test_exon_to_protein_pos_strand_by_exon_id(self):
        pos = exon_to_protein("ENSE00003659301")
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0].protein_id, "ENSP00000369497")
        self.assertEqual(pos[0].start, 106)
        self.assertEqual(pos[0].end, 142)
        self.assertEqual(pos[0].strand, "+")

    def test_exon_to_protein_neg_strand_by_exon_id(self):
        pos = exon_to_protein("ENSE00001719809")
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0].protein_id, "ENSP00000256078")
        self.assertEqual(pos[0].start, 38)
        self.assertEqual(pos[0].end, 97)
        self.assertEqual(pos[0].strand, "-")


class TestExonToTranscript(unittest.TestCase):
    def test_exon_to_transcript_pos_strand_by_exon_id(self):
        pos = exon_to_transcript("ENSE00003659301")
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0].transcript_id, "ENST00000380152")
        self.assertEqual(pos[0].start, 550)
        self.assertEqual(pos[0].end, 658)
        self.assertEqual(pos[0].strand, "+")

    def test_exon_to_transcript_neg_strand_by_exon_id(self):
        pos = exon_to_transcript("ENSE00001719809")
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0].transcript_id, "ENST00000256078")
        self.assertEqual(pos[0].start, 302)
        self.assertEqual(pos[0].end, 480)
        self.assertEqual(pos[0].strand, "-")


class TestGeneToCds(unittest.TestCase):
    def test_gene_to_cds_pos_strand_by_gene_id(self):
        pos = gene_to_cds("ENSG00000157404", 54658015)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0].transcript_id, "ENST00000288135")
        self.assertEqual(pos[0].start, 1)
        self.assertEqual(pos[0].end, 1)
        self.assertEqual(pos[0].strand, "+")

    def test_gene_to_cds_pos_strand_by_gene_id_2(self):
        pos = gene_to_cds("ENSG00000157404", 54695656)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0].transcript_id, "ENST00000288135")
        self.assertEqual(pos[0].start, 212)
        self.assertEqual(pos[0].end, 212)
        self.assertEqual(pos[0].strand, "+")

    def test_gene_to_cds_pos_strand_by_gene_id_3(self):
        pos = gene_to_cds("ENSG00000157404", 54695656, 54695664)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0].transcript_id, "ENST00000288135")
        self.assertEqual(pos[0].start, 212)
        self.assertEqual(pos[0].end, 220)
        self.assertEqual(pos[0].strand, "+")

    def test_gene_to_cds_pos_strand_by_gene_name(self):
        pos = gene_to_cds("KIT", 54658015)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0].transcript_id, "ENST00000288135")
        self.assertEqual(pos[0].start, 1)
        self.assertEqual(pos[0].end, 1)
        self.assertEqual(pos[0].strand, "+")

    def test_gene_to_cds_pos_strand_by_gene_name_2(self):
        pos = gene_to_cds("KIT", 54695656)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0].transcript_id, "ENST00000288135")
        self.assertEqual(pos[0].start, 212)
        self.assertEqual(pos[0].end, 212)
        self.assertEqual(pos[0].strand, "+")

    def test_gene_to_cds_pos_strand_by_gene_name_3(self):
        pos = gene_to_cds("KIT", 54695656, 54695664)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0].transcript_id, "ENST00000288135")
        self.assertEqual(pos[0].start, 212)
        self.assertEqual(pos[0].end, 220)
        self.assertEqual(pos[0].strand, "+")

    def test_gene_to_cds_neg_strand_by_gene_id(self):
        pos = gene_to_cds("ENSG00000133703", 25245384)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0].transcript_id, "ENST00000256078")
        self.assertEqual(pos[0].start, 1)
        self.assertEqual(pos[0].end, 1)
        self.assertEqual(pos[0].strand, "-")

    def test_gene_to_cds_neg_strand_by_gene_id_2(self):
        pos = gene_to_cds("ENSG00000133703", 25225740)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0].transcript_id, "ENST00000256078")
        self.assertEqual(pos[0].start, 324)
        self.assertEqual(pos[0].end, 324)
        self.assertEqual(pos[0].strand, "-")

    def test_gene_to_cds_neg_strand_by_gene_id_3(self):
        pos = gene_to_cds("ENSG00000133703", 25225740, 25225750)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0].transcript_id, "ENST00000256078")
        self.assertEqual(pos[0].start, 314)
        self.assertEqual(pos[0].end, 324)
        self.assertEqual(pos[0].strand, "-")

    def test_gene_to_cds_neg_strand_by_gene_name(self):
        pos = gene_to_cds("KRAS", 25245384)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0].transcript_id, "ENST00000256078")
        self.assertEqual(pos[0].start, 1)
        self.assertEqual(pos[0].end, 1)
        self.assertEqual(pos[0].strand, "-")

    def test_gene_to_cds_neg_strand_by_gene_name_2(self):
        pos = gene_to_cds("KRAS", 25225740)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0].transcript_id, "ENST00000256078")
        self.assertEqual(pos[0].start, 324)
        self.assertEqual(pos[0].end, 324)
        self.assertEqual(pos[0].strand, "-")

    def test_gene_to_cds_neg_strand_by_gene_name_3(self):
        pos = gene_to_cds("KRAS", 25225740, 25225750)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0].transcript_id, "ENST00000256078")
        self.assertEqual(pos[0].start, 314)
        self.assertEqual(pos[0].end, 324)
        self.assertEqual(pos[0].strand, "-")


class TestGeneToExon(unittest.TestCase):
    def test_gene_to_exon_pos_strand_by_gene_id(self):
        pos = gene_to_exon("ENSG00000157404", 54657918)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0].transcript_id, "ENST00000288135")
        self.assertEqual(pos[0].exon_id, "ENSE00001905199")
        self.assertEqual(pos[0].start, 54657918)
        self.assertEqual(pos[0].end, 54658081)
        self.assertEqual(pos[0].strand, "+")

    def test_gene_to_exon_pos_strand_by_gene_id_2(self):
        pos = gene_to_exon("ENSG00000157404", 54698301)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0].transcript_id, "ENST00000288135")
        self.assertEqual(pos[0].exon_id, "ENSE00001074448")
        self.assertEqual(pos[0].start, 54698284)
        self.assertEqual(pos[0].end, 54698565)
        self.assertEqual(pos[0].strand, "+")

    def test_gene_to_exon_pos_strand_by_gene_id_3(self):
        pos = gene_to_exon("ENSG00000157404", 54698301, 54698330)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0].transcript_id, "ENST00000288135")
        self.assertEqual(pos[0].exon_id, "ENSE00001074448")
        self.assertEqual(pos[0].start, 54698284)
        self.assertEqual(pos[0].end, 54698565)
        self.assertEqual(pos[0].strand, "+")

    def test_gene_to_exon_pos_strand_by_gene_name(self):
        pos = gene_to_exon("KIT", 54657918)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0].transcript_id, "ENST00000288135")
        self.assertEqual(pos[0].exon_id, "ENSE00001905199")
        self.assertEqual(pos[0].start, 54657918)
        self.assertEqual(pos[0].end, 54658081)
        self.assertEqual(pos[0].strand, "+")

    def test_gene_to_exon_pos_strand_by_gene_name_2(self):
        pos = gene_to_exon("KIT", 54698301)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0].transcript_id, "ENST00000288135")
        self.assertEqual(pos[0].exon_id, "ENSE00001074448")
        self.assertEqual(pos[0].start, 54698284)
        self.assertEqual(pos[0].end, 54698565)
        self.assertEqual(pos[0].strand, "+")

    def test_gene_to_exon_pos_strand_by_gene_name_3(self):
        pos = gene_to_exon("KIT", 54698301, 54698330)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0].transcript_id, "ENST00000288135")
        self.assertEqual(pos[0].exon_id, "ENSE00001074448")
        self.assertEqual(pos[0].start, 54698284)
        self.assertEqual(pos[0].end, 54698565)
        self.assertEqual(pos[0].strand, "+")

    def test_gene_to_exon_neg_strand_by_gene_id(self):
        pos = gene_to_exon("ENSG00000133703", 25250929)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0].transcript_id, "ENST00000256078")
        self.assertEqual(pos[0].exon_id, "ENSE00001189804")
        self.assertEqual(pos[0].start, 25250751)
        self.assertEqual(pos[0].end, 25250929)
        self.assertEqual(pos[0].strand, "-")

    def test_gene_to_exon_neg_strand_by_gene_id_2(self):
        pos = gene_to_exon("ENSG00000133703", 25227400)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0].transcript_id, "ENST00000256078")
        self.assertEqual(pos[0].exon_id, "ENSE00001719809")
        self.assertEqual(pos[0].start, 25227234)
        self.assertEqual(pos[0].end, 25227412)
        self.assertEqual(pos[0].strand, "-")

    def test_gene_to_exon_neg_strand_by_gene_id_3(self):
        pos = gene_to_exon("ENSG00000133703", 25227380, 25227400)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0].transcript_id, "ENST00000256078")
        self.assertEqual(pos[0].exon_id, "ENSE00001719809")
        self.assertEqual(pos[0].start, 25227234)
        self.assertEqual(pos[0].end, 25227412)
        self.assertEqual(pos[0].strand, "-")

    def test_gene_to_exon_neg_strand_by_gene_name(self):
        pos = gene_to_exon("KRAS", 25250929)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0].transcript_id, "ENST00000256078")
        self.assertEqual(pos[0].exon_id, "ENSE00001189804")
        self.assertEqual(pos[0].start, 25250751)
        self.assertEqual(pos[0].end, 25250929)
        self.assertEqual(pos[0].strand, "-")

    def test_gene_to_exon_neg_strand_by_gene_name_2(self):
        pos = gene_to_exon("KRAS", 25227400)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0].transcript_id, "ENST00000256078")
        self.assertEqual(pos[0].exon_id, "ENSE00001719809")
        self.assertEqual(pos[0].start, 25227234)
        self.assertEqual(pos[0].end, 25227412)
        self.assertEqual(pos[0].strand, "-")

    def test_gene_to_exon_neg_strand_by_gene_name_3(self):
        pos = gene_to_exon("KRAS", 25227380, 25227400)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0].transcript_id, "ENST00000256078")
        self.assertEqual(pos[0].exon_id, "ENSE00001719809")
        self.assertEqual(pos[0].start, 25227234)
        self.assertEqual(pos[0].end, 25227412)
        self.assertEqual(pos[0].strand, "-")


class TestGeneToProtein(unittest.TestCase):
    def test_gene_to_protein_pos_strand_by_gene_id(self):
        pos = gene_to_protein("ENSG00000157404", 54658015)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0].protein_id, "ENSP00000288135")
        self.assertEqual(pos[0].start, 1)
        self.assertEqual(pos[0].end, 1)
        self.assertEqual(pos[0].strand, "+")

    def test_gene_to_protein_pos_strand_by_gene_id_2(self):
        pos = gene_to_protein("ENSG00000157404", 54695656)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0].protein_id, "ENSP00000288135")
        self.assertEqual(pos[0].start, 71)
        self.assertEqual(pos[0].end, 71)
        self.assertEqual(pos[0].strand, "+")

    def test_gene_to_protein_pos_strand_by_gene_id_3(self):
        pos = gene_to_protein("ENSG00000157404", 54695656, 54695664)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0].protein_id, "ENSP00000288135")
        self.assertEqual(pos[0].start, 71)
        self.assertEqual(pos[0].end, 74)
        self.assertEqual(pos[0].strand, "+")

    def test_gene_to_protein_pos_strand_by_gene_name(self):
        pos = gene_to_protein("KIT", 54658015)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0].protein_id, "ENSP00000288135")
        self.assertEqual(pos[0].start, 1)
        self.assertEqual(pos[0].end, 1)
        self.assertEqual(pos[0].strand, "+")

    def test_gene_to_protein_pos_strand_by_gene_name_2(self):
        pos = gene_to_protein("KIT", 54695656)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0].protein_id, "ENSP00000288135")
        self.assertEqual(pos[0].start, 71)
        self.assertEqual(pos[0].end, 71)
        self.assertEqual(pos[0].strand, "+")

    def test_gene_to_protein_pos_strand_by_gene_name_3(self):
        pos = gene_to_protein("KIT", 54695656, 54695664)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0].protein_id, "ENSP00000288135")
        self.assertEqual(pos[0].start, 71)
        self.assertEqual(pos[0].end, 74)
        self.assertEqual(pos[0].strand, "+")

    def test_gene_to_protein_neg_strand_by_gene_id(self):
        pos = gene_to_protein("ENSG00000133703", 25245384)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0].protein_id, "ENSP00000256078")
        self.assertEqual(pos[0].start, 1)
        self.assertEqual(pos[0].end, 1)
        self.assertEqual(pos[0].strand, "-")

    def test_gene_to_protein_neg_strand_by_gene_id_2(self):
        pos = gene_to_protein("ENSG00000133703", 25227402)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0].protein_id, "ENSP00000256078")
        self.assertEqual(pos[0].start, 41)
        self.assertEqual(pos[0].end, 41)
        self.assertEqual(pos[0].strand, "-")

    def test_gene_to_protein_neg_strand_by_gene_id_3(self):
        pos = gene_to_protein("ENSG00000133703", 25227382, 25227402)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0].protein_id, "ENSP00000256078")
        self.assertEqual(pos[0].start, 41)
        self.assertEqual(pos[0].end, 48)
        self.assertEqual(pos[0].strand, "-")

    def test_gene_to_protein_neg_strand_by_gene_name(self):
        pos = gene_to_protein("KRAS", 25245384)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0].protein_id, "ENSP00000256078")
        self.assertEqual(pos[0].start, 1)
        self.assertEqual(pos[0].end, 1)
        self.assertEqual(pos[0].strand, "-")

    def test_gene_to_protein_neg_strand_by_gene_name_2(self):
        pos = gene_to_protein("KRAS", 25227402)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0].protein_id, "ENSP00000256078")
        self.assertEqual(pos[0].start, 41)
        self.assertEqual(pos[0].end, 41)
        self.assertEqual(pos[0].strand, "-")

    def test_gene_to_protein_neg_strand_by_gene_name_3(self):
        pos = gene_to_protein("KRAS", 25227382, 25227402)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0].protein_id, "ENSP00000256078")
        self.assertEqual(pos[0].start, 41)
        self.assertEqual(pos[0].end, 48)
        self.assertEqual(pos[0].strand, "-")


class TestGeneToTranscript(unittest.TestCase):
    def test_gene_to_transcript_pos_strand_by_gene_id(self):
        pos = gene_to_transcript("ENSG00000188554", 43170481)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0].transcript_id, "ENST00000341165")
        self.assertEqual(pos[0].start, 1)
        self.assertEqual(pos[0].end, 1)
        self.assertEqual(pos[0].strand, "+")

    def test_gene_to_transcript_pos_strand_by_gene_id_2(self):
        pos = gene_to_transcript("ENSG00000188554", 43189698)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0].transcript_id, "ENST00000341165")
        self.assertEqual(pos[0].start, 731)
        self.assertEqual(pos[0].end, 731)
        self.assertEqual(pos[0].strand, "+")

    def test_gene_to_transcript_pos_strand_by_gene_id_3(self):
        pos = gene_to_transcript("ENSG00000188554", 43189698, 43189784)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0].transcript_id, "ENST00000341165")
        self.assertEqual(pos[0].start, 731)
        self.assertEqual(pos[0].end, 817)
        self.assertEqual(pos[0].strand, "+")

    def test_gene_to_transcript_pos_strand_by_gene_name(self):
        pos = gene_to_transcript("NBR1", 43170481)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0].transcript_id, "ENST00000341165")
        self.assertEqual(pos[0].start, 1)
        self.assertEqual(pos[0].end, 1)
        self.assertEqual(pos[0].strand, "+")

    def test_gene_to_transcript_pos_strand_by_gene_name_2(self):
        pos = gene_to_transcript("NBR1", 43189698)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0].transcript_id, "ENST00000341165")
        self.assertEqual(pos[0].start, 731)
        self.assertEqual(pos[0].end, 731)
        self.assertEqual(pos[0].strand, "+")

    def test_gene_to_transcript_pos_strand_by_gene_name_3(self):
        pos = gene_to_transcript("NBR1", 43189698, 43189784)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0].transcript_id, "ENST00000341165")
        self.assertEqual(pos[0].start, 731)
        self.assertEqual(pos[0].end, 817)
        self.assertEqual(pos[0].strand, "+")

    def test_gene_to_transcript_neg_strand_by_gene_id(self):
        pos = gene_to_transcript("ENSG00000133703", 25250929)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0].transcript_id, "ENST00000256078")
        self.assertEqual(pos[0].start, 1)
        self.assertEqual(pos[0].end, 1)
        self.assertEqual(pos[0].strand, "-")

    def test_gene_to_transcript_neg_strand_by_gene_id_2(self):
        pos = gene_to_transcript("ENSG00000133703", 25225634)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0].transcript_id, "ENST00000256078")
        self.assertEqual(pos[0].start, 620)
        self.assertEqual(pos[0].end, 620)
        self.assertEqual(pos[0].strand, "-")

    def test_gene_to_transcript_neg_strand_by_gene_id_3(self):
        pos = gene_to_transcript("ENSG00000133703", 25225634, 25225651)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0].transcript_id, "ENST00000256078")
        self.assertEqual(pos[0].start, 603)
        self.assertEqual(pos[0].end, 620)
        self.assertEqual(pos[0].strand, "-")

    def test_gene_to_transcript_neg_strand_by_gene_name(self):
        pos = gene_to_transcript("KRAS", 25225634, 25225651)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0].transcript_id, "ENST00000256078")
        self.assertEqual(pos[0].start, 603)
        self.assertEqual(pos[0].end, 620)
        self.assertEqual(pos[0].strand, "-")

    def test_gene_to_transcript_neg_strand_by_gene_name_2(self):
        pos = gene_to_transcript("KRAS", 25225634, 25225651)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0].transcript_id, "ENST00000256078")
        self.assertEqual(pos[0].start, 603)
        self.assertEqual(pos[0].end, 620)
        self.assertEqual(pos[0].strand, "-")

    def test_gene_to_transcript_neg_strand_by_gene_name_3(self):
        pos = gene_to_transcript("KRAS", 25225634, 25225651)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0].transcript_id, "ENST00000256078")
        self.assertEqual(pos[0].start, 603)
        self.assertEqual(pos[0].end, 620)
        self.assertEqual(pos[0].strand, "-")


class TestProteinToCds(unittest.TestCase):
    def test_protein_to_cds_pos_strand_by_protein_id(self):
        pos = protein_to_cds("ENSP00000288135", 1)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0].transcript_id, "ENST00000288135")
        self.assertEqual(pos[0].start, 1)
        self.assertEqual(pos[0].end, 3)
        self.assertEqual(pos[0].strand, "+")

    def test_protein_to_cds_pos_strand_by_protein_id_2(self):
        pos = protein_to_cds("ENSP00000288135", 12)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0].transcript_id, "ENST00000288135")
        self.assertEqual(pos[0].start, 34)
        self.assertEqual(pos[0].end, 36)
        self.assertEqual(pos[0].strand, "+")

    def test_protein_to_cds_pos_strand_by_protein_id_3(self):
        pos = protein_to_cds("ENSP00000288135", 12, 15)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0].transcript_id, "ENST00000288135")
        self.assertEqual(pos[0].start, 34)
        self.assertEqual(pos[0].end, 45)
        self.assertEqual(pos[0].strand, "+")

    def test_protein_to_cds_neg_strand_by_protein_id(self):
        pos = protein_to_cds("ENSP00000308495", 1)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0].transcript_id, "ENST00000311936")
        self.assertEqual(pos[0].start, 1)
        self.assertEqual(pos[0].end, 3)
        self.assertEqual(pos[0].strand, "-")

    def test_protein_to_cds_neg_strand_by_protein_id_2(self):
        pos = protein_to_cds("ENSP00000308495", 123)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0].transcript_id, "ENST00000311936")
        self.assertEqual(pos[0].start, 367)
        self.assertEqual(pos[0].end, 369)
        self.assertEqual(pos[0].strand, "-")

    def test_protein_to_cds_neg_strand_by_protein_id_3(self):
        pos = protein_to_cds("ENSP00000308495", 123, 124)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0].transcript_id, "ENST00000311936")
        self.assertEqual(pos[0].start, 367)
        self.assertEqual(pos[0].end, 372)
        self.assertEqual(pos[0].strand, "-")


class TestProteinToExon(unittest.TestCase):
    def test_protein_to_exon_pos_strand_by_protein_id(self):
        pos = protein_to_exon("ENSP00000288135", 1)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0].transcript_id, "ENST00000288135")
        self.assertEqual(pos[0].exon_id, "ENSE00001905199")
        self.assertEqual(pos[0].start, 54657918)
        self.assertEqual(pos[0].end, 54658081)
        self.assertEqual(pos[0].strand, "+")

    def test_protein_to_exon_pos_strand_by_protein_id_2(self):
        pos = protein_to_exon("ENSP00000288135", 4)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0].transcript_id, "ENST00000288135")
        self.assertEqual(pos[0].exon_id, "ENSE00001905199")
        self.assertEqual(pos[0].start, 54657918)
        self.assertEqual(pos[0].end, 54658081)
        self.assertEqual(pos[0].strand, "+")

    def test_protein_to_exon_pos_strand_by_protein_id_3(self):
        pos = protein_to_exon("ENSP00000288135", 4, 6)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0].transcript_id, "ENST00000288135")
        self.assertEqual(pos[0].exon_id, "ENSE00001905199")
        self.assertEqual(pos[0].start, 54657918)
        self.assertEqual(pos[0].end, 54658081)
        self.assertEqual(pos[0].strand, "+")

    def test_protein_to_exon_neg_strand_by_protein_id(self):
        pos = protein_to_exon("ENSP00000308495", 1)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0].transcript_id, "ENST00000256078")
        self.assertEqual(pos[0].exon_id, "ENSE00000936617")
        self.assertEqual(pos[0].start, 25245274)
        self.assertEqual(pos[0].end, 25245395)
        self.assertEqual(pos[0].strand, "-")

    def test_protein_to_exon_neg_strand_by_protein_id_2(self):
        pos = protein_to_exon("ENSP00000308495", 4)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0].transcript_id, "ENST00000256078")
        self.assertEqual(pos[0].exon_id, "ENSE00000936617")
        self.assertEqual(pos[0].start, 25245274)
        self.assertEqual(pos[0].end, 25245395)
        self.assertEqual(pos[0].strand, "-")

    def test_protein_to_exon_neg_strand_by_protein_id_3(self):
        pos = protein_to_exon("ENSP00000308495", 4, 6)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0].transcript_id, "ENST00000256078")
        self.assertEqual(pos[0].exon_id, "ENSE00000936617")
        self.assertEqual(pos[0].start, 25245274)
        self.assertEqual(pos[0].end, 25245395)
        self.assertEqual(pos[0].strand, "-")


class TestProteinToGene(unittest.TestCase):
    def test_protein_to_gene_pos_strand_by_protein_id(self):
        pos = protein_to_gene("ENSP00000288135", 1)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0].gene_id, "ENSG00000157404")
        self.assertEqual(pos[0].start, 54658015)
        self.assertEqual(pos[0].end, 54658017)
        self.assertEqual(pos[0].strand, "+")

    def test_protein_to_gene_pos_strand_by_protein_id_2(self):
        pos = protein_to_gene("ENSP00000288135", 71)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0].gene_id, "ENSG00000157404")
        self.assertEqual(pos[0].start, 54695655)
        self.assertEqual(pos[0].end, 54695657)
        self.assertEqual(pos[0].strand, "+")

    def test_protein_to_gene_pos_strand_by_protein_id_3(self):
        pos = protein_to_gene("ENSP00000288135", 71, 74)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0].gene_id, "ENSG00000157404")
        self.assertEqual(pos[0].start, 54695655)
        self.assertEqual(pos[0].end, 54695666)
        self.assertEqual(pos[0].strand, "+")

    def test_protein_to_gene_neg_strand_by_protein_id(self):
        pos = protein_to_gene("ENSP00000256078", 1)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0].gene_id, "ENSG00000133703")
        self.assertEqual(pos[0].start, 25245382)
        self.assertEqual(pos[0].end, 25245384)
        self.assertEqual(pos[0].strand, "-")

    def test_protein_to_gene_neg_strand_by_protein_id_2(self):
        pos = protein_to_gene("ENSP00000256078", 6)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0].gene_id, "ENSG00000133703")
        self.assertEqual(pos[0].start, 25245367)
        self.assertEqual(pos[0].end, 25245369)
        self.assertEqual(pos[0].strand, "-")

    def test_protein_to_gene_neg_strand_by_protein_id_3(self):
        pos = protein_to_gene("ENSP00000256078", 6, 9)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0].gene_id, "ENSG00000133703")
        self.assertEqual(pos[0].start, 25245358)
        self.assertEqual(pos[0].end, 25245369)
        self.assertEqual(pos[0].strand, "-")


class TestProteinToTranscript(unittest.TestCase):
    def test_protein_to_transcript_pos_strand_by_protein_id(self):
        pos = protein_to_transcript("ENSP00000288135", 1)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0].transcript_id, "ENST00000288135")
        self.assertEqual(pos[0].start, 98)
        self.assertEqual(pos[0].end, 100)
        self.assertEqual(pos[0].strand, "+")

    def test_protein_to_transcript_pos_strand_by_protein_id_2(self):
        pos = protein_to_transcript("ENSP00000288135", 12)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0].transcript_id, "ENST00000288135")
        self.assertEqual(pos[0].start, 131)
        self.assertEqual(pos[0].end, 133)
        self.assertEqual(pos[0].strand, "+")

    def test_protein_to_transcript_pos_strand_by_protein_id_3(self):
        pos = protein_to_transcript("ENSP00000288135", 12, 21)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0].transcript_id, "ENST00000288135")
        self.assertEqual(pos[0].start, 131)
        self.assertEqual(pos[0].end, 160)
        self.assertEqual(pos[0].strand, "+")

    def test_protein_to_transcript_neg_strand_by_protein_id(self):
        pos = protein_to_transcript("ENSP00000260947", 1)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0].transcript_id, "ENST00000260947")
        self.assertEqual(pos[0].start, 115)
        self.assertEqual(pos[0].end, 117)
        self.assertEqual(pos[0].strand, "-")

    def test_protein_to_transcript_neg_strand_by_protein_id_2(self):
        pos = protein_to_transcript("ENSP00000260947", 71)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0].transcript_id, "ENST00000260947")
        self.assertEqual(pos[0].start, 325)
        self.assertEqual(pos[0].end, 327)
        self.assertEqual(pos[0].strand, "-")

    def test_protein_to_transcript_neg_strand_by_protein_id_3(self):
        pos = protein_to_transcript("ENSP00000260947", 71, 74)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0].transcript_id, "ENST00000260947")
        self.assertEqual(pos[0].start, 325)
        self.assertEqual(pos[0].end, 336)
        self.assertEqual(pos[0].strand, "-")


class TestTranscriptToCds(unittest.TestCase):
    def test_transcript_to_cds_pos_strand_by_transcript_id(self):
        pos = transcript_to_cds("ENST00000288135", 98)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0].transcript_id, "ENST00000288135")
        self.assertEqual(pos[0].start, 1)
        self.assertEqual(pos[0].end, 1)
        self.assertEqual(pos[0].strand, "+")

    def test_transcript_to_cds_pos_strand_by_transcript_id_2(self):
        pos = transcript_to_cds("ENST00000288135", 100)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0].transcript_id, "ENST00000288135")
        self.assertEqual(pos[0].start, 3)
        self.assertEqual(pos[0].end, 3)
        self.assertEqual(pos[0].strand, "+")

    def test_transcript_to_cds_pos_strand_by_transcript_id_3(self):
        pos = transcript_to_cds("ENST00000288135", 98, 100)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0].transcript_id, "ENST00000288135")
        self.assertEqual(pos[0].start, 1)
        self.assertEqual(pos[0].end, 3)
        self.assertEqual(pos[0].strand, "+")

    def test_transcript_to_cds_pos_strand_by_transcript_name(self):
        pos = transcript_to_cds("KIT-201", 98)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0].transcript_id, "ENST00000288135")
        self.assertEqual(pos[0].start, 1)
        self.assertEqual(pos[0].end, 1)
        self.assertEqual(pos[0].strand, "+")

    def test_transcript_to_cds_pos_strand_by_transcript_name_2(self):
        pos = transcript_to_cds("KIT-201", 100)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0].transcript_id, "ENST00000288135")
        self.assertEqual(pos[0].start, 3)
        self.assertEqual(pos[0].end, 3)
        self.assertEqual(pos[0].strand, "+")

    def test_transcript_to_cds_pos_strand_by_transcript_name_3(self):
        pos = transcript_to_cds("KIT-201", 98, 100)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0].transcript_id, "ENST00000288135")
        self.assertEqual(pos[0].start, 1)
        self.assertEqual(pos[0].end, 3)
        self.assertEqual(pos[0].strand, "+")

    def test_transcript_to_cds_neg_strand_by_transcript_id(self):
        pos = transcript_to_cds("ENST00000256078", 191)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0].transcript_id, "ENST00000256078")
        self.assertEqual(pos[0].start, 1)
        self.assertEqual(pos[0].end, 1)
        self.assertEqual(pos[0].strand, "-")

    def test_transcript_to_cds_neg_strand_by_transcript_id_2(self):
        pos = transcript_to_cds("ENST00000256078", 301)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0].transcript_id, "ENST00000256078")
        self.assertEqual(pos[0].start, 111)
        self.assertEqual(pos[0].end, 111)
        self.assertEqual(pos[0].strand, "-")

    def test_transcript_to_cds_neg_strand_by_transcript_id_3(self):
        pos = transcript_to_cds("ENST00000256078", 301, 323)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0].transcript_id, "ENST00000256078")
        self.assertEqual(pos[0].start, 111)
        self.assertEqual(pos[0].end, 133)
        self.assertEqual(pos[0].strand, "-")

    def test_transcript_to_cds_neg_strand_by_transcript_name(self):
        pos = transcript_to_cds("KRAS-201", 191)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0].transcript_id, "ENST00000256078")
        self.assertEqual(pos[0].start, 1)
        self.assertEqual(pos[0].end, 1)
        self.assertEqual(pos[0].strand, "-")

    def test_transcript_to_cds_neg_strand_by_transcript_name_2(self):
        pos = transcript_to_cds("KRAS-201", 301)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0].transcript_id, "ENST00000256078")
        self.assertEqual(pos[0].start, 111)
        self.assertEqual(pos[0].end, 111)
        self.assertEqual(pos[0].strand, "-")

    def test_transcript_to_cds_neg_strand_by_transcript_name_3(self):
        pos = transcript_to_cds("KRAS-201", 301, 323)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0].transcript_id, "ENST00000256078")
        self.assertEqual(pos[0].start, 111)
        self.assertEqual(pos[0].end, 133)
        self.assertEqual(pos[0].strand, "-")


class TestTranscriptToExon(unittest.TestCase):
    def test_transcript_to_exon_pos_strand_by_transcript_id(self):
        pos = transcript_to_exon("ENST00000380152", 1)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0].transcript_id, "ENST00000380152")
        self.assertEqual(pos[0].exon_id, "ENSE00001184784")
        self.assertEqual(pos[0].start, 32315474)
        self.assertEqual(pos[0].end, 32315667)
        self.assertEqual(pos[0].strand, "+")

    def test_transcript_to_exon_pos_strand_by_transcript_id_2(self):
        pos = transcript_to_exon("ENST00000380152", 500)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0].transcript_id, "ENST00000380152")
        self.assertEqual(pos[0].exon_id, "ENSE00003666217")
        self.assertEqual(pos[0].start, 32319077)
        self.assertEqual(pos[0].end, 32319325)
        self.assertEqual(pos[0].strand, "+")

    def test_transcript_to_exon_pos_strand_by_transcript_id_3(self):
        pos = transcript_to_exon("ENST00000380152", 500, 510)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0].transcript_id, "ENST00000380152")
        self.assertEqual(pos[0].exon_id, "ENSE00003666217")
        self.assertEqual(pos[0].start, 32319077)
        self.assertEqual(pos[0].end, 32319325)
        self.assertEqual(pos[0].strand, "+")

    def test_transcript_to_exon_pos_strand_by_transcript_name(self):
        pos = transcript_to_exon("BRCA2-201", 1)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0].transcript_id, "ENST00000380152")
        self.assertEqual(pos[0].exon_id, "ENSE00001184784")
        self.assertEqual(pos[0].start, 32315474)
        self.assertEqual(pos[0].end, 32315667)
        self.assertEqual(pos[0].strand, "+")

    def test_transcript_to_exon_pos_strand_by_transcript_name_2(self):
        pos = transcript_to_exon("BRCA2-201", 500)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0].transcript_id, "ENST00000380152")
        self.assertEqual(pos[0].exon_id, "ENSE00003666217")
        self.assertEqual(pos[0].start, 32319077)
        self.assertEqual(pos[0].end, 32319325)
        self.assertEqual(pos[0].strand, "+")

    def test_transcript_to_exon_pos_strand_by_transcript_name_3(self):
        pos = transcript_to_exon("BRCA2-201", 500, 510)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0].transcript_id, "ENST00000380152")
        self.assertEqual(pos[0].exon_id, "ENSE00003666217")
        self.assertEqual(pos[0].start, 32319077)
        self.assertEqual(pos[0].end, 32319325)
        self.assertEqual(pos[0].strand, "+")

    def test_transcript_to_exon_neg_strand_by_transcript_id(self):
        pos = transcript_to_exon("ENST00000256078", 1)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0].transcript_id, "ENST00000256078")
        self.assertEqual(pos[0].exon_id, "ENSE00001189804")
        self.assertEqual(pos[0].start, 25250751)
        self.assertEqual(pos[0].end, 25250929)
        self.assertEqual(pos[0].strand, "-")

    def test_transcript_to_exon_neg_strand_by_transcript_id_2(self):
        pos = transcript_to_exon("ENST00000256078", 400)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0].transcript_id, "ENST00000256078")
        self.assertEqual(pos[0].exon_id, "ENSE00001719809")
        self.assertEqual(pos[0].start, 25227234)
        self.assertEqual(pos[0].end, 25227412)
        self.assertEqual(pos[0].strand, "-")

    def test_transcript_to_exon_neg_strand_by_transcript_id_3(self):
        pos = transcript_to_exon("ENST00000256078", 400, 420)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0].transcript_id, "ENST00000256078")
        self.assertEqual(pos[0].exon_id, "ENSE00001719809")
        self.assertEqual(pos[0].start, 25227234)
        self.assertEqual(pos[0].end, 25227412)
        self.assertEqual(pos[0].strand, "-")

    def test_transcript_to_exon_neg_strand_by_transcript_name(self):
        pos = transcript_to_exon("KRAS-201", 1)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0].transcript_id, "ENST00000256078")
        self.assertEqual(pos[0].exon_id, "ENSE00001189804")
        self.assertEqual(pos[0].start, 25250751)
        self.assertEqual(pos[0].end, 25250929)
        self.assertEqual(pos[0].strand, "-")

    def test_transcript_to_exon_neg_strand_by_transcript_name_2(self):
        pos = transcript_to_exon("KRAS-201", 400)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0].transcript_id, "ENST00000256078")
        self.assertEqual(pos[0].exon_id, "ENSE00001719809")
        self.assertEqual(pos[0].start, 25227234)
        self.assertEqual(pos[0].end, 25227412)
        self.assertEqual(pos[0].strand, "-")

    def test_transcript_to_exon_neg_strand_by_transcript_name_3(self):
        pos = transcript_to_exon("KRAS-201", 400, 420)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0].transcript_id, "ENST00000256078")
        self.assertEqual(pos[0].exon_id, "ENSE00001719809")
        self.assertEqual(pos[0].start, 25227234)
        self.assertEqual(pos[0].end, 25227412)
        self.assertEqual(pos[0].strand, "-")


class TestTranscriptToGene(unittest.TestCase):
    def test_transcript_to_gene_pos_strand_by_transcript_id(self):
        pos = transcript_to_gene("ENST00000341165", 1)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0].gene_id, "ENSG00000188554")
        self.assertEqual(pos[0].start, 43170481)
        self.assertEqual(pos[0].end, 43170481)
        self.assertEqual(pos[0].strand, "+")

    def test_transcript_to_gene_pos_strand_by_transcript_id_2(self):
        pos = transcript_to_gene("ENST00000341165", 731)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0].gene_id, "ENSG00000188554")
        self.assertEqual(pos[0].start, 43189698)
        self.assertEqual(pos[0].end, 43189698)
        self.assertEqual(pos[0].strand, "+")

    def test_transcript_to_gene_pos_strand_by_transcript_id_3(self):
        pos = transcript_to_gene("ENST00000341165", 731, 817)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0].gene_id, "ENSG00000188554")
        self.assertEqual(pos[0].start, 43189698)
        self.assertEqual(pos[0].end, 43189784)
        self.assertEqual(pos[0].strand, "+")

    def test_transcript_to_gene_pos_strand_by_transcript_name(self):
        pos = transcript_to_gene("NBR1-201", 1)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0].gene_id, "ENSG00000188554")
        self.assertEqual(pos[0].start, 43170481)
        self.assertEqual(pos[0].end, 43170481)
        self.assertEqual(pos[0].strand, "+")

    def test_transcript_to_gene_pos_strand_by_transcript_name_2(self):
        pos = transcript_to_gene("NBR1-201", 731)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0].gene_id, "ENSG00000188554")
        self.assertEqual(pos[0].start, 43189698)
        self.assertEqual(pos[0].end, 43189698)
        self.assertEqual(pos[0].strand, "+")

    def test_transcript_to_gene_pos_strand_by_transcript_name_3(self):
        pos = transcript_to_gene("NBR1-201", 731, 817)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0].gene_id, "ENSG00000188554")
        self.assertEqual(pos[0].start, 43189698)
        self.assertEqual(pos[0].end, 43189784)
        self.assertEqual(pos[0].strand, "+")

    def test_transcript_to_gene_neg_strand_by_transcript_id(self):
        pos = transcript_to_gene("ENST00000256078", 1)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0].gene_id, "ENSG00000133703")
        self.assertEqual(pos[0].start, 25250929)
        self.assertEqual(pos[0].end, 25250929)
        self.assertEqual(pos[0].strand, "-")

    def test_transcript_to_gene_neg_strand_by_transcript_id_2(self):
        pos = transcript_to_gene("ENST00000256078", 466)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0].gene_id, "ENSG00000133703")
        self.assertEqual(pos[0].start, 25227248)
        self.assertEqual(pos[0].end, 25227248)
        self.assertEqual(pos[0].strand, "-")

    def test_transcript_to_gene_neg_strand_by_transcript_id_3(self):
        pos = transcript_to_gene("ENST00000256078", 466, 501)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0].gene_id, "ENSG00000133703")
        self.assertEqual(pos[0].start, 25225753)
        self.assertEqual(pos[0].end, 25227248)
        self.assertEqual(pos[0].strand, "-")

    def test_transcript_to_gene_neg_strand_by_transcript_name(self):
        pos = transcript_to_gene("KRAS-201", 1)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0].gene_id, "ENSG00000133703")
        self.assertEqual(pos[0].start, 25250929)
        self.assertEqual(pos[0].end, 25250929)
        self.assertEqual(pos[0].strand, "-")

    def test_transcript_to_gene_neg_strand_by_transcript_name_2(self):
        pos = transcript_to_gene("KRAS-201", 466)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0].gene_id, "ENSG00000133703")
        self.assertEqual(pos[0].start, 25227248)
        self.assertEqual(pos[0].end, 25227248)
        self.assertEqual(pos[0].strand, "-")

    def test_transcript_to_gene_neg_strand_by_transcript_name_3(self):
        pos = transcript_to_gene("KRAS-201", 466, 501)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0].gene_id, "ENSG00000133703")
        self.assertEqual(pos[0].start, 25225753)
        self.assertEqual(pos[0].end, 25227248)
        self.assertEqual(pos[0].strand, "-")


class TestTranscriptToProtein(unittest.TestCase):
    def test_transcript_to_protein_pos_strand_by_transcript_id(self):
        pos = transcript_to_protein("ENST00000288135", 98)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0].protein_id, "ENSP00000288135")
        self.assertEqual(pos[0].start, 1)
        self.assertEqual(pos[0].end, 1)
        self.assertEqual(pos[0].strand, "+")

    def test_transcript_to_protein_pos_strand_by_transcript_id_2(self):
        pos = transcript_to_protein("ENST00000288135", 128)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0].protein_id, "ENSP00000288135")
        self.assertEqual(pos[0].start, 11)
        self.assertEqual(pos[0].end, 11)
        self.assertEqual(pos[0].strand, "+")

    def test_transcript_to_protein_pos_strand_by_transcript_id_3(self):
        pos = transcript_to_protein("ENST00000288135", 128, 145)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0].protein_id, "ENSP00000288135")
        self.assertEqual(pos[0].start, 11)
        self.assertEqual(pos[0].end, 16)
        self.assertEqual(pos[0].strand, "+")

    def test_transcript_to_protein_pos_strand_by_transcript_name(self):
        pos = transcript_to_protein("KIT-201", 98)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0].protein_id, "ENSP00000288135")
        self.assertEqual(pos[0].start, 1)
        self.assertEqual(pos[0].end, 1)
        self.assertEqual(pos[0].strand, "+")

    def test_transcript_to_protein_pos_strand_by_transcript_name_2(self):
        pos = transcript_to_protein("KIT-201", 128)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0].protein_id, "ENSP00000288135")
        self.assertEqual(pos[0].start, 11)
        self.assertEqual(pos[0].end, 11)
        self.assertEqual(pos[0].strand, "+")

    def test_transcript_to_protein_pos_strand_by_transcript_name_3(self):
        pos = transcript_to_protein("KIT-201", 128, 145)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0].protein_id, "ENSP00000288135")
        self.assertEqual(pos[0].start, 11)
        self.assertEqual(pos[0].end, 16)
        self.assertEqual(pos[0].strand, "+")

    def test_transcript_to_protein_neg_strand_by_transcript_id(self):
        pos = transcript_to_protein("ENST00000260947", 115)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0].protein_id, "ENSP00000260947")
        self.assertEqual(pos[0].start, 1)
        self.assertEqual(pos[0].end, 1)
        self.assertEqual(pos[0].strand, "-")

    def test_transcript_to_protein_neg_strand_by_transcript_id_2(self):
        pos = transcript_to_protein("ENST00000260947", 346)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0].protein_id, "ENSP00000260947")
        self.assertEqual(pos[0].start, 78)
        self.assertEqual(pos[0].end, 78)
        self.assertEqual(pos[0].strand, "-")

    def test_transcript_to_protein_neg_strand_by_transcript_id_3(self):
        pos = transcript_to_protein("ENST00000260947", 345, 356)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0].protein_id, "ENSP00000260947")
        self.assertEqual(pos[0].start, 77)
        self.assertEqual(pos[0].end, 81)
        self.assertEqual(pos[0].strand, "-")

    def test_transcript_to_protein_neg_strand_by_transcript_name(self):
        pos = transcript_to_protein("BARD1-201", 115)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0].protein_id, "ENSP00000260947")
        self.assertEqual(pos[0].start, 1)
        self.assertEqual(pos[0].end, 1)
        self.assertEqual(pos[0].strand, "-")

    def test_transcript_to_protein_neg_strand_by_transcript_name_2(self):
        pos = transcript_to_protein("BARD1-201", 346)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0].protein_id, "ENSP00000260947")
        self.assertEqual(pos[0].start, 78)
        self.assertEqual(pos[0].end, 78)
        self.assertEqual(pos[0].strand, "-")

    def test_transcript_to_protein_neg_strand_by_transcript_name_3(self):
        pos = transcript_to_protein("BARD1-201", 345, 356)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0].protein_id, "ENSP00000260947")
        self.assertEqual(pos[0].start, 77)
        self.assertEqual(pos[0].end, 81)
        self.assertEqual(pos[0].strand, "-")
