import unittest

from ensembl_map.mapper import (
    cds_to_exon,
    cds_to_gene,
    cds_to_protein,
    cds_to_transcript,
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


class TestCdsToExon(unittest.TestCase):
    pass


class TestCdsToGene(unittest.TestCase):
    pass


class TestCdsToProtein(unittest.TestCase):
    def test_cds_to_protein_by_transcript_id(self):
        pos = cds_to_protein("ENST00000288135", 1)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0][0], "ENSP00000288135")
        self.assertEqual(pos[0][1], 1)
        self.assertEqual(pos[0][2], 1)
        self.assertEqual(pos[0][3], "+")

    def test_cds_to_protein_by_transcript_id_2(self):
        pos = cds_to_protein("ENST00000288135", 213)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0][0], "ENSP00000288135")
        self.assertEqual(pos[0][1], 71)
        self.assertEqual(pos[0][2], 71)
        self.assertEqual(pos[0][3], "+")

    def test_cds_to_protein_by_transcript_id_3(self):
        pos = cds_to_protein("ENST00000288135", 213, 221)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0][0], "ENSP00000288135")
        self.assertEqual(pos[0][1], 71)
        self.assertEqual(pos[0][2], 74)
        self.assertEqual(pos[0][3], "+")

    def test_cds_to_protein_by_transcript_name(self):
        pos = cds_to_protein("KIT-001", 1)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0][0], "ENSP00000288135")
        self.assertEqual(pos[0][1], 1)
        self.assertEqual(pos[0][2], 1)
        self.assertEqual(pos[0][3], "+")

    def test_cds_to_protein_by_transcript_name_2(self):
        pos = cds_to_protein("KIT-001", 213)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0][0], "ENSP00000288135")
        self.assertEqual(pos[0][1], 71)
        self.assertEqual(pos[0][2], 71)
        self.assertEqual(pos[0][3], "+")

    def test_cds_to_protein_by_transcript_name_3(self):
        pos = cds_to_protein("KIT-001", 213, 221)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0][0], "ENSP00000288135")
        self.assertEqual(pos[0][1], 71)
        self.assertEqual(pos[0][2], 74)
        self.assertEqual(pos[0][3], "+")


class TestCdsToTranscript(unittest.TestCase):
    def test_cds_to_transcript_by_transcript_id(self):
        pos = cds_to_transcript("ENST00000288135", 1)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0][0], "ENST00000288135")
        self.assertEqual(pos[0][1], 97)
        self.assertEqual(pos[0][2], 97)
        self.assertEqual(pos[0][3], "+")

    def test_cds_to_transcript_by_transcript_id_2(self):
        pos = cds_to_transcript("ENST00000288135", 213)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0][0], "ENST00000288135")
        self.assertEqual(pos[0][1], 309)
        self.assertEqual(pos[0][2], 309)
        self.assertEqual(pos[0][3], "+")

    def test_cds_to_transcript_by_transcript_id_3(self):
        pos = cds_to_transcript("ENST00000288135", 213, 221)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0][0], "ENST00000288135")
        self.assertEqual(pos[0][1], 309)
        self.assertEqual(pos[0][2], 317)
        self.assertEqual(pos[0][3], "+")

    def test_cds_to_transcript_by_transcript_name(self):
        pos = cds_to_transcript("KIT-001", 1)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0][0], "ENST00000288135")
        self.assertEqual(pos[0][1], 97)
        self.assertEqual(pos[0][2], 97)
        self.assertEqual(pos[0][3], "+")

    def test_cds_to_transcript_by_transcript_name_2(self):
        pos = cds_to_transcript("KIT-001", 213)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0][0], "ENST00000288135")
        self.assertEqual(pos[0][1], 309)
        self.assertEqual(pos[0][2], 309)
        self.assertEqual(pos[0][3], "+")

    def test_cds_to_transcript_by_transcript_name_3(self):
        pos = cds_to_transcript("KIT-001", 213, 221)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0][0], "ENST00000288135")
        self.assertEqual(pos[0][1], 309)
        self.assertEqual(pos[0][2], 317)
        self.assertEqual(pos[0][3], "+")


class TestExonToCds(unittest.TestCase):
    pass


class TestExonToGene(unittest.TestCase):
    pass


class TestExonToProtein(unittest.TestCase):
    pass


# TODO: For some reason `transcript_ids_of_exon_ids` does not match anything.
# class TestExonToTranscript(unittest.TestCase):
#     def test_exon_to_transcript_by_transcript_id(self):
#         pos = exon_to_transcript("ENSE00001798125")
#         self.assertIsInstance(pos, list)
#         self.assertEqual(pos[0][0], "ENST00000275493")
#         self.assertEqual(pos[0][1], 55214299)
#         self.assertEqual(pos[0][2], 55214433)
#         self.assertEqual(pos[0][3], "+")

#     def test_exon_to_transcript_by_transcript_id_2(self):
#         pos = exon_to_transcript("ENSE00000936617")
#         self.assertIsInstance(pos, list)
#         self.assertEqual(pos[0][0], "ENST00000275493")
#         self.assertEqual(pos[0][1], 25398329)
#         self.assertEqual(pos[0][2], 25398208)
#         self.assertEqual(pos[0][3], "+")


class TestGeneToCds(unittest.TestCase):
    pass


class TestGeneToExon(unittest.TestCase):
    pass


class TestGeneToProtein(unittest.TestCase):
    pass


class TestGeneToTranscript(unittest.TestCase):
    def test_gene_to_transcript_by_transcript_id(self):
        pos = gene_to_transcript("ENSG00000188554", 41322498)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0][0], "ENST00000341165")
        self.assertEqual(pos[0][1], 1)
        self.assertEqual(pos[0][2], 1)
        self.assertEqual(pos[0][3], "+")

    def test_gene_to_transcript_by_transcript_id_2(self):
        pos = gene_to_transcript("ENSG00000188554", 41341715)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0][0], "ENST00000341165")
        self.assertEqual(pos[0][1], 731)
        self.assertEqual(pos[0][2], 731)
        self.assertEqual(pos[0][3], "+")

    def test_gene_to_transcript_by_transcript_id_3(self):
        pos = gene_to_transcript("ENSG00000188554", 41341715, 41341801)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0][0], "ENST00000341165")
        self.assertEqual(pos[0][1], 731)
        self.assertEqual(pos[0][2], 817)
        self.assertEqual(pos[0][3], "+")

    def test_gene_to_transcript_by_transcript_name(self):
        pos = gene_to_transcript("NBR1", 41322498)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0][0], "ENST00000341165")
        self.assertEqual(pos[0][1], 1)
        self.assertEqual(pos[0][2], 1)
        self.assertEqual(pos[0][3], "+")

    def test_gene_to_transcript_by_transcript_name_2(self):
        pos = gene_to_transcript("NBR1", 41341715)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0][0], "ENST00000341165")
        self.assertEqual(pos[0][1], 731)
        self.assertEqual(pos[0][2], 731)
        self.assertEqual(pos[0][3], "+")

    def test_gene_to_transcript_by_transcript_name_3(self):
        pos = gene_to_transcript("NBR1", 41341715, 41341801)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0][0], "ENST00000341165")
        self.assertEqual(pos[0][1], 731)
        self.assertEqual(pos[0][2], 817)
        self.assertEqual(pos[0][3], "+")


class TestProteinToCds(unittest.TestCase):
    def test_protein_to_cds_by_transcript_id(self):
        pos = protein_to_cds("ENSP00000308495", 1)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0][0], "ENST00000311936")
        self.assertEqual(pos[0][1], 1)
        self.assertEqual(pos[0][2], 3)
        self.assertEqual(pos[0][3], "-")

    def test_protein_to_cds_by_transcript_id_2(self):
        pos = protein_to_cds("ENSP00000308495", 123)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0][0], "ENST00000311936")
        self.assertEqual(pos[0][1], 367)
        self.assertEqual(pos[0][2], 369)
        self.assertEqual(pos[0][3], "-")

    def test_protein_to_cds_by_transcript_id_3(self):
        pos = protein_to_cds("ENSP00000308495", 123, 124)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0][0], "ENST00000311936")
        self.assertEqual(pos[0][1], 367)
        self.assertEqual(pos[0][2], 372)
        self.assertEqual(pos[0][3], "-")


class TestProteinToExon(unittest.TestCase):
    pass


class TestProteinToGene(unittest.TestCase):
    pass


class TestProteinToTranscript(unittest.TestCase):
    pass


class TestTranscriptToCds(unittest.TestCase):
    def test_transcript_to_cds_by_transcript_id(self):
        pos = transcript_to_cds("ENST00000288135", 97)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0][0], "ENST00000288135")
        self.assertEqual(pos[0][1], 1)
        self.assertEqual(pos[0][2], 1)
        self.assertEqual(pos[0][3], "+")

    def test_transcript_to_cds_by_transcript_id_2(self):
        pos = transcript_to_cds("ENST00000288135", 100)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0][0], "ENST00000288135")
        self.assertEqual(pos[0][1], 4)
        self.assertEqual(pos[0][2], 4)
        self.assertEqual(pos[0][3], "+")

    def test_transcript_to_cds_by_transcript_id_3(self):
        pos = transcript_to_cds("ENST00000288135", 97, 100)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0][0], "ENST00000288135")
        self.assertEqual(pos[0][1], 1)
        self.assertEqual(pos[0][2], 4)
        self.assertEqual(pos[0][3], "+")

    def test_transcript_to_cds_by_transcript_name(self):
        pos = transcript_to_cds("KIT-001", 97)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0][0], "ENST00000288135")
        self.assertEqual(pos[0][1], 1)
        self.assertEqual(pos[0][2], 1)
        self.assertEqual(pos[0][3], "+")

    def test_transcript_to_cds_by_transcript_name_2(self):
        pos = transcript_to_cds("KIT-001", 100)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0][0], "ENST00000288135")
        self.assertEqual(pos[0][1], 4)
        self.assertEqual(pos[0][2], 4)
        self.assertEqual(pos[0][3], "+")

    def test_transcript_to_cds_by_transcript_name_3(self):
        pos = transcript_to_cds("KIT-001", 97, 100)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0][0], "ENST00000288135")
        self.assertEqual(pos[0][1], 1)
        self.assertEqual(pos[0][2], 4)
        self.assertEqual(pos[0][3], "+")


class TestTranscriptToExon(unittest.TestCase):
    pass


class TestTranscriptToGene(unittest.TestCase):
    pass


class TestTranscriptToProtein(unittest.TestCase):
    pass

