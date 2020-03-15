import unittest

from ensembl_map import set_ensembl_release
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


# TODO: Mock Cache object instead of needing to download annotations to run tests.
set_ensembl_release(69, "human")


class TestCdsToExon(unittest.TestCase):
    # TODO: def test_cds_to_exon_pos_strand_by_transcript_id(self):

    # TODO: def test_cds_to_exon_pos_strand_by_transcript_id_2(self):

    # TODO: def test_cds_to_exon_pos_strand_by_transcript_id_3(self):

    # TODO: def test_cds_to_exon_pos_strand_by_transcript_name(self):

    # TODO: def test_cds_to_exon_pos_strand_by_transcript_name_2(self):

    # TODO: def test_cds_to_exon_pos_strand_by_transcript_name_3(self):

    # TODO: def test_cds_to_exon_neg_strand_by_transcript_id(self):

    # TODO: def test_cds_to_exon_neg_strand_by_transcript_id_2(self):

    # TODO: def test_cds_to_exon_neg_strand_by_transcript_id_3(self):

    # TODO: def test_cds_to_exon_neg_strand_by_transcript_name(self):

    # TODO: def test_cds_to_exon_neg_strand_by_transcript_name_2(self):

    # TODO: def test_cds_to_exon_neg_strand_by_transcript_name_3(self):
    pass


class TestCdsToGene(unittest.TestCase):
    def test_cds_to_gene_pos_strand_by_transcript_id(self):
        pos = cds_to_gene("ENST00000288135", 1)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0][0], "ENSG00000157404")
        self.assertEqual(pos[0][1], 55524181)
        self.assertEqual(pos[0][2], 55524181)
        self.assertEqual(pos[0][3], "+")

    def test_cds_to_gene_pos_strand_by_transcript_id_2(self):
        pos = cds_to_gene("ENST00000288135", 213)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0][0], "ENSG00000157404")
        self.assertEqual(pos[0][1], 55561822)
        self.assertEqual(pos[0][2], 55561822)
        self.assertEqual(pos[0][3], "+")

    def test_cds_to_gene_pos_strand_by_transcript_id_3(self):
        pos = cds_to_gene("ENST00000288135", 213, 221)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0][0], "ENSG00000157404")
        self.assertEqual(pos[0][1], 55561822)
        self.assertEqual(pos[0][2], 55561830)
        self.assertEqual(pos[0][3], "+")

    def test_cds_to_gene_pos_strand_by_transcript_name(self):
        pos = cds_to_gene("KIT-001", 1)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0][0], "ENSG00000157404")
        self.assertEqual(pos[0][1], 55524181)
        self.assertEqual(pos[0][2], 55524181)
        self.assertEqual(pos[0][3], "+")

    def test_cds_to_gene_pos_strand_by_transcript_name_2(self):
        pos = cds_to_gene("KIT-001", 213)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0][0], "ENSG00000157404")
        self.assertEqual(pos[0][1], 55561822)
        self.assertEqual(pos[0][2], 55561822)
        self.assertEqual(pos[0][3], "+")

    def test_cds_to_gene_pos_strand_by_transcript_name_3(self):
        pos = cds_to_gene("KIT-001", 213, 221)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0][0], "ENSG00000157404")
        self.assertEqual(pos[0][1], 55561822)
        self.assertEqual(pos[0][2], 55561830)
        self.assertEqual(pos[0][3], "+")

    # TODO: def test_cds_to_gene_neg_strand_by_transcript_id(self):

    # TODO: def test_cds_to_gene_neg_strand_by_transcript_id_2(self):

    # TODO: def test_cds_to_gene_neg_strand_by_transcript_id_3(self):

    # TODO: def test_cds_to_gene_neg_strand_by_transcript_name(self):

    # TODO: def test_cds_to_gene_neg_strand_by_transcript_name_2(self):

    # TODO: def test_cds_to_gene_neg_strand_by_transcript_name_3(self):


class TestCdsToProtein(unittest.TestCase):
    def test_cds_to_protein_pos_strand_by_transcript_id(self):
        pos = cds_to_protein("ENST00000288135", 1)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0][0], "ENSP00000288135")
        self.assertEqual(pos[0][1], 1)
        self.assertEqual(pos[0][2], 1)
        self.assertEqual(pos[0][3], "+")

    def test_cds_to_protein_pos_strand_by_transcript_id_2(self):
        pos = cds_to_protein("ENST00000288135", 213)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0][0], "ENSP00000288135")
        self.assertEqual(pos[0][1], 71)
        self.assertEqual(pos[0][2], 71)
        self.assertEqual(pos[0][3], "+")

    def test_cds_to_protein_pos_strand_by_transcript_id_3(self):
        pos = cds_to_protein("ENST00000288135", 213, 221)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0][0], "ENSP00000288135")
        self.assertEqual(pos[0][1], 71)
        self.assertEqual(pos[0][2], 74)
        self.assertEqual(pos[0][3], "+")

    def test_cds_to_protein_pos_strand_by_transcript_name(self):
        pos = cds_to_protein("KIT-001", 1)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0][0], "ENSP00000288135")
        self.assertEqual(pos[0][1], 1)
        self.assertEqual(pos[0][2], 1)
        self.assertEqual(pos[0][3], "+")

    def test_cds_to_protein_pos_strand_by_transcript_name_2(self):
        pos = cds_to_protein("KIT-001", 213)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0][0], "ENSP00000288135")
        self.assertEqual(pos[0][1], 71)
        self.assertEqual(pos[0][2], 71)
        self.assertEqual(pos[0][3], "+")

    def test_cds_to_protein_pos_strand_by_transcript_name_3(self):
        pos = cds_to_protein("KIT-001", 213, 221)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0][0], "ENSP00000288135")
        self.assertEqual(pos[0][1], 71)
        self.assertEqual(pos[0][2], 74)
        self.assertEqual(pos[0][3], "+")

    # TODO: def test_cds_to_protein_neg_strand_by_transcript_id(self):

    # TODO: def test_cds_to_protein_neg_strand_by_transcript_id_2(self):

    # TODO: def test_cds_to_protein_neg_strand_by_transcript_id_3(self):

    # TODO: def test_cds_to_protein_neg_strand_by_transcript_name(self):

    # TODO: def test_cds_to_protein_neg_strand_by_transcript_name_2(self):

    # TODO: def test_cds_to_protein_neg_strand_by_transcript_name_3(self):


class TestCdsToTranscript(unittest.TestCase):
    def test_cds_to_transcript_pos_strand_by_transcript_id(self):
        pos = cds_to_transcript("ENST00000288135", 1)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0][0], "ENST00000288135")
        self.assertEqual(pos[0][1], 97)
        self.assertEqual(pos[0][2], 97)
        self.assertEqual(pos[0][3], "+")

    def test_cds_to_transcript_pos_strand_by_transcript_id_2(self):
        pos = cds_to_transcript("ENST00000288135", 213)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0][0], "ENST00000288135")
        self.assertEqual(pos[0][1], 309)
        self.assertEqual(pos[0][2], 309)
        self.assertEqual(pos[0][3], "+")

    def test_cds_to_transcript_pos_strand_by_transcript_id_3(self):
        pos = cds_to_transcript("ENST00000288135", 213, 221)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0][0], "ENST00000288135")
        self.assertEqual(pos[0][1], 309)
        self.assertEqual(pos[0][2], 317)
        self.assertEqual(pos[0][3], "+")

    def test_cds_to_transcript_pos_strand_by_transcript_name(self):
        pos = cds_to_transcript("KIT-001", 1)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0][0], "ENST00000288135")
        self.assertEqual(pos[0][1], 97)
        self.assertEqual(pos[0][2], 97)
        self.assertEqual(pos[0][3], "+")

    def test_cds_to_transcript_pos_strand_by_transcript_name_2(self):
        pos = cds_to_transcript("KIT-001", 213)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0][0], "ENST00000288135")
        self.assertEqual(pos[0][1], 309)
        self.assertEqual(pos[0][2], 309)
        self.assertEqual(pos[0][3], "+")

    def test_cds_to_transcript_pos_strand_by_transcript_name_3(self):
        pos = cds_to_transcript("KIT-001", 213, 221)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0][0], "ENST00000288135")
        self.assertEqual(pos[0][1], 309)
        self.assertEqual(pos[0][2], 317)
        self.assertEqual(pos[0][3], "+")

    # TODO: def test_cds_to_transcript_neg_strand_by_transcript_id(self):

    # TODO: def test_cds_to_transcript_neg_strand_by_transcript_id_2(self):

    # TODO: def test_cds_to_transcript_neg_strand_by_transcript_id_3(self):

    # TODO: def test_cds_to_transcript_neg_strand_by_transcript_name(self):

    # TODO: def test_cds_to_transcript_neg_strand_by_transcript_name_2(self):

    # TODO: def test_cds_to_transcript_neg_strand_by_transcript_name_3(self):


class TestExonToCds(unittest.TestCase):
    # TODO: def test_exon_to_cds_pos_strand_by_exon_id(self):

    # TODO: def test_exon_to_cds_pos_strand_by_exon_id_2(self):

    # TODO: def test_exon_to_cds_pos_strand_by_exon_id_3(self)::

    # TODO: def test_exon_to_cds_neg_strand_by_exon_id(self):

    # TODO: def test_exon_to_cds_neg_strand_by_exon_id_2(self):

    # TODO: def test_exon_to_cds_neg_strand_by_exon_id_3(self):
    pass


class TestExonToGene(unittest.TestCase):
    # TODO: def test_exon_to_gene_pos_strand_by_exon_id(self):

    # TODO: def test_exon_to_gene_pos_strand_by_exon_id_2(self):

    # TODO: def test_exon_to_gene_pos_strand_by_exon_id_3(self)::

    # TODO: def test_exon_to_gene_neg_strand_by_exon_id(self):

    # TODO: def test_exon_to_gene_neg_strand_by_exon_id_2(self):

    # TODO: def test_exon_to_gene_neg_strand_by_exon_id_3(self):
    pass


class TestExonToProtein(unittest.TestCase):
    # TODO: def test_exon_to_protein_pos_strand_by_exon_id(self):

    # TODO: def test_exon_to_protein_pos_strand_by_exon_id_2(self):

    # TODO: def test_exon_to_protein_pos_strand_by_exon_id_3(self)::

    # TODO: def test_exon_to_protein_neg_strand_by_exon_id(self):

    # TODO: def test_exon_to_protein_neg_strand_by_exon_id_2(self):

    # TODO: def test_exon_to_protein_neg_strand_by_exon_id_3(self):
    pass


class TestExonToTranscript(unittest.TestCase):
    # TODO: def test_exon_to_transcript_pos_strand_by_exon_id(self):

    # TODO: def test_exon_to_transcript_pos_strand_by_exon_id_2(self):

    # TODO: def test_exon_to_transcript_pos_strand_by_exon_id_3(self)::

    # TODO: def test_exon_to_transcript_neg_strand_by_exon_id(self):

    # TODO: def test_exon_to_transcript_neg_strand_by_exon_id_2(self):

    # TODO: def test_exon_to_transcript_neg_strand_by_exon_id_3(self):
    pass


class TestGeneToCds(unittest.TestCase):
    def test_gene_to_cds_pos_strand_by_gene_id(self):
        pos = gene_to_cds("ENSG00000157404", 55524181)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0][0], "ENST00000288135")
        self.assertEqual(pos[0][1], 1)
        self.assertEqual(pos[0][2], 1)
        self.assertEqual(pos[0][3], "+")

    def test_gene_to_cds_pos_strand_by_gene_id_2(self):
        pos = gene_to_cds("ENSG00000157404", 55561822)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0][0], "ENST00000288135")
        self.assertEqual(pos[0][1], 213)
        self.assertEqual(pos[0][2], 213)
        self.assertEqual(pos[0][3], "+")

    def test_gene_to_cds_pos_strand_by_gene_id_3(self):
        pos = gene_to_cds("ENSG00000157404", 55561822, 55561830)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0][0], "ENST00000288135")
        self.assertEqual(pos[0][1], 213)
        self.assertEqual(pos[0][2], 221)
        self.assertEqual(pos[0][3], "+")

    def test_gene_to_cds_pos_strand_by_gene_name(self):
        pos = gene_to_cds("KIT", 55524181)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0][0], "ENST00000288135")
        self.assertEqual(pos[0][1], 1)
        self.assertEqual(pos[0][2], 1)
        self.assertEqual(pos[0][3], "+")

    def test_gene_to_cds_pos_strand_by_gene_name_2(self):
        pos = gene_to_cds("KIT", 55561822)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0][0], "ENST00000288135")
        self.assertEqual(pos[0][1], 213)
        self.assertEqual(pos[0][2], 213)
        self.assertEqual(pos[0][3], "+")

    def test_gene_to_cds_pos_strand_by_gene_name_3(self):
        pos = gene_to_cds("KIT", 55561822, 55561830)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0][0], "ENST00000288135")
        self.assertEqual(pos[0][1], 213)
        self.assertEqual(pos[0][2], 221)
        self.assertEqual(pos[0][3], "+")

    # TODO: def test_gene_to_cds_neg_strand_by_gene_id(self):

    # TODO: def test_gene_to_cds_neg_strand_by_gene_id_2(self):

    # TODO: def test_gene_to_cds_neg_strand_by_gene_id_3(self):

    # TODO: def test_gene_to_cds_neg_strand_by_gene_name(self):

    # TODO: def test_gene_to_cds_neg_strand_by_gene_name_2(self):

    # TODO: def test_gene_to_cds_neg_strand_by_gene_name_3(self):


class TestGeneToExon(unittest.TestCase):
    # TODO: def test_gene_to_exon_pos_strand_by_gene_id(self):

    # TODO: def test_gene_to_exon_pos_strand_by_gene_id_2(self):

    # TODO: def test_gene_to_exon_pos_strand_by_gene_id_3(self):

    # TODO: def test_gene_to_exon_pos_strand_by_gene_name(self):

    # TODO: def test_gene_to_exon_pos_strand_by_gene_name_2(self):

    # TODO: def test_gene_to_exon_pos_strand_by_gene_name_3(self):

    # TODO: def test_gene_to_exon_neg_strand_by_gene_id(self):

    # TODO: def test_gene_to_exon_neg_strand_by_gene_id_2(self):

    # TODO: def test_gene_to_exon_neg_strand_by_gene_id_3(self):

    # TODO: def test_gene_to_exon_neg_strand_by_gene_name(self):

    # TODO: def test_gene_to_exon_neg_strand_by_gene_name_2(self):

    # TODO: def test_gene_to_exon_neg_strand_by_gene_name_3(self):
    pass


class TestGeneToProtein(unittest.TestCase):
    def test_gene_to_protein_pos_strand_by_gene_id(self):
        pos = gene_to_protein("ENSG00000157404", 55524181)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0][0], "ENSP00000288135")
        self.assertEqual(pos[0][1], 1)
        self.assertEqual(pos[0][2], 1)
        self.assertEqual(pos[0][3], "+")

    def test_gene_to_protein_pos_strand_by_gene_id_2(self):
        pos = gene_to_protein("ENSG00000157404", 55561822)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0][0], "ENSP00000288135")
        self.assertEqual(pos[0][1], 71)
        self.assertEqual(pos[0][2], 71)
        self.assertEqual(pos[0][3], "+")

    def test_gene_to_protein_pos_strand_by_gene_id_3(self):
        pos = gene_to_protein("ENSG00000157404", 55561822, 55561830)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0][0], "ENSP00000288135")
        self.assertEqual(pos[0][1], 71)
        self.assertEqual(pos[0][2], 74)
        self.assertEqual(pos[0][3], "+")

    def test_gene_to_protein_pos_strand_by_gene_name(self):
        pos = gene_to_protein("KIT", 55524181)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0][0], "ENSP00000288135")
        self.assertEqual(pos[0][1], 1)
        self.assertEqual(pos[0][2], 1)
        self.assertEqual(pos[0][3], "+")

    def test_gene_to_protein_pos_strand_by_gene_name_2(self):
        pos = gene_to_protein("KIT", 55561822)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0][0], "ENSP00000288135")
        self.assertEqual(pos[0][1], 71)
        self.assertEqual(pos[0][2], 71)
        self.assertEqual(pos[0][3], "+")

    def test_gene_to_protein_pos_strand_by_gene_name_3(self):
        pos = gene_to_protein("KIT", 55561822, 55561830)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0][0], "ENSP00000288135")
        self.assertEqual(pos[0][1], 71)
        self.assertEqual(pos[0][2], 74)
        self.assertEqual(pos[0][3], "+")

    # TODO: def test_gene_to_cds_neg_strand_by_gene_id(self):

    # TODO: def test_gene_to_cds_neg_strand_by_gene_id_2(self):

    # TODO: def test_gene_to_cds_neg_strand_by_gene_id_3(self):

    # TODO: def test_gene_to_cds_neg_strand_by_gene_name(self):

    # TODO: def test_gene_to_cds_neg_strand_by_gene_name_2(self):

    # TODO: def test_gene_to_cds_neg_strand_by_gene_name_3(self):


class TestGeneToTranscript(unittest.TestCase):
    def test_gene_to_transcript_pos_strand_by_gene_id(self):
        pos = gene_to_transcript("ENSG00000188554", 41322498)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0][0], "ENST00000341165")
        self.assertEqual(pos[0][1], 1)
        self.assertEqual(pos[0][2], 1)
        self.assertEqual(pos[0][3], "+")

    def test_gene_to_transcript_pos_strand_by_gene_id_2(self):
        pos = gene_to_transcript("ENSG00000188554", 41341715)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0][0], "ENST00000341165")
        self.assertEqual(pos[0][1], 731)
        self.assertEqual(pos[0][2], 731)
        self.assertEqual(pos[0][3], "+")

    def test_gene_to_transcript_pos_strand_by_gene_id_3(self):
        pos = gene_to_transcript("ENSG00000188554", 41341715, 41341801)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0][0], "ENST00000341165")
        self.assertEqual(pos[0][1], 731)
        self.assertEqual(pos[0][2], 817)
        self.assertEqual(pos[0][3], "+")

    def test_gene_to_transcript_pos_strand_by_gene_name(self):
        pos = gene_to_transcript("NBR1", 41322498)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0][0], "ENST00000341165")
        self.assertEqual(pos[0][1], 1)
        self.assertEqual(pos[0][2], 1)
        self.assertEqual(pos[0][3], "+")

    def test_gene_to_transcript_pos_strand_by_gene_name_2(self):
        pos = gene_to_transcript("NBR1", 41341715)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0][0], "ENST00000341165")
        self.assertEqual(pos[0][1], 731)
        self.assertEqual(pos[0][2], 731)
        self.assertEqual(pos[0][3], "+")

    def test_gene_to_transcript_pos_strand_by_gene_name_3(self):
        pos = gene_to_transcript("NBR1", 41341715, 41341801)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0][0], "ENST00000341165")
        self.assertEqual(pos[0][1], 731)
        self.assertEqual(pos[0][2], 817)
        self.assertEqual(pos[0][3], "+")

    # TODO: def test_gene_to_transcript_neg_strand_by_gene_id(self):

    # TODO: def test_gene_to_transcript_neg_strand_by_gene_id_2(self):

    # TODO: def test_gene_to_transcript_neg_strand_by_gene_id_3(self):

    # TODO: def test_gene_to_transcript_neg_strand_by_gene_name(self):

    # TODO: def test_gene_to_transcript_neg_strand_by_gene_name_2(self):

    # TODO: def test_gene_to_transcript_neg_strand_by_gene_name_3(self):


class TestProteinToCds(unittest.TestCase):
    # TODO: def test_protein_to_cds_pos_strand_by_protein_id(self):

    # TODO: def test_protein_to_cds_pos_strand_by_protein_id_2(self):

    # TODO: def test_protein_to_cds_pos_strand_by_protein_id_3(self):

    def test_protein_to_cds_neg_strand_by_protein_id(self):
        pos = protein_to_cds("ENSP00000308495", 1)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0][0], "ENST00000311936")
        self.assertEqual(pos[0][1], 1)
        self.assertEqual(pos[0][2], 3)
        self.assertEqual(pos[0][3], "-")

    def test_protein_to_cds_neg_strand_by_protein_id_2(self):
        pos = protein_to_cds("ENSP00000308495", 123)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0][0], "ENST00000311936")
        self.assertEqual(pos[0][1], 367)
        self.assertEqual(pos[0][2], 369)
        self.assertEqual(pos[0][3], "-")

    def test_protein_to_cds_neg_strand_by_protein_id_3(self):
        pos = protein_to_cds("ENSP00000308495", 123, 124)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0][0], "ENST00000311936")
        self.assertEqual(pos[0][1], 367)
        self.assertEqual(pos[0][2], 372)
        self.assertEqual(pos[0][3], "-")


class TestProteinToExon(unittest.TestCase):
    # TODO: def test_protein_to_exon_pos_strand_by_protein_id(self):

    # TODO: def test_protein_to_exon_pos_strand_by_protein_id_2(self):

    # TODO: def test_protein_to_exon_pos_strand_by_protein_id_3(self):

    # TODO: def test_protein_to_exon_neg_strand_by_protein_id(self):

    # TODO: def test_protein_to_exon_neg_strand_by_protein_id_2(self):

    # TODO: def test_protein_to_exon_neg_strand_by_protein_id_3(self):
    pass


class TestProteinToGene(unittest.TestCase):
    def test_protein_to_gene_pos_strand_by_protein_id(self):
        pos = protein_to_gene("ENSP00000288135", 1)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0][0], "ENSG00000157404")
        self.assertEqual(pos[0][1], 55524181)
        self.assertEqual(pos[0][2], 55524183)
        self.assertEqual(pos[0][3], "+")

    def test_protein_to_gene_pos_strand_by_protein_id_2(self):
        pos = protein_to_gene("ENSP00000288135", 71)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0][0], "ENSG00000157404")
        self.assertEqual(pos[0][1], 55561820)
        self.assertEqual(pos[0][2], 55561822)
        self.assertEqual(pos[0][3], "+")

    def test_protein_to_gene_pos_strand_by_protein_id_3(self):
        pos = protein_to_gene("ENSP00000288135", 71, 74)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0][0], "ENSG00000157404")
        self.assertEqual(pos[0][1], 55561820)
        self.assertEqual(pos[0][2], 55561831)
        self.assertEqual(pos[0][3], "+")

    # TODO: def test_protein_to_gene_neg_strand_by_protein_id(self):

    # TODO: def test_protein_to_gene_neg_strand_by_protein_id_2(self):

    # TODO: def test_protein_to_gene_neg_strand_by_protein_id_3(self):


class TestProteinToTranscript(unittest.TestCase):
    # TODO: def test_protein_to_transcript_pos_strand_by_protein_id(self):

    # TODO: def test_protein_to_transcript_pos_strand_by_protein_id_2(self):

    # TODO: def test_protein_to_transcript_pos_strand_by_protein_id_3(self):

    def test_protein_to_transcript_neg_strand_by_protein_id(self):
        pos = protein_to_transcript("ENSP00000260947", 1)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0][0], "ENST00000260947")
        self.assertEqual(pos[0][1], 135)
        self.assertEqual(pos[0][2], 137)
        self.assertEqual(pos[0][3], "-")

    def test_protein_to_transcript_neg_strand_by_protein_id_2(self):
        pos = protein_to_transcript("ENSP00000260947", 71)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0][0], "ENST00000260947")
        self.assertEqual(pos[0][1], 345)
        self.assertEqual(pos[0][2], 347)
        self.assertEqual(pos[0][3], "-")

    def test_protein_to_transcript_neg_strand_by_protein_id_3(self):
        pos = protein_to_transcript("ENSP00000260947", 71, 74)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0][0], "ENST00000260947")
        self.assertEqual(pos[0][1], 345)
        self.assertEqual(pos[0][2], 356)
        self.assertEqual(pos[0][3], "-")


class TestTranscriptToCds(unittest.TestCase):
    def test_transcript_to_cds_pos_strand_by_transcript_id(self):
        pos = transcript_to_cds("ENST00000288135", 97)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0][0], "ENST00000288135")
        self.assertEqual(pos[0][1], 1)
        self.assertEqual(pos[0][2], 1)
        self.assertEqual(pos[0][3], "+")

    def test_transcript_to_cds_pos_strand_by_transcript_id_2(self):
        pos = transcript_to_cds("ENST00000288135", 100)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0][0], "ENST00000288135")
        self.assertEqual(pos[0][1], 4)
        self.assertEqual(pos[0][2], 4)
        self.assertEqual(pos[0][3], "+")

    def test_transcript_to_cds_pos_strand_by_transcript_id_3(self):
        pos = transcript_to_cds("ENST00000288135", 97, 100)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0][0], "ENST00000288135")
        self.assertEqual(pos[0][1], 1)
        self.assertEqual(pos[0][2], 4)
        self.assertEqual(pos[0][3], "+")

    def test_transcript_to_cds_pos_strand_by_transcript_name(self):
        pos = transcript_to_cds("KIT-001", 97)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0][0], "ENST00000288135")
        self.assertEqual(pos[0][1], 1)
        self.assertEqual(pos[0][2], 1)
        self.assertEqual(pos[0][3], "+")

    def test_transcript_to_cds_pos_strand_by_transcript_name_2(self):
        pos = transcript_to_cds("KIT-001", 100)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0][0], "ENST00000288135")
        self.assertEqual(pos[0][1], 4)
        self.assertEqual(pos[0][2], 4)
        self.assertEqual(pos[0][3], "+")

    def test_transcript_to_cds_pos_strand_by_transcript_name_3(self):
        pos = transcript_to_cds("KIT-001", 97, 100)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0][0], "ENST00000288135")
        self.assertEqual(pos[0][1], 1)
        self.assertEqual(pos[0][2], 4)
        self.assertEqual(pos[0][3], "+")

    # TODO: def test_transcript_to_cds_neg_strand_by_transcript_id(self):

    # TODO: def test_transcript_to_cds_neg_strand_by_transcript_id_2(self):

    # TODO: def test_transcript_to_cds_neg_strand_by_transcript_id_3(self):

    # TODO: def test_transcript_to_cds_neg_strand_by_transcript_name(self):

    # TODO: def test_transcript_to_cds_neg_strand_by_transcript_name_2(self):

    # TODO: def test_transcript_to_cds_neg_strand_by_transcript_name_3(self):


class TestTranscriptToExon(unittest.TestCase):
    # TODO: def test_transcript_to_exon_pos_strand_by_transcript_id(self):

    # TODO: def test_transcript_to_exon_pos_strand_by_transcript_id_2(self):

    # TODO: def test_transcript_to_exon_pos_strand_by_transcript_id_3(self):

    # TODO: def test_transcript_to_exon_pos_strand_by_transcript_name(self):

    # TODO: def test_transcript_to_exon_pos_strand_by_transcript_name_2(self):

    # TODO: def test_transcript_to_exon_pos_strand_by_transcript_name_3(self):

    # TODO: def test_transcript_to_exon_neg_strand_by_transcript_id(self):

    # TODO: def test_transcript_to_exon_neg_strand_by_transcript_id_2(self):

    # TODO: def test_transcript_to_exon_neg_strand_by_transcript_id_3(self):

    # TODO: def test_transcript_to_exon_neg_strand_by_transcript_name(self):

    # TODO: def test_transcript_to_exon_neg_strand_by_transcript_name_2(self):

    # TODO: def test_transcript_to_exon_neg_strand_by_transcript_name_3(self):
    pass


class TestTranscriptToGene(unittest.TestCase):
    def test_transcript_to_gene_pos_strand_by_transcript_id(self):
        pos = transcript_to_gene("ENST00000341165", 1)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0][0], "ENSG00000188554")
        self.assertEqual(pos[0][1], 41322498)
        self.assertEqual(pos[0][2], 41322498)
        self.assertEqual(pos[0][3], "+")

    def test_transcript_to_gene_pos_strand_by_transcript_id_2(self):
        pos = transcript_to_gene("ENST00000341165", 731)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0][0], "ENSG00000188554")
        self.assertEqual(pos[0][1], 41341715)
        self.assertEqual(pos[0][2], 41341715)
        self.assertEqual(pos[0][3], "+")

    def test_transcript_to_gene_pos_strand_by_transcript_id_3(self):
        pos = transcript_to_gene("ENST00000341165", 731, 817)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0][0], "ENSG00000188554")
        self.assertEqual(pos[0][1], 41341715)
        self.assertEqual(pos[0][2], 41341801)
        self.assertEqual(pos[0][3], "+")

    def test_transcript_to_gene_pos_strand_by_transcript_name(self):
        pos = transcript_to_gene("NBR1-001", 1)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0][0], "ENSG00000188554")
        self.assertEqual(pos[0][1], 41322498)
        self.assertEqual(pos[0][2], 41322498)
        self.assertEqual(pos[0][3], "+")

    def test_transcript_to_gene_pos_strand_by_transcript_name_2(self):
        pos = transcript_to_gene("NBR1-001", 731)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0][0], "ENSG00000188554")
        self.assertEqual(pos[0][1], 41341715)
        self.assertEqual(pos[0][2], 41341715)
        self.assertEqual(pos[0][3], "+")

    def test_transcript_to_gene_pos_strand_by_transcript_name_3(self):
        pos = transcript_to_gene("NBR1-001", 731, 817)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0][0], "ENSG00000188554")
        self.assertEqual(pos[0][1], 41341715)
        self.assertEqual(pos[0][2], 41341801)
        self.assertEqual(pos[0][3], "+")

    # TODO: def test_transcript_to_gene_neg_strand_by_transcript_id(self):

    # TODO: def test_transcript_to_gene_neg_strand_by_transcript_id_2(self):

    # TODO: def test_transcript_to_gene_neg_strand_by_transcript_id_3(self):

    # TODO: def test_transcript_to_gene_neg_strand_by_transcript_name(self):

    # TODO: def test_transcript_to_gene_neg_strand_by_transcript_name_2(self):

    # TODO: def test_transcript_to_gene_neg_strand_by_transcript_name_3(self):


class TestTranscriptToProtein(unittest.TestCase):
    # TODO: def test_transcript_to_protein_pos_strand_by_transcript_id(self):

    # TODO: def test_transcript_to_protein_pos_strand_by_transcript_id_2(self):

    # TODO: def test_transcript_to_protein_pos_strand_by_transcript_id_3(self):

    # TODO: def test_transcript_to_protein_pos_strand_by_transcript_name(self):

    # TODO: def test_transcript_to_protein_pos_strand_by_transcript_name_2(self):

    # TODO: def test_transcript_to_protein_pos_strand_by_transcript_name_3(self):

    def test_transcript_to_protein_neg_strand_by_transcript_id(self):
        pos = transcript_to_protein("ENST00000260947", 136)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0][0], "ENSP00000260947")
        self.assertEqual(pos[0][1], 1)
        self.assertEqual(pos[0][2], 1)
        self.assertEqual(pos[0][3], "-")

    def test_transcript_to_protein_neg_strand_by_transcript_id_2(self):
        pos = transcript_to_protein("ENST00000260947", 346)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0][0], "ENSP00000260947")
        self.assertEqual(pos[0][1], 71)
        self.assertEqual(pos[0][2], 71)
        self.assertEqual(pos[0][3], "-")

    def test_transcript_to_protein_neg_strand_by_transcript_id_3(self):
        pos = transcript_to_protein("ENST00000260947", 345, 356)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0][0], "ENSP00000260947")
        self.assertEqual(pos[0][1], 71)
        self.assertEqual(pos[0][2], 74)
        self.assertEqual(pos[0][3], "-")

    def test_transcript_to_protein_neg_strand_by_transcript_name(self):
        pos = transcript_to_protein("BARD1-001", 136)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0][0], "ENSP00000260947")
        self.assertEqual(pos[0][1], 1)
        self.assertEqual(pos[0][2], 1)
        self.assertEqual(pos[0][3], "-")

    def test_transcript_to_protein_neg_strand_by_transcript_name_2(self):
        pos = transcript_to_protein("BARD1-001", 346)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0][0], "ENSP00000260947")
        self.assertEqual(pos[0][1], 71)
        self.assertEqual(pos[0][2], 71)
        self.assertEqual(pos[0][3], "-")

    def test_transcript_to_protein_neg_strand_by_transcript_name_3(self):
        pos = transcript_to_protein("BARD1-001", 345, 356)
        self.assertIsInstance(pos, list)
        self.assertEqual(pos[0][0], "ENSP00000260947")
        self.assertEqual(pos[0][1], 71)
        self.assertEqual(pos[0][2], 74)
        self.assertEqual(pos[0][3], "-")
