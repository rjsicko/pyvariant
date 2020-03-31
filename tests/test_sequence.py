import unittest

from ensembl_map import set_ensembl_release
from ensembl_map.sequence import cds_sequence, protein_sequence, transcript_sequence


# TODO: Mock Cache object instead of needing to download annotations to run tests.
set_ensembl_release(release=99, species="human", download_if_missing=True)


class TestSequence(unittest.TestCase):
    def test_cds_sequence(self):
        self.assertEqual(cds_sequence("ENST00000288135", 1), "A")

    def test_cds_sequence_2(self):
        self.assertEqual(cds_sequence("ENST00000288135", 1, 3), "ATG")

    def test_protein_sequence(self):
        self.assertEqual(protein_sequence("ENSP00000288135", 1), "M")

    def test_protein_sequence_2(self):
        self.assertEqual(protein_sequence("ENSP00000288135", 3, 7), "GARGA")

    def test_transcript_sequence(self):
        self.assertEqual(transcript_sequence("ENST00000288135", 1), "C")

    def test_transcript_sequence_2(self):
        self.assertEqual(transcript_sequence("ENST00000288135", 166, 168), "CTC")
