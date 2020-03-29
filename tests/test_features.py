import unittest
from unittest.mock import MagicMock

from ensembl_map.features import CDS, Exon, Gene, Protein, Transcript


class TestCds(unittest.TestCase):
    def setUp(self):
        self.obj = MagicMock()
        self.cds = CDS(self.obj, 1, 2)

    def test_cds_transcript(self):
        self.assertEqual(self.cds.transcript, self.obj)

    def test_cds_transcript_id(self):
        self.assertEqual(self.cds.transcript_id, self.obj.transcript_id)

    def test_cds_transcript_id_none(self):
        self.cds._pyensembl_obj = None  # causes AttributeError
        self.assertEqual(self.cds.transcript_id, None)

    def test_cds_transcript_name(self):
        self.assertEqual(self.cds.transcript_name, self.obj.transcript_name)

    def test_cds_transcript_name_none(self):
        self.cds._pyensembl_obj = None  # causes AttributeError
        self.assertEqual(self.cds.transcript_name, None)

    def test_cds_transcript_sequence(self):
        self.assertEqual(self.cds.sequence, self.obj.coding_sequence.__getitem__())

    def test_cds_transcript_to_tuple(self):
        self.assertEqual(self.cds.to_tuple(), (self.obj.transcript_id, 1, 2))


class TestExon(unittest.TestCase):
    def setUp(self):
        self.eobj = MagicMock()
        self.tobj = MagicMock()
        self.exon = Exon(self.tobj, self.eobj, 1, 2)

    def test_exon_exon(self):
        self.assertEqual(self.exon.exon, self.eobj)

    def test_exon_exon_id(self):
        self.assertEqual(self.exon.exon_id, self.eobj.exon_id)

    def test_exon_exon_id_none(self):
        self.exon._exon = None  # causes AttributeError
        self.assertEqual(self.exon.exon_id, None)

