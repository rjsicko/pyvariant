import unittest

from ensembl_map.cache import set_ensembl_release
from ensembl_map.transcript import get_transcripts


# TODO: Mock Cache object instead of needing to download annotations to run tests.
set_ensembl_release(release=99, species="human", download_if_missing=True)


class TestTranscript(unittest.TestCase):
    def test_from_cds(self):
        ret_by_id = get_transcripts("ENST00000288135", "cds")
        ret_by_name = get_transcripts("KIT-201", "cds")
        self.assertTrue("ENST00000288135" in [i.transcript_id for i in ret_by_id])
        self.assertTrue("ENST00000288135" in [i.transcript_id for i in ret_by_name])
        self.assertEqual(ret_by_id, ret_by_name)

    def test_from_exon(self):
        ret_by_id = get_transcripts("ENSE00003659301", "exon")
        self.assertTrue("ENST00000380152" in [i.transcript_id for i in ret_by_id])

    def test_from_gene(self):
        ret_by_id = get_transcripts("ENSG00000157404", "gene")
        ret_by_name = get_transcripts("KIT", "gene")
        self.assertTrue("ENST00000288135" in [i.transcript_id for i in ret_by_id])
        self.assertTrue("ENST00000288135" in [i.transcript_id for i in ret_by_name])
        self.assertEqual(ret_by_id, ret_by_name)

    def test_from_protein(self):
        ret_by_id = get_transcripts("ENSP00000288135", "protein")
        self.assertTrue("ENST00000288135" in [i.transcript_id for i in ret_by_id])

    def test_from_transcript(self):
        ret_by_id = get_transcripts("ENST00000288135", "transcript")
        ret_by_name = get_transcripts("KIT-201", "transcript")
        self.assertTrue("ENST00000288135" in [i.transcript_id for i in ret_by_id])
        self.assertTrue("ENST00000288135" in [i.transcript_id for i in ret_by_name])
        self.assertEqual(ret_by_id, ret_by_name)
