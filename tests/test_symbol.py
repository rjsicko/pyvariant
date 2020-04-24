import unittest

from ensembl_map.ensembl import set_ensembl_release
from ensembl_map.symbol import (
    get_exons,
    get_exon_ids,
    get_genes,
    get_gene_ids,
    get_protein_ids,
    get_transcripts,
    get_transcript_ids,
)


# TODO: Mock Cache object instead of needing to download annotations to run tests.
set_ensembl_release(release=99, species="human", download_if_missing=True)


class TestExon(unittest.TestCase):
    def test_exon_ids_from_cds(self):
        ret_by_id = get_exon_ids("ENST00000288135", "cds")
        ret_by_name = get_exon_ids("KIT-201", "cds")
        self.assertTrue("ENSE00001074448" in ret_by_id)
        self.assertTrue("ENSE00001074448" in ret_by_name)
        self.assertEqual(ret_by_id, ret_by_name)

    def test_exon_ids_from_exon(self):
        ret_by_id = get_exon_ids("ENSE00001074448", "exon")
        self.assertTrue("ENSE00001074448" in ret_by_id)

    def test_exon_ids_from_gene(self):
        ret_by_id = get_exon_ids("ENSG00000157404", "gene")
        ret_by_name = get_exon_ids("KIT", "gene")
        self.assertTrue("ENSE00001074448" in ret_by_id)
        self.assertTrue("ENSE00001074448" in ret_by_name)
        self.assertEqual(ret_by_id, ret_by_name)

    def test_exon_ids_from_protein(self):
        ret_by_id = get_exon_ids("ENSP00000288135", "protein")
        self.assertTrue("ENSE00001074448" in ret_by_id)

    def test_exon_ids_from_transcript(self):
        ret_by_id = get_exon_ids("ENST00000288135", "transcript")
        ret_by_name = get_exon_ids("KIT-201", "transcript")
        self.assertTrue("ENSE00001074448" in ret_by_id)
        self.assertTrue("ENSE00001074448" in ret_by_name)
        self.assertEqual(ret_by_id, ret_by_name)

    def test_exon_from_cds(self):
        ret_by_id = get_exons("ENST00000288135", "cds")
        ret_by_name = get_exons("KIT-201", "cds")
        self.assertTrue("ENSE00001074448" in [i.exon_id for i in ret_by_id])
        self.assertTrue("ENSE00001074448" in [i.exon_id for i in ret_by_name])
        self.assertEqual(ret_by_id, ret_by_name)

    def test_exon_from_exon(self):
        ret_by_id = get_exons("ENSE00001074448", "exon")
        self.assertTrue("ENSE00001074448" in [i.exon_id for i in ret_by_id])

    def test_exon_from_gene(self):
        ret_by_id = get_exons("ENSG00000157404", "gene")
        ret_by_name = get_exons("KIT", "gene")
        self.assertTrue("ENSE00001074448" in [i.exon_id for i in ret_by_id])
        self.assertTrue("ENSE00001074448" in [i.exon_id for i in ret_by_name])
        self.assertEqual(ret_by_id, ret_by_name)

    def test_exon_from_protein(self):
        ret_by_id = get_exons("ENSP00000288135", "protein")
        self.assertTrue("ENSE00001074448" in [i.exon_id for i in ret_by_id])

    def test_exon_from_transcript(self):
        ret_by_id = get_exons("ENST00000288135", "transcript")
        ret_by_name = get_exons("KIT-201", "transcript")
        self.assertTrue("ENSE00001074448" in [i.exon_id for i in ret_by_id])
        self.assertTrue("ENSE00001074448" in [i.exon_id for i in ret_by_name])
        self.assertEqual(ret_by_id, ret_by_name)


class TestGene(unittest.TestCase):
    def test_gene_ids_from_cds(self):
        ret_by_id = get_gene_ids("ENST00000288135", "cds")
        ret_by_name = get_gene_ids("KIT-201", "cds")
        self.assertTrue("ENSG00000157404" in ret_by_id)
        self.assertTrue("ENSG00000157404" in ret_by_name)
        self.assertEqual(ret_by_id, ret_by_name)

    def test_gene_ids_from_exon(self):
        ret_by_id = get_gene_ids("ENSE00001074448", "exon")
        self.assertTrue("ENSG00000157404" in ret_by_id)

    def test_gene_ids_from_gene(self):
        ret_by_id = get_gene_ids("ENSG00000157404", "gene")
        ret_by_name = get_gene_ids("KIT", "gene")
        self.assertTrue("ENSG00000157404" in ret_by_id)
        self.assertTrue("ENSG00000157404" in ret_by_name)
        self.assertEqual(ret_by_id, ret_by_name)

    def test_gene_ids_from_protein(self):
        ret_by_id = get_gene_ids("ENSP00000288135", "protein")
        self.assertTrue("ENSG00000157404" in ret_by_id)

    def test_gene_ids_from_transcript(self):
        ret_by_id = get_gene_ids("ENST00000288135", "transcript")
        ret_by_name = get_gene_ids("KIT-201", "transcript")
        self.assertTrue("ENSG00000157404" in ret_by_id)
        self.assertTrue("ENSG00000157404" in ret_by_name)
        self.assertEqual(ret_by_id, ret_by_name)

    def test_gene_from_cds(self):
        ret_by_id = get_genes("ENST00000288135", "cds")
        ret_by_name = get_genes("KIT-201", "cds")
        self.assertTrue("ENSG00000157404" in [i.gene_id for i in ret_by_id])
        self.assertTrue("ENSG00000157404" in [i.gene_id for i in ret_by_name])
        self.assertEqual(ret_by_id, ret_by_name)

    def test_gene_from_exon(self):
        ret_by_id = get_genes("ENSE00001074448", "exon")
        self.assertTrue("ENSG00000157404" in [i.gene_id for i in ret_by_id])

    def test_gene_from_gene(self):
        ret_by_id = get_genes("ENSG00000157404", "gene")
        ret_by_name = get_genes("KIT", "gene")
        self.assertTrue("ENSG00000157404" in [i.gene_id for i in ret_by_id])
        self.assertTrue("ENSG00000157404" in [i.gene_id for i in ret_by_name])
        self.assertEqual(ret_by_id, ret_by_name)

    def test_gene_from_protein(self):
        ret_by_id = get_genes("ENSP00000288135", "protein")
        self.assertTrue("ENSG00000157404" in [i.gene_id for i in ret_by_id])

    def test_gene_from_transcript(self):
        ret_by_id = get_genes("ENST00000288135", "transcript")
        ret_by_name = get_genes("KIT-201", "transcript")
        self.assertTrue("ENSG00000157404" in [i.gene_id for i in ret_by_id])
        self.assertTrue("ENSG00000157404" in [i.gene_id for i in ret_by_name])
        self.assertEqual(ret_by_id, ret_by_name)


class TestProtein(unittest.TestCase):
    def test_protein_ids_from_cds(self):
        ret_by_id = get_protein_ids("ENST00000288135", "cds")
        ret_by_name = get_protein_ids("KIT-201", "transcript")
        self.assertTrue("ENSP00000288135" in ret_by_id)
        self.assertTrue("ENSP00000288135" in ret_by_name)
        self.assertEqual(ret_by_id, ret_by_name)

    def test_protein_ids_from_exon(self):
        ret_by_id = get_protein_ids("ENSE00001074448", "exon")
        self.assertTrue("ENSP00000288135" in ret_by_id)

    def test_protein_ids_from_gene(self):
        ret_by_id = get_protein_ids("ENSG00000157404", "gene")
        ret_by_name = get_protein_ids("KIT", "gene")
        self.assertTrue("ENSP00000288135" in ret_by_id)
        self.assertTrue("ENSP00000288135" in ret_by_name)
        self.assertEqual(ret_by_id, ret_by_name)

    def test_protein_ids_from_protein(self):
        ret_by_id = get_protein_ids("ENSP00000288135", "protein")
        self.assertTrue("ENSP00000288135" in ret_by_id)

    def test_protein_ids_from_transcript(self):
        ret_by_id = get_protein_ids("ENST00000288135", "transcript")
        ret_by_name = get_protein_ids("KIT-201", "transcript")
        self.assertTrue("ENSP00000288135" in ret_by_id)
        self.assertTrue("ENSP00000288135" in ret_by_name)
        self.assertEqual(ret_by_id, ret_by_name)


class TestTranscript(unittest.TestCase):
    def test_transcript_ids_from_cds(self):
        ret_by_id = get_transcript_ids("ENST00000288135", "cds")
        ret_by_name = get_transcript_ids("KIT-201", "cds")
        self.assertTrue("ENST00000288135" in ret_by_id)
        self.assertTrue("ENST00000288135" in ret_by_name)
        self.assertEqual(ret_by_id, ret_by_name)

    def test_transcript_ids_from_exon(self):
        ret_by_id = get_transcript_ids("ENSE00001074448", "exon")
        self.assertTrue("ENST00000288135" in ret_by_id)

    def test_transcript_ids_from_gene(self):
        ret_by_id = get_transcript_ids("ENSG00000157404", "gene")
        ret_by_name = get_transcript_ids("KIT", "gene")
        self.assertTrue("ENST00000288135" in ret_by_id)
        self.assertTrue("ENST00000288135" in ret_by_name)
        self.assertEqual(ret_by_id, ret_by_name)

    def test_transcript_ids_from_protein(self):
        ret_by_id = get_transcript_ids("ENSP00000288135", "protein")
        self.assertTrue("ENST00000288135" in ret_by_id)

    def test_transcript_ids_from_transcript(self):
        ret_by_id = get_transcript_ids("ENST00000288135", "transcript")
        ret_by_name = get_transcript_ids("KIT-201", "transcript")
        self.assertTrue("ENST00000288135" in ret_by_id)
        self.assertTrue("ENST00000288135" in ret_by_name)
        self.assertEqual(ret_by_id, ret_by_name)

    def test_transcript_from_cds(self):
        ret_by_id = get_transcripts("ENST00000288135", "cds")
        ret_by_name = get_transcripts("KIT-201", "cds")
        self.assertTrue("ENST00000288135" in [i.transcript_id for i in ret_by_id])
        self.assertTrue("ENST00000288135" in [i.transcript_id for i in ret_by_name])
        self.assertEqual(ret_by_id, ret_by_name)

    def test_transcript_from_exon(self):
        ret_by_id = get_transcripts("ENSE00003659301", "exon")
        self.assertTrue("ENST00000380152" in [i.transcript_id for i in ret_by_id])

    def test_transcript_from_gene(self):
        ret_by_id = get_transcripts("ENSG00000157404", "gene")
        ret_by_name = get_transcripts("KIT", "gene")
        self.assertTrue("ENST00000288135" in [i.transcript_id for i in ret_by_id])
        self.assertTrue("ENST00000288135" in [i.transcript_id for i in ret_by_name])
        self.assertEqual(ret_by_id, ret_by_name)

    def test_transcript_from_protein(self):
        ret_by_id = get_transcripts("ENSP00000288135", "protein")
        self.assertTrue("ENST00000288135" in [i.transcript_id for i in ret_by_id])

    def test_transcript_from_transcript(self):
        ret_by_id = get_transcripts("ENST00000288135", "transcript")
        ret_by_name = get_transcripts("KIT-201", "transcript")
        self.assertTrue("ENST00000288135" in [i.transcript_id for i in ret_by_id])
        self.assertTrue("ENST00000288135" in [i.transcript_id for i in ret_by_name])
        self.assertEqual(ret_by_id, ret_by_name)
