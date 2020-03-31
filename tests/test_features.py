import unittest
from unittest.mock import MagicMock

from ensembl_map.features import CDS, Exon, Gene, Protein, Transcript


class TestCds(unittest.TestCase):
    def setUp(self):
        self.obj = MagicMock()
        self.feature = CDS.load(self.obj, 1, 2)

    def test_biotype(self):
        self.assertEqual(self.feature.biotype, self.obj.biotype)

    def test_contig(self):
        self.assertEqual(self.feature.contig, self.obj.contig)

    def test_sequence(self):
        self.assertEqual(self.feature.sequence, self.obj.coding_sequence.__getitem__())

    def test_strand(self):
        self.assertEqual(self.feature.strand, self.obj.strand)

    def test_to_tuple(self):
        self.assertEqual(self.feature.to_tuple(), (self.obj.transcript_id, 1, 2))

    def test_transcript_id(self):
        self.assertEqual(self.feature.transcript_id, self.obj.transcript_id)

    def test_transcript_name(self):
        self.assertEqual(self.feature.transcript_name, self.obj.transcript_name)


class TestExon(unittest.TestCase):
    def setUp(self):
        self.obj_e = MagicMock()
        self.obj_e.exon_id = "test"
        self.obj_t = MagicMock()
        self.obj_t.exons = [self.obj_e]
        self.feature = Exon.load(self.obj_t, 1, 2, "test")

    def test_biotype(self):
        self.assertEqual(self.feature.biotype, self.obj_t.biotype)

    def test_contig(self):
        self.assertEqual(self.feature.contig, self.obj_t.contig)

    def test_exon_id(self):
        self.assertEqual(self.feature.exon_id, self.obj_e.exon_id)

    def test_exon_id_position(self):
        self.assertEqual(self.feature.index, 1)

    def test_strand(self):
        self.assertEqual(self.feature.strand, self.obj_t.strand)

    def test_to_tuple(self):
        self.assertEqual(self.feature.to_tuple(), (self.obj_e.exon_id, 1, 2))

    def test_transcript_id(self):
        self.assertEqual(self.feature.transcript_id, self.obj_t.transcript_id)

    def test_transcript_name(self):
        self.assertEqual(self.feature.transcript_name, self.obj_t.transcript_name)


class TestGene(unittest.TestCase):
    def setUp(self):
        self.obj = MagicMock()
        self.feature = Gene.load(self.obj, 1, 2)

    def test_biotype(self):
        self.assertEqual(self.feature.biotype, self.obj.gene.biotype)

    def test_contig(self):
        self.assertEqual(self.feature.contig, self.obj.gene.contig)

    def test_gene_id(self):
        self.assertEqual(self.feature.gene_id, self.obj.gene.gene_id)

    def test_gene_name(self):
        self.assertEqual(self.feature.gene_name, self.obj.gene.gene_name)

    def test_strand(self):
        self.assertEqual(self.feature.strand, self.obj.gene.strand)

    def test_to_tuple(self):
        self.assertEqual(self.feature.to_tuple(), (self.obj.gene.gene_id, 1, 2))


class TestProtein(unittest.TestCase):
    def setUp(self):
        self.obj = MagicMock()
        self.feature = Protein.load(self.obj, 1, 2)

    def test_biotype(self):
        self.assertEqual(self.feature.biotype, self.obj.biotype)

    def test_contig(self):
        self.assertEqual(self.feature.contig, self.obj.contig)

    def test_protein_id(self):
        self.assertEqual(self.feature.protein_id, self.obj.protein_id)

    def test_sequence(self):
        self.assertEqual(self.feature.sequence, self.obj.protein_sequence.__getitem__())

    def test_strand(self):
        self.assertEqual(self.feature.strand, self.obj.strand)

    def test_to_tuple(self):
        self.assertEqual(self.feature.to_tuple(), (self.obj.protein_id, 1, 2))


class TestTranscript(unittest.TestCase):
    def setUp(self):
        self.obj = MagicMock()
        self.feature = Transcript.load(self.obj, 1, 2)

    def test_biotype(self):
        self.assertEqual(self.feature.biotype, self.obj.biotype)

    def test_contig(self):
        self.assertEqual(self.feature.contig, self.obj.contig)

    def test_sequence(self):
        self.assertEqual(self.feature.sequence, self.obj.sequence.__getitem__())

    def test_strand(self):
        self.assertEqual(self.feature.strand, self.obj.strand)

    def test_to_tuple(self):
        self.assertEqual(self.feature.to_tuple(), (self.obj.transcript_id, 1, 2))

    def test_transcript_id(self):
        self.assertEqual(self.feature.transcript_id, self.obj.transcript_id)

    def test_transcript_name(self):
        self.assertEqual(self.feature.transcript_name, self.obj.transcript_name)
