import unittest
from unittest.mock import MagicMock

from ensembl_map.features import CDS, Exon, Gene, Protein, Transcript


class TestCds(unittest.TestCase):
    def setUp(self):
        self.obj = MagicMock()
        self.cds = CDS(self.obj, 1, 2)

    def test_biotype(self):
        self.assertEqual(self.cds.biotype, self.obj.biotype)

    def test_contig(self):
        self.assertEqual(self.cds.contig, self.obj.contig)

    def test_strand(self):
        self.assertEqual(self.cds.strand, self.obj.strand)

    def test_transcript(self):
        self.assertEqual(self.cds.transcript, self.obj)

    def test_transcript_id(self):
        self.assertEqual(self.cds.transcript_id, self.obj.transcript_id)

    def test_transcript_id_none(self):
        self.cds._pyensembl_obj = None  # causes AttributeError
        self.assertEqual(self.cds.transcript_id, None)

    def test_transcript_name(self):
        self.assertEqual(self.cds.transcript_name, self.obj.transcript_name)

    def test_transcript_name_none(self):
        self.cds._pyensembl_obj = None  # causes AttributeError
        self.assertEqual(self.cds.transcript_name, None)

    def test_sequence(self):
        self.assertEqual(self.cds.sequence, self.obj.coding_sequence.__getitem__())

    def test_to_tuple(self):
        self.assertEqual(self.cds.to_tuple(), (self.obj.transcript_id, 1, 2))


class TestExon(unittest.TestCase):
    def setUp(self):
        self.obj_t = MagicMock()
        self.obj_e = MagicMock()
        self.feature = Exon(self.obj_t, self.obj_e, 1, 2)

    def test_biotype(self):
        self.assertEqual(self.feature.biotype, self.obj_t.biotype)

    def test_contig(self):
        self.assertEqual(self.feature.contig, self.obj_t.contig)

    def test_strand(self):
        self.assertEqual(self.feature.strand, self.obj_t.strand)

    def test_exon(self):
        self.assertEqual(self.feature.exon, self.obj_e)

    def test_exon_id(self):
        self.assertEqual(self.feature.exon_id, self.obj_e.exon_id)

    def test_exon_id_none(self):
        self.feature._exon = None  # causes AttributeError
        self.assertEqual(self.feature.exon_id, None)

    # def test_exon_id_position(self):
    #     self.assertEqual(self.feature.exon_position, None)

    def test_transcript(self):
        self.assertEqual(self.feature.transcript, self.obj_t)

    def test_transcript_id(self):
        self.assertEqual(self.feature.transcript_id, self.obj_t.transcript_id)

    def test_transcript_id_none(self):
        self.feature._pyensembl_obj = None  # causes AttributeError
        self.assertEqual(self.feature.transcript_id, None)

    def test_transcript_name(self):
        self.assertEqual(self.feature.transcript_name, self.obj_t.transcript_name)

    def test_transcript_name_none(self):
        self.feature._pyensembl_obj = None  # causes AttributeError
        self.assertEqual(self.feature.transcript_name, None)

    def test_to_tuple(self):
        self.assertEqual(self.feature.to_tuple(), (self.obj_e.exon_id, 1, 2))


class TestGene(unittest.TestCase):
    def setUp(self):
        self.obj = MagicMock()
        self.feature = Gene(self.obj, 1, 2)

    def test_biotype(self):
        self.assertEqual(self.feature.biotype, self.obj.biotype)

    def test_contig(self):
        self.assertEqual(self.feature.contig, self.obj.contig)

    def test_strand(self):
        self.assertEqual(self.feature.strand, self.obj.strand)

    def test_gene(self):
        self.assertEqual(self.feature.gene, self.obj)

    def test_gene_id(self):
        self.assertEqual(self.feature.gene_id, self.obj.gene_id)

    def test_gene_id_none(self):
        self.feature._pyensembl_obj = None  # causes AttributeError
        self.assertEqual(self.feature.gene_id, None)

    def test_gene_name(self):
        self.assertEqual(self.feature.gene_name, self.obj.gene_name)

    def test_gene_name_none(self):
        self.feature._pyensembl_obj = None  # causes AttributeError
        self.assertEqual(self.feature.gene_name, None)

    def test_to_tuple(self):
        self.assertEqual(self.feature.to_tuple(), (self.obj.gene_id, 1, 2))


class TestProtein(unittest.TestCase):
    def setUp(self):
        self.obj = MagicMock()
        self.feature = Protein(self.obj, 1, 2)

    def test_biotype(self):
        self.assertEqual(self.feature.biotype, self.obj.biotype)

    def test_contig(self):
        self.assertEqual(self.feature.contig, self.obj.contig)

    def test_strand(self):
        self.assertEqual(self.feature.strand, self.obj.strand)

    def test_protein_id(self):
        self.assertEqual(self.feature.protein_id, self.obj.protein_id)

    def test_sequence(self):
        self.assertEqual(self.feature.sequence, self.obj.protein_sequence.__getitem__())

    def test_transcript(self):
        self.assertEqual(self.feature.transcript, self.obj)

    def test_to_tuple(self):
        self.assertEqual(self.feature.to_tuple(), (self.obj.protein_id, 1, 2))


class TestTranscript(unittest.TestCase):
    def setUp(self):
        self.obj = MagicMock()
        self.feature = Transcript(self.obj, 1, 2)

    def test_biotype(self):
        self.assertEqual(self.feature.biotype, self.obj.biotype)

    def test_contig(self):
        self.assertEqual(self.feature.contig, self.obj.contig)

    def test_strand(self):
        self.assertEqual(self.feature.strand, self.obj.strand)

    def test_transcript(self):
        self.assertEqual(self.feature.transcript, self.obj)

    def test_transcript_id(self):
        self.assertEqual(self.feature.transcript_id, self.obj.transcript_id)

    def test_transcript_id_none(self):
        self.feature._pyensembl_obj = None  # causes AttributeError
        self.assertEqual(self.feature.transcript_id, None)

    def test_transcript_name(self):
        self.assertEqual(self.feature.transcript_name, self.obj.transcript_name)

    def test_transcript_name_none(self):
        self.feature._pyensembl_obj = None  # causes AttributeError
        self.assertEqual(self.feature.transcript_name, None)

    def test_sequence(self):
        self.assertEqual(self.feature.sequence, self.obj.sequence.__getitem__())

    def test_to_tuple(self):
        self.assertEqual(self.feature.to_tuple(), (self.obj.transcript_id, 1, 2))
