import os.path

CACHE_DIR = os.path.join(os.path.dirname(__file__), "cache")
TEST_DATA = os.path.join(CACHE_DIR, "homo_sapiens", "GRCh38", "100")
CANONICAL_TRANSCRIPT = os.path.join(TEST_DATA, "canonical_transcript.txt")
CONTIG_ALIAS = os.path.join(TEST_DATA, "contig_alias.tsv")
EXON_ALIAS = os.path.join(TEST_DATA, "exon_alias.tsv")
GENE_ALIAS = os.path.join(TEST_DATA, "gene_alias.tsv")
PROTEIN_ALIAS = os.path.join(TEST_DATA, "protein_alias.tsv")
TRANSCRIPT_ALIAS = os.path.join(TEST_DATA, "transcript_alias.tsv")
