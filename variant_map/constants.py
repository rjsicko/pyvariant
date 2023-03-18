"""Collection of constants used throughout the package."""
import os.path

# Package constants
NAME = "variant_map"

# Default species/release
DEFAULT_SPECIES = "homo_sapiens"
DEFAULT_ENSEMBL_RELEASE = 100

# Variable set the cache directory
CACHE_DIR_ENV = "VARIANT_MAP_CACHE"

# Literals used throughout the package
CONTIG_ID = "contig_id"
EXON_ID = "exon_id"
GENE_ID = "gene_id"
GENE_NAME = "gene_name"
PROTEIN_ID = "protein_id"
TRANSCRIPT_ID = "transcript_id"
TRANSCRIPT_NAME = "transcript_name"

# Position type literals
CDNA = "cdna"
CDS = "CDS"
DNA = "dna"
EXON = "exon"
PROTEIN = "protein"
RNA = "rna"
STOP_CODON = "stop_codon"

# Variant type literals
DELETION = "deletion"
DELINS = "delins"
DUPLICATION = "duplication"
FRAMESHIFT = "frameshift"
FUSION = "fusion"
INSERTION = "insertion"
SUBTITUTION = "substitution"

# Fallback FASTA file used if none are provided
EMPTY_FASTA = os.path.join(os.path.dirname(__file__), "data", "empty.fa")
