#!/usr/bin/env python
import os.path

from constants import CACHE_DIR

from ensembl_map.cache import EnsemblCache

RESTRICT_GENES = [
    "AGBL1",
    "ALDOA",
    "BRCA2",
    "CHEK2",
    "KIT",
    "KMT2A",
    "MEN1",
    "MLL",
    "PTEN",
    "POU5F1",
    "SMARCA4",
    "SSX4",
    "TCF19",
    "TERT",
]


def write_canonical_transcript_file():
    with open(os.path.join(CACHE_DIR, "canonical_transcript.txt"), "w") as fh:
        print("ENST00000000233", file=fh)


def write_contig_alias_file():
    with open(os.path.join(CACHE_DIR, "contig_alias.tsv"), "w") as fh:
        print("chr4\t4", file=fh)


def write_exon_alias_file():
    with open(os.path.join(CACHE_DIR, "exon_alias.tsv"), "w"):
        pass


def write_gene_alias_file():
    with open(os.path.join(CACHE_DIR, "gene_alias.tsv"), "w"):
        pass


def write_protein_alias_file():
    with open(os.path.join(CACHE_DIR, "protein_alias.tsv"), "w"):
        pass


def write_transcript_alias_file():
    with open(os.path.join(CACHE_DIR, "transcript_alias.tsv"), "w") as fh:
        print("NM_000244\tENST00000394374", file=fh)
        print("NM_007194\tENST00000404276", file=fh)
        print("NM_001128849\tENST00000646693", file=fh)
        print("NM_000314\tENST00000371953", file=fh)
        print("NM_000314\tENST00000645317", file=fh)
        print("NM_001004491\tENST00000366480", file=fh)


if __name__ == "__main__":
    write_canonical_transcript_file()
    write_contig_alias_file()
    write_exon_alias_file()
    write_gene_alias_file()
    write_protein_alias_file()
    write_transcript_alias_file()
    for release in [69, 100]:
        cache = EnsemblCache("homo_sapiens", release, cache_dir=CACHE_DIR)
        cache.make(clean=False, restrict_genes=RESTRICT_GENES)
