import pandas as pd
from constants import (
    CANONICAL_TRANSCRIPT,
    CONTIG_ALIAS,
    EXON_ALIAS,
    GENE_ALIAS,
    PROTEIN_ALIAS,
    TRANSCRIPT_ALIAS,
)
from pyfaidx import Fasta

from variant_map.ensembl_release import EnsemblRelease


def test_init():
    obj = EnsemblRelease(
        species="homo_sapiens",
        release=100,
        cache_dir="",
        canonical_transcript=CANONICAL_TRANSCRIPT,
        contig_alias=CONTIG_ALIAS,
        exon_alias=EXON_ALIAS,
        gene_alias=GENE_ALIAS,
        protein_alias=PROTEIN_ALIAS,
        transcript_alias=TRANSCRIPT_ALIAS,
    )
    assert isinstance(obj.df, pd.DataFrame)
    assert isinstance(obj.cds, list)
    assert not obj.cds
    assert isinstance(obj.dna, list)
    assert isinstance(obj.dna[0], Fasta)
    assert isinstance(obj.peptide, list)
    assert isinstance(obj.peptide[0], Fasta)
    assert isinstance(obj.rna, list)
    assert isinstance(obj.rna[0], Fasta)
    assert isinstance(obj.canonical_transcript, list)
    assert isinstance(obj.contig_alias, dict)
    assert isinstance(obj.exon_alias, dict)
    assert isinstance(obj.gene_alias, dict)
    assert isinstance(obj.protein_alias, dict)
    assert isinstance(obj.transcript_alias, dict)


def test_unique_instances():
    a = EnsemblRelease(species="homo_sapiens", release=100)
    b = EnsemblRelease(species="homo_sapiens", release=100)
    c = EnsemblRelease(species="homo_sapiens", release=69)
    assert a is b
    assert a is not c
