import pytest
from constants import (
    CANONICAL_TRANSCRIPT,
    CONTIG_ALIAS,
    EXON_ALIAS,
    GENE_ALIAS,
    PROTEIN_ALIAS,
    TRANSCRIPT_ALIAS,
)

from ensembl_map.ensembl_release import EnsemblRelease


@pytest.fixture(scope="session")
def ensembl69():
    return EnsemblRelease(
        species="homo_sapiens",
        release=69,
        cache_dir="",
        canonical_transcript="",
        contig_alias="",
        exon_alias="",
        gene_alias="",
        protein_alias="",
        transcript_alias="",
    )


@pytest.fixture(scope="session")
def ensembl100():
    return EnsemblRelease(
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
