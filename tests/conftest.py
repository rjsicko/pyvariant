import pytest
from constants import (
    TEST_ENS100_CANONICAL_TRANSCRIPT,
    TEST_ENS100_CONTIG_ALIAS,
    TEST_ENS100_EXON_ALIAS,
    TEST_ENS100_GENE_ALIAS,
    TEST_ENS100_PROTEIN_ALIAS,
    TEST_ENS100_TRANSCRIPT_ALIAS,
)

from pyvariant.ensembl_release import EnsemblRelease


@pytest.fixture(scope="session")
def ensembl69():
    return EnsemblRelease(
        "homo_sapiens",
        69,
        cache_dir="",
        canonical_transcript=[],
        contig_alias={},
        exon_alias={},
        gene_alias={},
        protein_alias={},
        transcript_alias={},
    )


@pytest.fixture(scope="session")
def ensembl100():
    return EnsemblRelease(
        "homo_sapiens",
        100,
        cache_dir="",
        canonical_transcript=TEST_ENS100_CANONICAL_TRANSCRIPT,
        contig_alias=TEST_ENS100_CONTIG_ALIAS,
        exon_alias=TEST_ENS100_EXON_ALIAS,
        gene_alias=TEST_ENS100_GENE_ALIAS,
        protein_alias=TEST_ENS100_PROTEIN_ALIAS,
        transcript_alias=TEST_ENS100_TRANSCRIPT_ALIAS,
    )
