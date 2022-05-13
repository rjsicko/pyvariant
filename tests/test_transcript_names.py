import pytest

from coordinate_mapper.transcript_names import TranscriptNames

from . import REFSEQ_FILE


@pytest.fixture()
def test_transcript_names():
    TranscriptNames.load(REFSEQ_FILE)


def test_data(test_transcript_names):
    assert "nm_00123" in TranscriptNames._refseq_to_ens
    assert TranscriptNames._refseq_to_ens["nm_00123"] == {"ENST00000487221"}


def test_to_ensembl_from_refseq(test_transcript_names):
    assert TranscriptNames.to_ensembl("NM_00123") == ["ENST00000487221"]


def test_to_ensembl_missing(test_transcript_names):
    assert TranscriptNames.to_ensembl("MISSING") == []
