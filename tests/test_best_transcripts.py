import pytest

from coordinate_mapper.best_transcripts import BestTranscript

from . import BEST_FILE


@pytest.fixture(autouse=True)
def reset_after_each_test():
    BestTranscript.load(BEST_FILE)


def test_load():
    assert "TEST_GENE" in BestTranscript._data
    assert BestTranscript._data["TEST_GENE"] == "TEST_TRANSCRIPT"


def test_is_best():
    assert BestTranscript.is_best("TEST_TRANSCRIPT")
    assert not BestTranscript.is_best("TRANSCRIPT_DOES_NOT_EXIST")
