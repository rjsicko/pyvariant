import pytest
from coordinate_mapper.exceptions import (
    CdsOutOfRange,
    ExonOutOfRange,
    GeneOutOfRange,
    OutOfRangeErrorBase,
    TranscriptOutOfRange,
)


def test_out_of_range_error_inherits_from_value_error(mocker):
    with pytest.raises(ValueError):
        raise OutOfRangeErrorBase()


def test_cds_out_of_range_error_inherits_from_value_error(mocker):
    with pytest.raises(ValueError):
        raise CdsOutOfRange(mocker.MagicMock(), -1)


def test_exon_out_of_range_error_inherits_from_value_error(mocker):
    with pytest.raises(ValueError):
        raise ExonOutOfRange(mocker.MagicMock(), -1)


def test_gene_out_of_range_error_inherits_from_value_error(mocker):
    with pytest.raises(ValueError):
        raise GeneOutOfRange(mocker.MagicMock(), -1)


def test_transcript_out_of_range_error_inherits_from_value_error(mocker):
    with pytest.raises(ValueError):
        raise TranscriptOutOfRange(mocker.MagicMock(), -1)
