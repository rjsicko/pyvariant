import unittest
from unittest.mock import MagicMock

from ensembl_map.exceptions import (
    ConvertError,
    OutOfRangeErrorBase,
    CdsOutOfRange,
    ExonOutOfRange,
    TranscriptOutOfRange,
)


class TestExceptions(unittest.TestCase):
    def test_convert_error_inherits_from_value_error(self):
        with self.assertRaises(ValueError):
            raise ConvertError()

    def test_out_of_range_error_inherits_from_value_error(self):
        transcript = MagicMock()
        with self.assertRaises(ValueError):
            raise OutOfRangeErrorBase(transcript, -1)

    def test_cds_out_of_range_error_inherits_from_value_error(self):
        transcript = MagicMock()
        with self.assertRaises(ValueError):
            raise CdsOutOfRange(transcript, -1)

    def test_exon_out_of_range_error_inherits_from_value_error(self):
        transcript = MagicMock()
        with self.assertRaises(ValueError):
            raise ExonOutOfRange(transcript, -1)

    def test_transcript_out_of_range_error_inherits_from_value_error(self):
        transcript = MagicMock()
        with self.assertRaises(ValueError):
            raise TranscriptOutOfRange(transcript, -1)
