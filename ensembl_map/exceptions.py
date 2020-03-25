class ConvertError(ValueError):
    """Base error for `convert` function errors."""

    pass


class OutOfRangeErrorBase(ConvertError):
    """Base error for 'OutOfRange' errors."""

    def __init__(self, transcript, position):
        self.transcript = transcript
        self.position = position
        self.feature_type = None
        self.range = (None, None)

    def __str__(self):
        return f"{self.position} is outside {self.feature_type} {self.range}"


class CdsOutOfRange(OutOfRangeErrorBase):
    def __init__(self, transcript, position):
        self.transcript = transcript
        self.position = position
        self.feature_type = "CDS"
        self.range = (1, len(self.transcript.coding_sequence))


class ExonOutOfRange(OutOfRangeErrorBase):
    def __init__(self, transcript, position):
        self.transcript = transcript
        self.position = position
        self.feature_type = "exons"
        self.range = transcript.exon_intervals


class TranscriptOutOfRange(OutOfRangeErrorBase):
    def __init__(self, transcript, position):
        self.transcript = transcript
        self.position = position
        self.feature_type = "transcript"
        self.range = (1, len(self.transcript.sequence))
