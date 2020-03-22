class OutOfRangeError(Exception):
    def __init__(self, transcript, position, feature_type):
        self.transcript = transcript
        self.position = position
        self.feature_type = feature_type

    def __str__(self):
        return f"{self.position} is outside {self.feature_type} (1, {len(self.transcript.coding_sequence)})"
