class OutOfRangeError(Exception):
    def __init__(self, tobj, position, feature_type):
        self.tobj = tobj
        self.position = position
        self.feature_type = feature_type

    def __str__(self):
        return f"{self.position} is outside {self.feature_type} (1, {len(self.tobj.coding_sequence)})"
