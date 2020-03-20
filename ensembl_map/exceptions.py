class ConvertError(Exception):
    def __init__(self, from_type, to_type):
        self.from_type = from_type
        self.to_type = to_type

    def __str__(self):
        return f"Cannot convert {self.from_type} directly to {self.to_type}"


class OutOfRangeError(Exception):
    def __init__(self, tobj, position, feature_type):
        self.tobj = tobj
        self.position = position
        self.feature_type = feature_type

    def __str__(self):
        return f"{self.position} is outside {self.feature_type} (1, {len(self.tobj.coding_sequence)})"
