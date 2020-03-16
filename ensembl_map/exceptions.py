class OutsideCdsError(Exception):
    def __init__(self, tobj, position):
        self.tobj = tobj
        self.position = position

    def __str__(self):
        return f"{self.position} is outside CDS (1, {len(self.tobj.coding_sequence)})"


class OutsideTranscriptError(Exception):
    def __init__(self, tobj, position):
        self.tobj = tobj
        self.position = position

    def __str__(self):
        return f"{self.position} is outside transcript (1, {len(self.tobj.sequence)})"
