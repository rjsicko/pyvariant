class CdsRangeError(Exception):
    def __init__(self, tobj, position):
        self.tobj = tobj
        self.position = position

    def __str__(self):
        return (
            f"{self.position} is not a CDS postion in {self.tobj.transcript_id}: "
            f"{self.tobj.coding_sequence_position_ranges}"
        )


class ExonRangeError(Exception):
    def __init__(self, tobj, position):
        self.tobj = tobj
        self.position = position

    def __str__(self):
        return (
            f"{self.position} is not an exon position in {self.tobj.transcript_id}: "
            f"{self.tobj.exons}"
        )
