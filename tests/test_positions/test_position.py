import pickle

from pyvariant.positions import _Position


def test_pickle():
    position = _Position(
        contig_id="5", start=1282623, start_offset=0, end=1293313, end_offset=0, strand="-"
    )
    position2 = pickle.loads(pickle.dumps(position))
    assert position2.contig_id == "5"
    assert position2.start == 1282623
    assert position2.start_offset == 0
    assert position2.end == 1293313
    assert position2.end_offset == 0
    assert position2.strand == "-"
