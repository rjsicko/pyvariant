import coordinate_mapper
from pyensembl import EnsemblRelease


def test_pyensembl():
    assert isinstance(coordinate_mapper.pyensembl, EnsemblRelease)
