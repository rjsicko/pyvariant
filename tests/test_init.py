from pyensembl import EnsemblRelease

import coordinate_mapper


def test_pyensembl():
    assert isinstance(coordinate_mapper.pyensembl, EnsemblRelease)
