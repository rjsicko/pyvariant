import ensembl_map


def test_init_ensembl_release():
    assert hasattr(ensembl_map, "EnsemblRelease")
