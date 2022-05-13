from ensembl_map.data_cache.data_cache import DataCache
import pytest

CACHE_DIR = "/home/matt/Downloads/ensembl_map_data"  # DEBUG


@pytest.fixture
def data_cache_ens100():
    return DataCache(species="homo_sapiens", reference="GRCh38", release=100, cache_dir=CACHE_DIR)


def test_load_data(data_cache_ens100):
    data_cache_ens100.load_data()
