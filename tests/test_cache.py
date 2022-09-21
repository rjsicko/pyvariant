import os

import pytest

from ensembl_map.cache import Cache

from . import CACHE_DIR

NO_CACHE = not os.environ.get("CACHETEST")


@pytest.fixture
def ensembl100_cache():
    return Cache(species="homo_sapiens", reference="GRCh38", release=100, cache_dir=CACHE_DIR)


@pytest.mark.skipif(NO_CACHE, reason="use 'CACHETEST' to run all cache tests")
def test_download_all(ensembl100_cache):
    ensembl100_cache.download_all(force=True)
