import os

import pytest
from appdirs import user_data_dir

from ensembl_map.cache import get_cache_dir
from ensembl_map.constants import ENSEMBL_MAP_CACHE_VAR


@pytest.fixture
def cache_env_var():
    os.environ[ENSEMBL_MAP_CACHE_VAR] = "TEST"
    yield
    del os.environ[ENSEMBL_MAP_CACHE_VAR]


def test_get_cache_dir():
    assert get_cache_dir() == user_data_dir()


def test_get_cache_dir_env_var(cache_env_var):
    assert get_cache_dir() == "TEST"
