import os

import pytest
from appdirs import user_data_dir

from ensembl_map.cache import (
    get_cache_dir,
    normalize_release,
    normalize_species,
    reference_by_release,
)
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


def test_normalize_release():
    assert normalize_release(100) == 100
    assert normalize_release("100") == 100


def test_normalize_species():
    assert normalize_species("homo_sapiens") == "homo_sapiens"
    assert normalize_species("HOMO_SAPIENS") == "homo_sapiens"
    assert normalize_species("homo sapiens") == "homo_sapiens"
    assert normalize_species("homo-sapiens") == "homo_sapiens"


def test_reference_by_release():
    assert reference_by_release(54) == "GRCh36"
    assert reference_by_release(55) == "GRCh37"
    assert reference_by_release(75) == "GRCh37"
    assert reference_by_release(76) == "GRCh38"
    assert reference_by_release(100) == "GRCh38"
    with pytest.raises(ValueError):
        reference_by_release(0)
