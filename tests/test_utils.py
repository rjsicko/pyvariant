import gzip
import os.path

import pytest
from Bio.bgzf import BgzfWriter
from coordinate_mapper.utils import assert_valid_position, bgzip, is_bgzipped, is_ensembl_id


@pytest.fixture
def bgzip_file(tmpdir) -> str:
    path = tmpdir / "test.txt.gz"
    with BgzfWriter(path.strpath, compresslevel=9) as fh:
        fh.write(b"")

    return path.strpath


@pytest.fixture
def gzip_file(tmpdir) -> str:
    path = tmpdir / "test.txt.gz"
    with gzip.open(path.strpath, "wb") as fh:
        fh.write(b"")

    return path.strpath


@pytest.fixture
def uncompressed_file(tmpdir) -> str:
    path = tmpdir / "test.txt"
    with open(path.strpath, "w") as fh:
        fh.write("")

    return path.strpath


def test_assert_valid_postion_lt_one_1():
    with pytest.raises(ValueError):
        assert_valid_position(0)


def test_assert_valid_postion_lt_one_2():
    with pytest.raises(ValueError):
        assert_valid_position(0, 1)


def test_assert_valid_postion_lt_one_3():
    with pytest.raises(ValueError):
        assert_valid_position(1, 0)


def test_assert_valid_postion_compare():
    with pytest.raises(ValueError):
        assert_valid_position(100, 1)


def test_bgzip_gzip(gzip_file):
    bgzip(gzip_file)
    assert os.path.exists(gzip_file)


def test_bgzip_uncompressed(uncompressed_file):
    compressed_file = uncompressed_file + ".gz"
    bgzip(uncompressed_file)
    assert os.path.exists(compressed_file)


def test_is_bgzipped_false(gzip_file):
    assert is_bgzipped(gzip_file) is False


def test_is_bgzipped_true(bgzip_file):
    assert is_bgzipped(bgzip_file) is True


def test_is_ensembl_id_true():
    assert is_ensembl_id("ENSG00000000000") is True


def test_is_ensembl_id_true_1():
    assert is_ensembl_id("ENSG00000000000.0") is True


def test_is_ensembl_id_false():
    assert is_ensembl_id("blah") is False


def test_is_ensembl_id_false_2():
    assert is_ensembl_id("ENSG0") is False
