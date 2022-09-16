import gzip
import os.path

import pytest
from Bio.bgzf import BgzfWriter

from ensembl_map.utils import bgzip, is_bgzipped


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

