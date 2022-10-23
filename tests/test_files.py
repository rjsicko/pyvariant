import gzip
import os.path

import pytest
from Bio.bgzf import BgzfWriter
from pyfaidx import Fasta

from ensembl_map.files import bgzip, is_bgzipped, read_fasta, tsv_to_dict, txt_to_list

from . import CANONICAL_TRANSCRIPT, CONTIG_ALIAS, TEST_DATA


@pytest.fixture
def bgzip_file(tmpdir) -> str:
    path = tmpdir / "test.txt.gz"
    with BgzfWriter(path.strpath, compresslevel=9) as fh:
        fh.write(b"test")

    return path.strpath


@pytest.fixture
def gzip_file(tmpdir) -> str:
    path = tmpdir / "test.txt.gz"
    with gzip.open(path.strpath, "wb") as fh:
        fh.write(b"test")

    return path.strpath


@pytest.fixture
def uncompressed_file(tmpdir) -> str:
    path = tmpdir / "test.txt"
    with open(path.strpath, "w") as fh:
        fh.write("test")

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


def test_read_fasta():
    fasta = os.path.join(TEST_DATA, "Homo_sapiens.GRCh38.cdna.all.fa.gz")
    result = read_fasta(fasta)
    assert isinstance(result, Fasta)
    assert result["ENST00000643777"]


def test_read_fasta_empty():
    result = read_fasta("")
    assert isinstance(result, Fasta)


def test_tsv_to_dict():
    result = tsv_to_dict(CONTIG_ALIAS)
    assert result == {"chr12": ["12"]}


def test_tsv_to_dict_empty():
    result = tsv_to_dict("")
    assert result == {}


def test_txt_to_list():
    result = txt_to_list(CANONICAL_TRANSCRIPT)
    assert result == ["ENST00000000233"]


def test_txt_to_list_empty():
    result = txt_to_list("")
    assert result == []
