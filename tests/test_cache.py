import os

import numpy
import pandas
import pytest

from ensembl_map.cache import Cache

from . import CACHE_DIR

NO_CACHE = not os.environ.get("CACHETEST")


@pytest.fixture
def ensembl100():
    return Cache(species="homo_sapiens", reference="GRCh38", release=100, cache_dir=CACHE_DIR)


@pytest.mark.skipif(NO_CACHE, reason="use 'CACHETEST' to run all cache tests")
def test_download_all(ensembl100):
    ensembl100.download_all(force=True)


def test_assign_protein_id(ensembl100):
    columns = ["feature", "transcript_id", "protein_id"]
    data = [
        ["CDS", "ENST00000123456", "ENSP00000654321"],
        ["CDS", "ENST00000412167", "ENSP00000390987"],
        ["stop_codon", "ENST00000412167", numpy.nan],
        ["cdna", "ENST00000412167", numpy.nan],
    ]
    in_df = pandas.DataFrame(data, columns=columns)
    out_df = ensembl100._assign_protein_id(in_df)
    assert out_df["protein_id"].to_list() == [
        "ENSP00000654321",
        "ENSP00000390987",
        "ENSP00000390987",
        "ENSP00000390987",
    ]
