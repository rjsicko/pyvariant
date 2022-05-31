#!/usr/bin/env python
import os
from ensembl_map.cache import Cache

path = os.path.join(os.path.dirname(__file__), "..", "tests", "cache")
if not os.path.exists(path):
    os.mkdir(path)

Cache("homo_sapiens", "GRCh38", 100, path).download_all()
Cache("homo_sapiens", "GRCh37", 69, path).download_all()
