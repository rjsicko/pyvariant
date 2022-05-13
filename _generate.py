#!/usr/bin/env python
import sys
from itertools import product

data1 = [
    ("contig_id", "seqname"),
    ("exon_id", "exon_id"),
    ("gene_id", "gene_id"),
    ("gene_name", "gene_name"),
    ("protein_id", "protein_id"),
    ("transcript_id", "transcript_id"),
    ("transcript_name", "transcript_name"),
]
data2 = [
    "contig_id",
    "exon_id",
    "gene_id",
    "gene_name",
    "protein_id",
    "transcript_id",
    "transcript_name",
]

template1 = """\
def {a}s_of_{c}(self, value: Union[List[str], str]) -> List[str]:
    return self._query(value, "{c}", "{b}")
"""

template2 = """\
def test_{a}s_of_{c}(ensembl100):
    result = ensembl100.{a}s_of_{c}({d})
    assert isinstance(result, list)
    assert {e} in result
"""

for i, c in product(data1, data2):
    a, b = i
    print(template1.format(a=a, b=b, c=c), file=sys.stdout)
    print(template2.format(a=a, b=b, c=c, d=c.upper(), e=a.upper()), file=sys.stderr)
