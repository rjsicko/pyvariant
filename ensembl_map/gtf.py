# Adapted from: https://gist.github.com/slowkow/8101481
import gzip
import re
from collections import defaultdict
from typing import Any, Dict, Iterator, List, Union

import pandas as pd

GTF_COLUMNS = {
    "seqname": str,
    "source": str,
    "feature": str,
    "start": int,
    "end": int,
    "score": int,
    "strand": str,
    "frame": int,
}
INFO_TYPES = {"exon_number": int}
NONE_VALUES = ["", ".", "NA"]
R_SEMICOLON = re.compile(r"\s*;\s*")
R_COMMA = re.compile(r"\s*,\s*")
R_KEYVALUE = re.compile(r"(\s+|\s*=\s*)")


def to_dataframe(filename: str) -> pd.DataFrame:
    """Parse a GTF file into a pandas.DataFrame."""
    # Each column is a list stored as a value in this dict
    result = defaultdict(list)

    for i, line in enumerate(_lines(filename)):
        for key in line.keys():
            # This key has not been seen yet, so set it to None for all previous lines
            if key not in result:
                result[key] = [None] * i

        # Ensure this row has some value for each column
        for key in result.keys():
            result[key].append(line.get(key, None))

    return pd.DataFrame(result)


def _lines(filename: str) -> Iterator[Dict[str, Any]]:
    """Iterate over the GTF file and parse each non-header line into a dict."""
    for line in _readlines(filename):
        if line.startswith("#"):
            continue
        else:
            yield _parse(line)


def _parse(line: str) -> Dict[str, Any]:
    """Parse a single GTF line into a dict."""
    result = {}

    fields = line.rstrip().split("\t")

    for n, col in enumerate(GTF_COLUMNS):
        result[col] = _get_value(fields[n], GTF_COLUMNS[col])

    # INFO field consists of "key1=value;key2=value;..."
    info = [x for x in re.split(R_SEMICOLON, fields[8]) if x.strip()]

    for n, info in enumerate(info, 1):
        # It should be key="value"
        try:
            key, _, value = re.split(R_KEYVALUE, info, 1)
        # But sometimes it is just "value"
        except ValueError:
            key = "INFO{}".format(n)
            value = info
        # Ignore the field if there is no value
        if value:
            result[key] = _get_value(value, INFO_TYPES.get(key, str))

    return result


def _readlines(filename: str):
    """Iterate over each line in an optionally gzipped file."""
    if filename.endswith(".gz"):
        openf = gzip.open
        mode = "rt"
    else:
        openf = open
        mode = "r"

    with openf(filename, mode) as fh:
        for line in fh:
            yield line


def _get_value(value: str, value_type: type) -> Union[List, str, None]:
    """Parse a value from a GTF INFO field."""
    parsed: Union[List, str, None]

    if not value:
        return None

    # Strip double and single quotes
    value = value.strip("\"'")

    # Some values should be treated as None
    if value in NONE_VALUES:
        return None

    # Return a list if the value has a comma
    if "," in value:
        parsed = []
        for i in re.split(R_COMMA, value):
            # Try to coerce each value into the expected type
            try:
                parsed.append(value_type(i))
            except ValueError:
                parsed.append(i)
    else:
        # Try to coerce the value into the expected type
        try:
            parsed = value_type(value)
        except ValueError:
            parsed = value

    return parsed
