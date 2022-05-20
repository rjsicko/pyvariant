#!/usr/bin/env python
"""
DESCRIPTION
===========
This module is for downloading data from BioMart.

AUTOMATIC DOWNLOAD
==================
This script can automatically download the required data from BioMart, given an Ensembl release.

MANUAL DOWNLOAD
===============
Alternatively, BioMart data can be manually downloaded instead. Specific columns are expected, so
here are instructions on how to export BioMart files with the required data.

Using Ensembl 100 as an example:

Step 1) go to Ensembl 100 archive: http://apr2020.archive.ensembl.org/index.html

Step 2) click on 'BioMart'

Step 3) select the correct dataset and attributes

    Dataset:
    * "- CHOOSE DATABASE -" -> "Ensembl Genes 100"
    * "- CHOOSE DATASET -" -> "Human Genes (GRCh38.p13)"

    Attributes:
        GENE:
            Ensembl:
            * TODO
        EXTERNAL:
            External References:
            * TODO

NOTE: If you see "Validation Error: Too many attributes selected for External References" when you
try to get results, do multiple queries instead, each with 1 or 2 External References, then merge
the exported files.

Step 4) click "Results" then next to "Export all results to" select "File" and "TSV" then check
"Unique results only" and finally click "Go" to download the results

OPTIONS
=======
"""
import os.path
import re

import requests

# paths to the BioMart query data, grouped by Ensembl release
XML_DATA = os.path.join(os.path.dirname(__file__), "data", "biomart_query.xml")

# defaults for command line arguments
DEFAULT_OUTPUT = "ens{release}_biomart_export.tsv"
DEFAULT_RELEASE = 100


def format_xml(xml_file: str) -> str:
    """Parse a BioMart query in XML format into a single line."""
    # read the XML data in from a file
    with open(xml_file, "r") as fh:
        data = fh.read()

    # collapse the XML into one line
    data = re.sub(r"\n\s{0,}", "", data)

    return data


def query_biomart(release: int, output: str):
    """Use the BioMart REST URL to export gene annotations."""
    xml = format_xml(XML_DATA)
    url = f"http://www.ensembl.org/biomart/martservice?query={xml}"
    response = requests.get(url)
    if response.ok:
        with open(output, "w") as fh:
            print(response.content.decode("utf-8"), file=fh)
    else:
        raise Exception(response.content)
