import gzip
import os.path
import shutil
import sys
from tempfile import TemporaryDirectory
from typing import Callable, Dict, List

from Bio.bgzf import BgzfWriter, _bgzf_magic
from pyfaidx import Fasta

from .utils import strip_version

EMPTY_FASTA = os.path.join(os.path.dirname(__file__), "data", "empty.fa")


def bgzip(path: str) -> str:
    """Compress a file in BGZF format."""
    chunk_size = 4096  # chunk size is set to balance runtime and memory usage
    openfunc: Callable
    openmode: str
    remove_original: bool

    # if the input is gzip-compressed, recompress it with bgzip using the same file name
    if path.endswith(".gz") or path.endswith(".gzip"):
        openfunc = gzip.open
        openmode = "rb"
        output = path
        remove_original = False
    # if the input is not compressed, compress it with bgzip and append ".gz" to the file name
    else:
        openfunc = open
        openmode = "r"
        output = path + ".gz"
        remove_original = True

    # compress to a temporary file first to avoid overwritting the file as it's being read
    with TemporaryDirectory() as tempdir:
        tempfile = os.path.join(tempdir, os.path.basename(output))
        print(f"Compressing {path} with bgzip (this may take some time)...", file=sys.stderr)
        print(f"Writing compressed output to temporary file '{tempfile}'", file=sys.stderr)
        with openfunc(path, openmode) as inf:
            with BgzfWriter(tempfile, compresslevel=9) as outf:
                chunk = inf.read(chunk_size)
                while chunk:
                    outf.write(chunk)
                    chunk = inf.read(chunk_size)
        # move the temporary file to the same directory as the input
        shutil.move(tempfile, output)

    # remove the original to save space
    if remove_original:
        os.remove(path)

    return output


def is_bgzipped(path: str) -> bool:
    """Check if a file is compressed in BGZF format."""
    with open(path, "rb") as f:
        file_start = f.read(len(_bgzf_magic))
        return file_start == _bgzf_magic


def read_fasta(path: str) -> Fasta:
    """Parse a FASTA/FASTQ file into a 'pyfaidx.Fasta' object."""
    if not path:
        path = EMPTY_FASTA

    return Fasta(
        path,
        key_function=strip_version,
        as_raw=True,
        sequence_always_upper=True,
        build_index=False,
        rebuild=False,
    )


def tsv_to_dict(path: str) -> Dict[str, List[str]]:
    """Parse a TSV of one-to-one mappings."""
    result: Dict = {}

    if path:
        # load data from the TSV file, keep only unique values
        with open(path, "r") as fh:
            for line in fh:
                if line:
                    alias, name = line.strip().split("\t")[:2]
                    result.setdefault(alias, set())
                    result[alias].add(name)

        # convert values to sorted lists
        for alias in result:
            result[alias] = sorted(result[alias])

    return result


def txt_to_list(path: str) -> List[str]:
    """Parse a text file into a list of unique strings."""
    result = set()

    if path:
        with open(path, "r") as fh:
            for line in fh:
                if line:
                    result.add(line.strip())

    return sorted(result)
