import gzip
import os.path
import shutil
from tempfile import TemporaryDirectory
from typing import Callable, Dict, List

from Bio.bgzf import BgzfWriter, _bgzf_magic
from logzero import logger


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
        logger.info(f"Compressing {path} with bgzip (this may take some time)...")
        logger.debug(f"Writing compressed output to temporary file '{tempfile}'")
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


def txt_to_list(path: str, label: str) -> List[str]:
    """Parse a text file into a list of unique strings."""
    result = set()

    if path:
        logger.debug(f"Loading {label} from {path}")
        with open(path, "r") as fh:
            for line in fh:
                if line:
                    result.add(line.strip())
    else:
        logger.debug(f"No {label} to load")

    return sorted(result)


def tsv_to_dict(path: str, label: str) -> Dict[str, List[str]]:
    """Parse a TSV of one-to-one mappings."""
    result: Dict = {}

    if path:
        # load data from the TSV file, keep only unique values
        logger.debug(f"Loading {label} from {path}")
        with open(path, "r") as fh:
            for line in fh:
                if line:
                    alias, name = line.strip().split("\t")[:2]
                    result.setdefault(alias, set())
                    result[alias].add(name)

        # convert values to sorted lists
        for alias in result:
            result[alias] = sorted(result[alias])
    else:
        logger.debug(f"No {label} to load")

    return result
