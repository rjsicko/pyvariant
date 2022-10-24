from Bio.Seq import Seq


def reverse_complement(sequence: str) -> str:
    """Return the reverse complement of a given nucleotide sequence."""
    return str(Seq(sequence).reverse_complement())


def strip_version(key: str) -> str:
    """Strip the version number from a transcript symbol.

    Examples:
        >>> _strip_version('NM_000546.5')
        'NM_000546'
        >>> _strip_version('ENST00000357191.1')
        'ENST00000357191'
    """
    return key.rsplit(".", 1)[0]
