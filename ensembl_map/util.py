import re


def is_ensembl_id(feature):
    """The given string looks like and Ensembl ID."""
    return bool(re.match(r"ENS[A-Z]\d{8}", feature))

