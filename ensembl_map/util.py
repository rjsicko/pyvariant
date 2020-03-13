import re


def is_ensembl_id(feature):
    """The given string looks like and Ensembl ID."""
    return bool(re.match(r"ENS[A-Z]\d{8}", feature))


def in_ranges(number, ranges):
    """Assert that a number falls inside one or more ranges."""
    for i, j in ranges:
        if i <= number <= j:
            return True
    else:
        return False
