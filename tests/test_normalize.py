from ensembl_map.normalize import normalize_cds, normalize_protein, normalize_transcript


def test_normalize_cds():
    """
    Shifts position 3' from 2398 to 2410:
        2400 - CTATGTCCGGG - 2410
                        ^
        2400 - CTATGTCCGGG - 2410
                            ^
    """
    result = normalize_cds("ENST00000275493", 2408)
    assert result.sequence == "G"
    assert result.start == 2410
    assert result.end == 2410


def test_normalize_protein():
    """
    Shifts position 3' from 293 to 294:
        290 - CYANNTFGSAN - 300
                    ^
        290 - CYANNTFGSAN - 300
                    ^
    """
    result = normalize_protein("ENSP00000288135", 293)
    assert result.sequence == "N"
    assert result.start == 294
    assert result.end == 294


def test_normalize_transcript():
    """
    Collapses repeat from "AGAG" to "AG":
        68 - GCTCGCGGCGC - 78
                ^^^^
        68 - GCTCGCGGCGC - 78
                    ^^
    """
    result = normalize_transcript("ENST00000288135", 71, 74)
    assert result.sequence == "CG"
    assert result.start == 73
    assert result.end == 74
