import pytest

from variant_map.parser import parse, protein_letters_3to1

PARSE_CASES = [
    # cDNA/RNA deletion
    (
        "LRG_199t1:r.10del",
        {
            "breakpoint1": {
                "feature": "LRG_199t1",
                "start": 10,
                "start_offset": None,
                "start_seq": None,
                "end": 10,
                "end_offset": None,
                "end_seq": None,
                "strand": None,
                "refseq": None,
                "altseq": None,
                "prefix": "r",
                "position_type": "rna",
                "suffix": "del",
                "variant_type": "deletion",
            },
            "breakpoint2": {
                "feature": None,
                "start": None,
                "start_offset": None,
                "start_seq": None,
                "end": None,
                "end_offset": None,
                "end_seq": None,
                "strand": None,
                "refseq": None,
                "altseq": None,
                "prefix": None,
                "position_type": None,
                "suffix": None,
                "variant_type": None,
            },
        },
    ),
    # cDNA/RNA deletion (multiple bases)
    (
        "LRG_2t1:r.1034_1036del",
        {
            "breakpoint1": {
                "feature": "LRG_2t1",
                "start": 1034,
                "start_offset": None,
                "start_seq": None,
                "end": 1036,
                "end_offset": None,
                "end_seq": None,
                "strand": None,
                "refseq": None,
                "altseq": None,
                "prefix": "r",
                "position_type": "rna",
                "suffix": "del",
                "variant_type": "deletion",
            },
            "breakpoint2": {
                "feature": None,
                "start": None,
                "start_offset": None,
                "start_seq": None,
                "end": None,
                "end_offset": None,
                "end_seq": None,
                "strand": None,
                "refseq": None,
                "altseq": None,
                "prefix": None,
                "position_type": None,
                "suffix": None,
                "variant_type": None,
            },
        },
    ),
    # cDNA/RNA deletion-insertion
    (
        "LRG_2t1:c.775delinsga",
        {
            "breakpoint1": {
                "feature": "LRG_2t1",
                "start": 775,
                "start_offset": None,
                "start_seq": None,
                "end": 775,
                "end_offset": None,
                "end_seq": None,
                "strand": None,
                "refseq": None,
                "altseq": "GA",
                "prefix": "c",
                "position_type": "cdna",
                "suffix": "delins",
                "variant_type": "delins",
            },
            "breakpoint2": {
                "feature": None,
                "start": None,
                "start_offset": None,
                "start_seq": None,
                "end": None,
                "end_offset": None,
                "end_seq": None,
                "strand": None,
                "refseq": None,
                "altseq": None,
                "prefix": None,
                "position_type": None,
                "suffix": None,
                "variant_type": None,
            },
        },
    ),
    # cDNA/RNA deletion-insertion (multiple bases)
    (
        "LRG_2t1:r.775_777delinsc",
        {
            "breakpoint1": {
                "feature": "LRG_2t1",
                "start": 775,
                "start_offset": None,
                "start_seq": None,
                "end": 777,
                "end_offset": None,
                "end_seq": None,
                "strand": None,
                "refseq": None,
                "altseq": "C",
                "prefix": "r",
                "position_type": "rna",
                "suffix": "delins",
                "variant_type": "delins",
            },
            "breakpoint2": {
                "feature": None,
                "start": None,
                "start_offset": None,
                "start_seq": None,
                "end": None,
                "end_offset": None,
                "end_seq": None,
                "strand": None,
                "refseq": None,
                "altseq": None,
                "prefix": None,
                "position_type": None,
                "suffix": None,
                "variant_type": None,
            },
        },
    ),
    # cDNA/RNA duplication
    (
        "LRG_2t1:r.7dup",
        {
            "breakpoint1": {
                "feature": "LRG_2t1",
                "start": 7,
                "start_offset": None,
                "start_seq": None,
                "end": 7,
                "end_offset": None,
                "end_seq": None,
                "strand": None,
                "refseq": None,
                "altseq": None,
                "prefix": "r",
                "position_type": "rna",
                "suffix": "dup",
                "variant_type": "duplication",
            },
            "breakpoint2": {
                "feature": None,
                "start": None,
                "start_offset": None,
                "start_seq": None,
                "end": None,
                "end_offset": None,
                "end_seq": None,
                "strand": None,
                "refseq": None,
                "altseq": None,
                "prefix": None,
                "position_type": None,
                "suffix": None,
                "variant_type": None,
            },
        },
    ),
    # cDNA/RNA duplication (multiple bases)
    (
        "LRG_2t1:r.6_8dup",
        {
            "breakpoint1": {
                "feature": "LRG_2t1",
                "start": 6,
                "start_offset": None,
                "start_seq": None,
                "end": 8,
                "end_offset": None,
                "end_seq": None,
                "strand": None,
                "refseq": None,
                "altseq": None,
                "prefix": "r",
                "position_type": "rna",
                "suffix": "dup",
                "variant_type": "duplication",
            },
            "breakpoint2": {
                "feature": None,
                "start": None,
                "start_offset": None,
                "start_seq": None,
                "end": None,
                "end_offset": None,
                "end_seq": None,
                "strand": None,
                "refseq": None,
                "altseq": None,
                "prefix": None,
                "position_type": None,
                "suffix": None,
                "variant_type": None,
            },
        },
    ),
    # cDNA/RNA insertion
    (
        "LRG_199t1:r.426_427insa",
        {
            "breakpoint1": {
                "feature": "LRG_199t1",
                "start": 426,
                "start_offset": None,
                "start_seq": None,
                "end": 427,
                "end_offset": None,
                "end_seq": None,
                "strand": None,
                "refseq": None,
                "altseq": "A",
                "prefix": "r",
                "position_type": "rna",
                "suffix": "ins",
                "variant_type": "insertion",
            },
            "breakpoint2": {
                "feature": None,
                "start": None,
                "start_offset": None,
                "start_seq": None,
                "end": None,
                "end_offset": None,
                "end_seq": None,
                "strand": None,
                "refseq": None,
                "altseq": None,
                "prefix": None,
                "position_type": None,
                "suffix": None,
                "variant_type": None,
            },
        },
    ),
    # cDNA/RNA substitution
    (
        "NM_004006.3:r.76a>c",
        {
            "breakpoint1": {
                "feature": "NM_004006.3",
                "start": 76,
                "start_offset": None,
                "start_seq": None,
                "end": 76,
                "end_offset": None,
                "end_seq": None,
                "strand": None,
                "refseq": "A",
                "altseq": "C",
                "prefix": "r",
                "position_type": "rna",
                "suffix": ">",
                "variant_type": "substitution",
            },
            "breakpoint2": {
                "feature": None,
                "start": None,
                "start_offset": None,
                "start_seq": None,
                "end": None,
                "end_offset": None,
                "end_seq": None,
                "strand": None,
                "refseq": None,
                "altseq": None,
                "prefix": None,
                "position_type": None,
                "suffix": None,
                "variant_type": None,
            },
        },
    ),
    # cDNA/RNA substitution (intron)
    (
        "NM_004006.3:r.76+4a>c",
        {
            "breakpoint1": {
                "feature": "NM_004006.3",
                "start": 76,
                "start_offset": 4,
                "start_seq": None,
                "end": 76,
                "end_offset": 4,
                "end_seq": None,
                "strand": None,
                "refseq": "A",
                "altseq": "C",
                "prefix": "r",
                "position_type": "rna",
                "suffix": ">",
                "variant_type": "substitution",
            },
            "breakpoint2": {
                "feature": None,
                "start": None,
                "start_offset": None,
                "start_seq": None,
                "end": None,
                "end_offset": None,
                "end_seq": None,
                "strand": None,
                "refseq": None,
                "altseq": None,
                "prefix": None,
                "position_type": None,
                "suffix": None,
                "variant_type": None,
            },
        },
    ),
    # cDNA/RNA substitution (promoter)
    (
        "TERT:c.-124C>T",
        {
            "breakpoint1": {
                "feature": "TERT",
                "start": 1,
                "start_offset": -124,
                "start_seq": None,
                "end": 1,
                "end_offset": -124,
                "end_seq": None,
                "strand": None,
                "refseq": "C",
                "altseq": "T",
                "prefix": "c",
                "position_type": "cdna",
                "suffix": ">",
                "variant_type": "substitution",
            },
            "breakpoint2": {
                "feature": None,
                "start": None,
                "start_offset": None,
                "start_seq": None,
                "end": None,
                "end_offset": None,
                "end_seq": None,
                "strand": None,
                "refseq": None,
                "altseq": None,
                "prefix": None,
                "position_type": None,
                "suffix": None,
                "variant_type": None,
            },
        },
    ),
    # DNA deletion
    (
        "NC_000023.11:g.33344591del",
        {
            "breakpoint1": {
                "feature": "NC_000023.11",
                "start": 33344591,
                "start_offset": None,
                "start_seq": None,
                "end": 33344591,
                "end_offset": None,
                "end_seq": None,
                "strand": None,
                "refseq": None,
                "altseq": None,
                "prefix": "g",
                "position_type": "dna",
                "suffix": "del",
                "variant_type": "deletion",
            },
            "breakpoint2": {
                "feature": None,
                "start": None,
                "start_offset": None,
                "start_seq": None,
                "end": None,
                "end_offset": None,
                "end_seq": None,
                "strand": None,
                "refseq": None,
                "altseq": None,
                "prefix": None,
                "position_type": None,
                "suffix": None,
                "variant_type": None,
            },
        },
    ),
    # DNA deletion (multiple bases)
    (
        "NC_000023.11:g.33344590_33344592del",
        {
            "breakpoint1": {
                "feature": "NC_000023.11",
                "start": 33344590,
                "start_offset": None,
                "start_seq": None,
                "end": 33344592,
                "end_offset": None,
                "end_seq": None,
                "strand": None,
                "refseq": None,
                "altseq": None,
                "prefix": "g",
                "position_type": "dna",
                "suffix": "del",
                "variant_type": "deletion",
            },
            "breakpoint2": {
                "feature": None,
                "start": None,
                "start_offset": None,
                "start_seq": None,
                "end": None,
                "end_offset": None,
                "end_seq": None,
                "strand": None,
                "refseq": None,
                "altseq": None,
                "prefix": None,
                "position_type": None,
                "suffix": None,
                "variant_type": None,
            },
        },
    ),
    # DNA deletion-insertion
    (
        "NC_000023.11:g.32386323delinsGA",
        {
            "breakpoint1": {
                "feature": "NC_000023.11",
                "start": 32386323,
                "start_offset": None,
                "start_seq": None,
                "end": 32386323,
                "end_offset": None,
                "end_seq": None,
                "strand": None,
                "refseq": None,
                "altseq": "GA",
                "prefix": "g",
                "position_type": "dna",
                "suffix": "delins",
                "variant_type": "delins",
            },
            "breakpoint2": {
                "feature": None,
                "start": None,
                "start_offset": None,
                "start_seq": None,
                "end": None,
                "end_offset": None,
                "end_seq": None,
                "strand": None,
                "refseq": None,
                "altseq": None,
                "prefix": None,
                "position_type": None,
                "suffix": None,
                "variant_type": None,
            },
        },
    ),
    # DNA deletion-insertion (multiple bases)
    (
        "NC_000023.11:g.32386323_32386325delinsC",
        {
            "breakpoint1": {
                "feature": "NC_000023.11",
                "start": 32386323,
                "start_offset": None,
                "start_seq": None,
                "end": 32386325,
                "end_offset": None,
                "end_seq": None,
                "strand": None,
                "refseq": None,
                "altseq": "C",
                "prefix": "g",
                "position_type": "dna",
                "suffix": "delins",
                "variant_type": "delins",
            },
            "breakpoint2": {
                "feature": None,
                "start": None,
                "start_offset": None,
                "start_seq": None,
                "end": None,
                "end_offset": None,
                "end_seq": None,
                "strand": None,
                "refseq": None,
                "altseq": None,
                "prefix": None,
                "position_type": None,
                "suffix": None,
                "variant_type": None,
            },
        },
    ),
    # DNA duplication
    (
        "NC_000023.11:g.32343183dup",
        {
            "breakpoint1": {
                "feature": "NC_000023.11",
                "start": 32343183,
                "start_offset": None,
                "start_seq": None,
                "end": 32343183,
                "end_offset": None,
                "end_seq": None,
                "strand": None,
                "refseq": None,
                "altseq": None,
                "prefix": "g",
                "position_type": "dna",
                "suffix": "dup",
                "variant_type": "duplication",
            },
            "breakpoint2": {
                "feature": None,
                "start": None,
                "start_offset": None,
                "start_seq": None,
                "end": None,
                "end_offset": None,
                "end_seq": None,
                "strand": None,
                "refseq": None,
                "altseq": None,
                "prefix": None,
                "position_type": None,
                "suffix": None,
                "variant_type": None,
            },
        },
    ),
    # DNA duplication (multiple bases)
    (
        "NC_000023.11:g.31060227_33274278dup",
        {
            "breakpoint1": {
                "feature": "NC_000023.11",
                "start": 31060227,
                "start_offset": None,
                "start_seq": None,
                "end": 33274278,
                "end_offset": None,
                "end_seq": None,
                "strand": None,
                "refseq": None,
                "altseq": None,
                "prefix": "g",
                "position_type": "dna",
                "suffix": "dup",
                "variant_type": "duplication",
            },
            "breakpoint2": {
                "feature": None,
                "start": None,
                "start_offset": None,
                "start_seq": None,
                "end": None,
                "end_offset": None,
                "end_seq": None,
                "strand": None,
                "refseq": None,
                "altseq": None,
                "prefix": None,
                "position_type": None,
                "suffix": None,
                "variant_type": None,
            },
        },
    ),
    # DNA insertion
    (
        "NC_000023.10:g.32867861_32867862insT",
        {
            "breakpoint1": {
                "feature": "NC_000023.10",
                "start": 32867861,
                "start_offset": None,
                "start_seq": None,
                "end": 32867862,
                "end_offset": None,
                "end_seq": None,
                "strand": None,
                "refseq": None,
                "altseq": "T",
                "prefix": "g",
                "position_type": "dna",
                "suffix": "ins",
                "variant_type": "insertion",
            },
            "breakpoint2": {
                "feature": None,
                "start": None,
                "start_offset": None,
                "start_seq": None,
                "end": None,
                "end_offset": None,
                "end_seq": None,
                "strand": None,
                "refseq": None,
                "altseq": None,
                "prefix": None,
                "position_type": None,
                "suffix": None,
                "variant_type": None,
            },
        },
    ),
    # DNA substitution
    (
        "NC_000023.10:g.33038255C>A",
        {
            "breakpoint1": {
                "feature": "NC_000023.10",
                "start": 33038255,
                "start_offset": None,
                "start_seq": None,
                "end": 33038255,
                "end_offset": None,
                "end_seq": None,
                "strand": None,
                "refseq": "C",
                "altseq": "A",
                "prefix": "g",
                "position_type": "dna",
                "suffix": ">",
                "variant_type": "substitution",
            },
            "breakpoint2": {
                "feature": None,
                "start": None,
                "start_offset": None,
                "start_seq": None,
                "end": None,
                "end_offset": None,
                "end_seq": None,
                "strand": None,
                "refseq": None,
                "altseq": None,
                "prefix": None,
                "position_type": None,
                "suffix": None,
                "variant_type": None,
            },
        },
    ),
    # Protein deletion
    (
        "NP_003997.2:p.Val7del",
        {
            "breakpoint1": {
                "feature": "NP_003997.2",
                "start": 7,
                "start_offset": None,
                "start_seq": "V",
                "end": 7,
                "end_offset": None,
                "end_seq": None,
                "strand": None,
                "refseq": None,
                "altseq": None,
                "prefix": "p",
                "position_type": "protein",
                "suffix": "del",
                "variant_type": "deletion",
            },
            "breakpoint2": {
                "feature": None,
                "start": None,
                "start_offset": None,
                "start_seq": None,
                "end": None,
                "end_offset": None,
                "end_seq": None,
                "strand": None,
                "refseq": None,
                "altseq": None,
                "prefix": None,
                "position_type": None,
                "suffix": None,
                "variant_type": None,
            },
        },
    ),
    # Protein deletion (multiple bases)
    (
        "NP_003997.2:p.Lys23_Val25del",
        {
            "breakpoint1": {
                "feature": "NP_003997.2",
                "start": 23,
                "start_offset": None,
                "start_seq": "K",
                "end": 25,
                "end_offset": None,
                "end_seq": "V",
                "strand": None,
                "refseq": None,
                "altseq": None,
                "prefix": "p",
                "position_type": "protein",
                "suffix": "del",
                "variant_type": "deletion",
            },
            "breakpoint2": {
                "feature": None,
                "start": None,
                "start_offset": None,
                "start_seq": None,
                "end": None,
                "end_offset": None,
                "end_seq": None,
                "strand": None,
                "refseq": None,
                "altseq": None,
                "prefix": None,
                "position_type": None,
                "suffix": None,
                "variant_type": None,
            },
        },
    ),
    # Protein deletion-insertion
    (
        "NP_003997.2:p.Cys28delinsTrpVal",
        {
            "breakpoint1": {
                "feature": "NP_003997.2",
                "start": 28,
                "start_offset": None,
                "start_seq": "C",
                "end": 28,
                "end_offset": None,
                "end_seq": None,
                "strand": None,
                "refseq": None,
                "altseq": "WV",
                "prefix": "p",
                "position_type": "protein",
                "suffix": "delins",
                "variant_type": "delins",
            },
            "breakpoint2": {
                "feature": None,
                "start": None,
                "start_offset": None,
                "start_seq": None,
                "end": None,
                "end_offset": None,
                "end_seq": None,
                "strand": None,
                "refseq": None,
                "altseq": None,
                "prefix": None,
                "position_type": None,
                "suffix": None,
                "variant_type": None,
            },
        },
    ),
    # Protein deletion-insertion (multiple bases)
    (
        "NP_003997.2:p.Cys28_Lys29delinsTrp",
        {
            "breakpoint1": {
                "feature": "NP_003997.2",
                "start": 28,
                "start_offset": None,
                "start_seq": "C",
                "end": 29,
                "end_offset": None,
                "end_seq": "K",
                "strand": None,
                "refseq": None,
                "altseq": "W",
                "prefix": "p",
                "position_type": "protein",
                "suffix": "delins",
                "variant_type": "delins",
            },
            "breakpoint2": {
                "feature": None,
                "start": None,
                "start_offset": None,
                "start_seq": None,
                "end": None,
                "end_offset": None,
                "end_seq": None,
                "strand": None,
                "refseq": None,
                "altseq": None,
                "prefix": None,
                "position_type": None,
                "suffix": None,
                "variant_type": None,
            },
        },
    ),
    # Protein duplication
    (
        "NP_003997.2:p.Val7dup",
        {
            "breakpoint1": {
                "feature": "NP_003997.2",
                "start": 7,
                "start_offset": None,
                "start_seq": "V",
                "end": 7,
                "end_offset": None,
                "end_seq": None,
                "strand": None,
                "refseq": None,
                "altseq": None,
                "prefix": "p",
                "position_type": "protein",
                "suffix": "dup",
                "variant_type": "duplication",
            },
            "breakpoint2": {
                "feature": None,
                "start": None,
                "start_offset": None,
                "start_seq": None,
                "end": None,
                "end_offset": None,
                "end_seq": None,
                "strand": None,
                "refseq": None,
                "altseq": None,
                "prefix": None,
                "position_type": None,
                "suffix": None,
                "variant_type": None,
            },
        },
    ),
    # Protein duplication (multiple bases)
    (
        "NP_003997.2:p.Lys23_Val25dup",
        {
            "breakpoint1": {
                "feature": "NP_003997.2",
                "start": 23,
                "start_offset": None,
                "start_seq": "K",
                "end": 25,
                "end_offset": None,
                "end_seq": "V",
                "strand": None,
                "refseq": None,
                "altseq": None,
                "prefix": "p",
                "position_type": "protein",
                "suffix": "dup",
                "variant_type": "duplication",
            },
            "breakpoint2": {
                "feature": None,
                "start": None,
                "start_offset": None,
                "start_seq": None,
                "end": None,
                "end_offset": None,
                "end_seq": None,
                "strand": None,
                "refseq": None,
                "altseq": None,
                "prefix": None,
                "position_type": None,
                "suffix": None,
                "variant_type": None,
            },
        },
    ),
    # Protein frameshift
    (
        "NP_003997.2:p.Arg97fs",
        {
            "breakpoint1": {
                "feature": "NP_003997.2",
                "start": 97,
                "start_offset": None,
                "start_seq": "R",
                "end": 97,
                "end_offset": None,
                "end_seq": None,
                "strand": None,
                "refseq": None,
                "altseq": None,
                "prefix": "p",
                "position_type": "protein",
                "suffix": "fs",
                "variant_type": "frameshift",
            },
            "breakpoint2": {
                "feature": None,
                "start": None,
                "start_offset": None,
                "start_seq": None,
                "end": None,
                "end_offset": None,
                "end_seq": None,
                "strand": None,
                "refseq": None,
                "altseq": None,
                "prefix": None,
                "position_type": None,
                "suffix": None,
                "variant_type": None,
            },
        },
    ),
    # Protein insertion
    (
        "NP_003997.2:p.His4_Gln5insAla",
        {
            "breakpoint1": {
                "feature": "NP_003997.2",
                "start": 4,
                "start_offset": None,
                "start_seq": "H",
                "end": 5,
                "end_offset": None,
                "end_seq": "Q",
                "strand": None,
                "refseq": None,
                "altseq": "A",
                "prefix": "p",
                "position_type": "protein",
                "suffix": "ins",
                "variant_type": "insertion",
            },
            "breakpoint2": {
                "feature": None,
                "start": None,
                "start_offset": None,
                "start_seq": None,
                "end": None,
                "end_offset": None,
                "end_seq": None,
                "strand": None,
                "refseq": None,
                "altseq": None,
                "prefix": None,
                "position_type": None,
                "suffix": None,
                "variant_type": None,
            },
        },
    ),
    # Protein position
    (
        "KIT:p.42",
        {
            "breakpoint1": {
                "feature": "KIT",
                "start": 42,
                "start_offset": None,
                "start_seq": None,
                "end": 42,
                "end_offset": None,
                "end_seq": None,
                "strand": None,
                "refseq": None,
                "altseq": None,
                "prefix": "p",
                "position_type": "protein",
                "suffix": None,
                "variant_type": None,
            },
            "breakpoint2": {
                "feature": None,
                "start": None,
                "start_offset": None,
                "start_seq": None,
                "end": None,
                "end_offset": None,
                "end_seq": None,
                "strand": None,
                "refseq": None,
                "altseq": None,
                "prefix": None,
                "position_type": None,
                "suffix": None,
                "variant_type": None,
            },
        },
    ),
    # Protein substitution
    (
        "LRG_199p1:p.Trp24Cys",
        {
            "breakpoint1": {
                "feature": "LRG_199p1",
                "start": 24,
                "start_offset": None,
                "start_seq": None,
                "end": 24,
                "end_offset": None,
                "end_seq": None,
                "strand": None,
                "refseq": "W",
                "altseq": "C",
                "prefix": "p",
                "position_type": "protein",
                "suffix": None,
                "variant_type": "substitution",
            },
            "breakpoint2": {
                "feature": None,
                "start": None,
                "start_offset": None,
                "start_seq": None,
                "end": None,
                "end_offset": None,
                "end_seq": None,
                "strand": None,
                "refseq": None,
                "altseq": None,
                "prefix": None,
                "position_type": None,
                "suffix": None,
                "variant_type": None,
            },
        },
    ),
    # DNA fusion (non-HGVS)
    (
        "(chrX,chrX):fusion(g.77589872,g.86714548)",
        {
            "breakpoint1": {
                "feature": "chrX",
                "start": 77589872,
                "start_offset": None,
                "start_seq": None,
                "end": 77589872,
                "end_offset": None,
                "end_seq": None,
                "strand": None,
                "refseq": None,
                "altseq": None,
                "prefix": "g",
                "position_type": "dna",
                "suffix": None,
                "variant_type": "fusion",
            },
            "breakpoint2": {
                "feature": "chrX",
                "start": 86714548,
                "start_offset": None,
                "start_seq": None,
                "end": 86714548,
                "end_offset": None,
                "end_seq": None,
                "strand": None,
                "refseq": None,
                "altseq": None,
                "prefix": "g",
                "position_type": "dna",
                "suffix": None,
                "variant_type": "fusion",
            },
        },
    ),
    # Exon fusion (non-HGVS)
    (
        "(ABCA5,PPP4R1L):fusion(e.1,e.16)",
        {
            "breakpoint1": {
                "feature": "ABCA5",
                "start": 1,
                "start_offset": None,
                "start_seq": None,
                "end": 1,
                "end_offset": None,
                "end_seq": None,
                "strand": None,
                "refseq": None,
                "altseq": None,
                "prefix": "e",
                "position_type": "exon",
                "suffix": None,
                "variant_type": "fusion",
            },
            "breakpoint2": {
                "feature": "PPP4R1L",
                "start": 16,
                "start_offset": None,
                "start_seq": None,
                "end": 16,
                "end_offset": None,
                "end_seq": None,
                "strand": None,
                "refseq": None,
                "altseq": None,
                "prefix": "e",
                "position_type": "exon",
                "suffix": None,
                "variant_type": "fusion",
            },
        },
    ),
    # RNA fusion
    (
        "NM_152263.2:r.-115_775::NM_002609.3:r.1580_1924",
        {
            "breakpoint1": {
                "feature": "NM_152263.2",
                "start": 1,
                "start_offset": -115,
                "start_seq": None,
                "end": 775,
                "end_offset": None,
                "end_seq": None,
                "strand": None,
                "refseq": None,
                "altseq": None,
                "prefix": "r",
                "position_type": "rna",
                "suffix": None,
                "variant_type": "fusion",
            },
            "breakpoint2": {
                "feature": "NM_002609.3",
                "start": 1580,
                "start_offset": None,
                "start_seq": None,
                "end": 1924,
                "end_offset": None,
                "end_seq": None,
                "strand": None,
                "refseq": None,
                "altseq": None,
                "prefix": "r",
                "position_type": "rna",
                "suffix": None,
                "variant_type": "fusion",
            },
        },
    ),
]


@pytest.mark.parametrize("string, parsed", PARSE_CASES, ids=[i[0] for i in PARSE_CASES])
def test_parse(string, parsed):
    assert parse(string) == parsed


def test_protein_letters_3to1():
    assert protein_letters_3to1("Val") == "V"
    assert protein_letters_3to1("V") == "V"
    assert protein_letters_3to1("HisGln") == "HQ"
    assert protein_letters_3to1("HisQ") == "HQ"
    assert protein_letters_3to1("HisQVal") == "HQV"
    assert protein_letters_3to1("HQV") == "HQV"
    assert protein_letters_3to1("") == ""
