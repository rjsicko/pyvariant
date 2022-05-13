from coordinate_mapper import set_ensembl_release
from coordinate_mapper.all_gene_names import GeneNames


def test_normalize_gene_id():
    assert GeneNames.normalize("ENSG00000000419") == "ENSG00000000419"


def test_normalize_gene_id_2():
    assert GeneNames.normalize("ENSG00000204645") == "ENSG00000204645"


def test_normalize_gene_name():
    assert GeneNames.normalize("DPM1") == "DPM1"


def test_normalize_gene_name_change():
    assert GeneNames.normalize("ABP1") == "AOC1"


def test_to_id_different_gene_names():
    # Gene name changed bewteen 69 and 100
    # ENSG00000118058 MLL
    # ENSG00000118058 KMT2A
    set_ensembl_release(69)
    assert GeneNames.to_id("MLL") == "ENSG00000118058"
    assert GeneNames.to_id("KMT2A") == "ENSG00000118058"


def test_to_id_different_gene_ids():
    # Gene ID changed bewteen 69 and 100
    # ENSG00000166748 AGBL1
    # ENSG00000273540 AGBL1
    set_ensembl_release(69)
    assert GeneNames.to_id("AGBL1") == "ENSG00000166748"
    set_ensembl_release(100)
    assert GeneNames.to_id("AGBL1") == "ENSG00000273540"


def test_to_id_different_gene_ids_2():
    # Gene ID changed bewteen 69 and 100
    # ENSG00000204645 SSX4
    # ENSG00000268009 SSX4
    set_ensembl_release(69)
    assert GeneNames.to_id("SSX4") == "ENSG00000204645"
    set_ensembl_release(100)
    assert GeneNames.to_id("SSX4") == "ENSG00000268009"


def test_to_name_different_gene_names():
    # Gene name changed bewteen 69 and 100
    # ENSG00000118058 MLL
    # ENSG00000118058 KMT2A
    set_ensembl_release(69)
    assert GeneNames.to_name("ENSG00000118058") == "MLL"
    set_ensembl_release(100)
    assert GeneNames.to_name("ENSG00000118058") == "KMT2A"


def test_to_name_different_gene_ids():
    # Gene ID changed bewteen 69 and 100
    # ENSG00000204645 SSX4
    # ENSG00000268009 SSX4
    set_ensembl_release(69)
    assert GeneNames.to_name("ENSG00000204645") == "SSX4"
    set_ensembl_release(100)
    assert GeneNames.to_name("ENSG00000268009") == "SSX4"
