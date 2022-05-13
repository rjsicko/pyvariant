from ensembl_map.ensembl_download import download, make_gtf_filename, make_gtf_url


def test_make_gtf_filename():
    assert make_gtf_filename("homo_sapiens", "GRCh38", 100) == "Homo_sapiens.GRCh38.100.gtf.gz"


def test_make_gtf_url():
    assert (
        make_gtf_url("homo_sapiens", "GRCh38", 100)
        == "ftp://ftp.ensembl.org/pub/release-100/gtf/homo_sapiens/Homo_sapiens.GRCh38.100.gtf.gz"
    )
