import pytest
from pyfaidx import Fasta

from coordinate_mapper import set_ensembl_release
from coordinate_mapper.ensembl import Ensembl, Genome
from coordinate_mapper.map import transcript_to_contig

from . import GENOME_FILE


def test_cache_clear():
    """If the cache is not cleared when changing the Ensembl release, running the same query
    will return the coordinates on the old release.
    """
    set_ensembl_release(100)
    pos = transcript_to_contig("ENST00000341165", 1)
    assert EnsemblRelease.release == 100
    assert EnsemblRelease.species.latin_name == "homo_sapiens"
    assert pos[0].start == 43170481

    set_ensembl_release(69)
    pos = transcript_to_contig("ENST00000341165", 1)
    assert EnsemblRelease.release == 69
    assert EnsemblRelease.species.latin_name == "homo_sapiens"
    assert pos[0].start == 41322498


def test_set_ensembl_release_by_release():
    set_ensembl_release(69)
    assert EnsemblRelease.release == 69
    assert EnsemblRelease.species.latin_name == "homo_sapiens"


def test_set_ensembl_release_by_release_2():
    set_ensembl_release(79)
    assert EnsemblRelease.release == 79
    assert EnsemblRelease.species.latin_name == "homo_sapiens"


def test_set_ensembl_release_by_name_hg19_1():
    set_ensembl_release("GRCh37")
    assert 54 < EnsemblRelease.release <= 75
    assert EnsemblRelease.species.latin_name == "homo_sapiens"


def test_set_ensembl_release_by_name_hg19_2():
    set_ensembl_release("hg19")
    assert 54 < EnsemblRelease.release <= 75
    assert EnsemblRelease.species.latin_name == "homo_sapiens"


def test_set_ensembl_release_by_name_hg19_3():
    set_ensembl_release("hg19a")
    assert 54 < EnsemblRelease.release <= 75
    assert EnsemblRelease.species.latin_name == "homo_sapiens"


def test_set_ensembl_release_by_name_hg38_1():
    set_ensembl_release("GRCh38")
    assert 75 < EnsemblRelease.release
    assert EnsemblRelease.species.latin_name == "homo_sapiens"


def test_set_ensembl_release_by_name_hg38_2():
    set_ensembl_release("hg38")
    assert 75 < EnsemblRelease.release
    assert EnsemblRelease.species.latin_name == "homo_sapiens"


def test_set_ensembl_release_by_name_hg38_3():
    set_ensembl_release("hg38_no_alt")
    assert 75 < EnsemblRelease.release
    assert EnsemblRelease.species.latin_name == "homo_sapiens"


def test_set_ensembl_release_by_name_hg18_1():
    set_ensembl_release("hg18")
    assert EnsemblRelease.release <= 54
    assert EnsemblRelease.species.latin_name == "homo_sapiens"


def test_genome_load():
    Genome.load(GENOME_FILE)
    assert isinstance(Genome._data, Fasta)
    assert "test" in Genome._data


def test_genome_sequence():
    Genome.load(GENOME_FILE)
    assert Genome.sequence("test", 1, 6) == "ATGTCC"
    with pytest.raises(ValueError):
        _ = Genome.sequence("test", 0)
        _ = Genome.sequence("test", 1, 999)
        _ = Genome.sequence("test", -6)
