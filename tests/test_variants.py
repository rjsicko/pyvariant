import pytest
from constants import CACHE_DIR

from ensembl_map.core import EnsemblRelease
from ensembl_map.variants import (
    CdsDeletion,
    CdsDelins,
    CdsDuplication,
    CdsFusion,
    CdsInsertion,
    CdsSubstitution,
    DnaDeletion,
    DnaDelins,
    DnaDuplication,
    DnaFusion,
    DnaInsertion,
    DnaSubstitution,
    ExonFusion,
    ProteinDeletion,
    ProteinDelins,
    ProteinDuplication,
    ProteinFrameShift,
    ProteinInsertion,
    ProteinSubstitution,
    RnaDeletion,
    RnaDelins,
    RnaDuplication,
    RnaFusion,
    RnaInsertion,
    RnaSubstitution,
)


@pytest.fixture(scope="module")
def ensembl69():
    return EnsemblRelease(
        species="homo_sapiens",
        release=69,
        cache_dir=CACHE_DIR,
        canonical_transcript="",
        contig_alias="",
        exon_alias="",
        gene_alias="",
        protein_alias="",
        transcript_alias="",
    )


# class TestCdsDeletion:
@pytest.fixture
def variant():
    return CdsDeletion(reference="VHL", reference_type="gene", start=478, end=480, ref_seq="GAG")


def test_hgvs(variant):
    assert variant.hgvs() == "VHL:c.478_480del"


def test_to_cds(variant):
    result = variant.to_cds()
    assert len(result) == 1
    assert result[0].reference == "ENST00000256474"
    assert result[0].start == 478
    assert result[0].end == 480
    assert result[0].strand == "+"
    assert result[0].ref_seq == "GAG"


def test_to_contig(variant):
    result = variant.to_contig()
    assert len(result) == 1
    assert result[0].reference == "3"
    assert result[0].start == 10191485
    assert result[0].end == 10191487
    assert result[0].strand == "+"
    assert result[0].ref_seq == "GAG"


def test_to_gene(variant):
    result = variant.to_gene()
    assert len(result) == 1
    assert result[0].reference == "VHL"
    assert result[0].start == 10191485
    assert result[0].end == 10191487
    assert result[0].strand == "+"
    assert result[0].ref_seq == "GAG"


def test_to_protein(variant):
    result = variant.to_protein()
    assert len(result) == 1
    assert result[0].reference == "ENSP00000256474"
    assert result[0].start == 160
    assert result[0].end == 160
    assert result[0].strand == "+"
    assert result[0].ref_seq == "E"


def test_to_transcript(variant):
    result = variant.to_rna()
    assert len(result) == 1
    assert result[0].reference == "ENST00000256474"
    assert result[0].start == 1318
    assert result[0].end == 1320
    assert result[0].strand == "+"
    assert result[0].ref_seq == "GAG"


def test_values(variant):
    assert variant.values() == {
        "reference": "VHL",
        "reference_type": "gene",
        "start": 478,
        "end": 480,
        "alt_seq": "",
        "ref_seq": "GAG",
        "strand": "+",
        "start_offset": 0,
        "start_seq": "G",
        "end_offset": 0,
        "end_seq": "G",
        "user_repr": None,
        "graphkb_id": None,
    }


# class TestCdsDeletionError:
def test_to_protein_inter_codon_error(variant):
    return CdsDeletion(
        reference="ENST00000256474", reference_type="transcript", start=479, end=481, ref_seq="AGA"
    )
    with pytest.raises(NotImplementedError):
        variant.to_protein()


def test_to_protein_ref_mismatch_error(variant):
    return CdsDeletion(
        reference="ENST00000256474", reference_type="transcript", start=478, end=480, ref_seq="NNN"
    )
    with pytest.raises(RefMismatchError):
        variant.to_protein()


# class TestCdsDeletionGeneName:
@pytest.fixture
def variant():
    return CdsDeletion(reference="VHL", reference_type="gene", start=478, end=480, ref_seq="GAG")


def test_to_cds(variant):
    result = variant.to_cds(gene_name=True)
    assert len(result) == 1
    assert result[0].reference == "VHL"
    assert result[0].start == 478
    assert result[0].end == 480
    assert result[0].strand == "+"
    assert result[0].ref_seq == "GAG"


def test_to_contig(variant):
    result = variant.to_contig(gene_name=True)
    assert len(result) == 1
    assert result[0].reference == "VHL"
    assert result[0].start == 10191485
    assert result[0].end == 10191487
    assert result[0].strand == "+"
    assert result[0].ref_seq == "GAG"


def test_to_gene(variant):
    result = variant.to_gene(gene_name=True)
    assert len(result) == 1
    assert result[0].reference == "VHL"
    assert result[0].start == 10191485
    assert result[0].end == 10191487
    assert result[0].strand == "+"
    assert result[0].ref_seq == "GAG"


def test_to_protein(variant):
    result = variant.to_protein(gene_name=True)
    assert len(result) == 1
    assert result[0].reference == "VHL"
    assert result[0].start == 160
    assert result[0].end == 160
    assert result[0].strand == "+"
    assert result[0].ref_seq == "E"


def test_to_transcript(variant):
    result = variant.to_rna(gene_name=True)
    assert len(result) == 1
    assert result[0].reference == "VHL"
    assert result[0].start == 1318
    assert result[0].end == 1320
    assert result[0].strand == "+"
    assert result[0].ref_seq == "GAG"


# class TestCdsDeletionOffsetSingle:
# NM_000264.3:c.-8456delG
@pytest.fixture
def variant():
    return CdsDeletion(
        reference="ENST00000331920",
        reference_type="transcript",
        start=1,
        start_offset=-8456,
        ref_seq="G",
    )


def test_hgvs(variant):
    assert variant.hgvs() == "ENST00000331920:c.-8456del"


def test_to_contig(variant):
    result = variant.to_contig()
    assert len(result) == 1
    assert result[0].reference == "9"
    assert result[0].start == 98279099
    assert result[0].end == 98279099
    assert result[0].strand == "+"
    assert result[0].ref_seq == "C"


# class TestCdsDeletionOffsetMultiple:
# VHL:c.-1_20del
@pytest.fixture
def variant():
    return CdsDeletion(
        reference="VHL",
        reference_type="gene",
        start=1,
        start_offset=-1,
        end=20,
        end_offset=0,
        ref_seq="AATGCCCCGGAGGGCGGAGAA",
    )


def test_hgvs(variant):
    assert variant.hgvs() == "VHL:c.-1_20del"


def test_to_contig(variant):
    result = variant.to_contig()
    assert len(result) == 1
    assert result[0].reference == "3"
    assert result[0].start == 10183531
    assert result[0].end == 10183551
    assert result[0].strand == "+"
    assert result[0].ref_seq == "AATGCCCCGGAGGGCGGAGAA"


# class TestCdsDeletionToProteinFrameShift:
@pytest.fixture
def variant():
    return CdsDeletion(
        reference="ENST00000257430",
        reference_type="transcript",
        start=3927,
        end=3931,
        ref_seq="AAAGA",
    )


def test_to_protein(variant):
    result = variant.to_protein()
    assert len(result) == 1
    assert result[0].reference == "ENSP00000257430"
    assert result[0].start == 1309
    assert result[0].end == 1309
    assert result[0].strand == "+"
    assert result[0].ref_seq == "E"
    assert result[0].hgvs() == "ENSP00000257430:p.E1309fs"


# class TestCdsDeletionRefSeq:
# NM_000264.3:c.-8456delG
@pytest.fixture
def variant():
    return CdsDeletion(reference="NM_000264.3", reference_type="transcript", start=1, ref_seq="A")


def test_hgvs(variant):
    assert variant.hgvs() == "NM_000264.3:c.1del"


def test_to_cds(variant):
    result = variant.to_cds()
    assert len(result) == 1
    assert result[0].reference == "ENST00000331920"
    assert result[0].start == 1
    assert result[0].strand == "-"
    assert result[0].ref_seq == "A"
    assert result[0].hgvs() == "ENST00000331920:c.1del"


# class TestDnaDeletion:
@pytest.fixture
def variant():
    return DnaDeletion(
        reference="17", reference_type="contig", start=7579471, end=7579473, ref_seq="GGG"
    )


def test_hgvs(variant):
    assert variant.hgvs() == "17:g.7579471_7579473del"


def test_to_cds(variant):
    result = variant.to_cds()
    assert len(result) == 6
    assert result[0].reference == "ENST00000269305"
    assert result[0].start == 214
    assert result[0].end == 216
    assert result[0].strand == "-"
    assert result[0].ref_seq == "CCC"


def test_to_contig(variant):
    result = variant.to_contig()
    assert len(result) == 1
    assert result[0].reference == "17"
    assert result[0].start == 7579471
    assert result[0].end == 7579473
    assert result[0].strand == "+"
    assert result[0].ref_seq == "GGG"


def test_to_gene(variant):
    result = variant.to_gene()
    assert len(result) == 1
    assert result[0].reference == "TP53"
    assert result[0].start == 7579471
    assert result[0].end == 7579473
    assert result[0].strand == "-"
    assert result[0].ref_seq == "CCC"


def test_to_protein(variant):
    result = variant.to_protein()
    assert len(result) == 6
    assert result[0].reference == "ENSP00000269305"
    assert result[0].start == 72
    assert result[0].end == 72
    assert result[0].strand == "-"
    assert result[0].ref_seq == "P"


def test_to_transcript(variant):
    result = variant.to_rna()
    assert len(result) == 9
    assert result[0].reference == "ENST00000269305"
    assert result[0].start == 404
    assert result[0].end == 406
    assert result[0].strand == "-"
    assert result[0].ref_seq == "CCC"


# class TestDnaDeletionGeneName:
@pytest.fixture
def variant():
    return DnaDeletion(
        reference="17", reference_type="contig", start=7579471, end=7579473, ref_seq="GGG"
    )


def test_to_cds(variant):
    result = variant.to_cds(gene_name=True)
    assert len(result) == 6
    assert result[0].reference == "TP53"
    assert result[0].start == 214
    assert result[0].end == 216
    assert result[0].strand == "-"
    assert result[0].ref_seq == "CCC"


def test_to_contig(variant):
    result = variant.to_contig(gene_name=True)
    assert len(result) == 1
    assert result[0].reference == "TP53"
    assert result[0].start == 7579471
    assert result[0].end == 7579473
    assert result[0].strand == "+"
    assert result[0].ref_seq == "GGG"


def test_to_gene(variant):
    result = variant.to_gene(gene_name=True)
    assert len(result) == 1
    assert result[0].reference == "TP53"
    assert result[0].start == 7579471
    assert result[0].end == 7579473
    assert result[0].strand == "-"
    assert result[0].ref_seq == "CCC"


def test_to_protein(variant):
    result = variant.to_protein(gene_name=True)
    assert len(result) == 1
    assert result[0].reference == "TP53"
    assert result[0].start == 72
    assert result[0].end == 72
    assert result[0].strand == "-"
    assert result[0].ref_seq == "P"


def test_to_transcript(variant):
    result = variant.to_rna(gene_name=True)
    assert len(result) == 9
    assert result[0].reference == "TP53"
    assert result[0].start == 404
    assert result[0].end == 406
    assert result[0].strand == "-"
    assert result[0].ref_seq == "CCC"


# class TestProteinDeletion:
@pytest.fixture
def variant():
    return ProteinDeletion(reference="RET", reference_type="gene", start=632, end=633, ref_seq="EL")


def test_hgvs(variant):
    assert variant.hgvs() == "RET:p.E632_L633del"


def test_to_cds(variant):
    result = variant.to_cds()
    assert len(result) == 2
    assert result[0].reference == "ENST00000340058"
    assert result[0].start == 1894
    assert result[0].end == 1899
    assert result[0].strand == "+"
    assert result[0].ref_seq == "GAGCTG"


def test_to_contig(variant):
    result = variant.to_contig()
    assert len(result) == 1
    assert result[0].reference == "10"
    assert result[0].start == 43609942
    assert result[0].end == 43609947
    assert result[0].strand == "+"
    assert result[0].ref_seq == "GAGCTG"


def test_to_gene(variant):
    result = variant.to_gene()
    assert len(result) == 1
    assert result[0].reference == "RET"
    assert result[0].start == 43609942
    assert result[0].end == 43609947
    assert result[0].strand == "+"
    assert result[0].ref_seq == "GAGCTG"


def test_to_protein(variant):
    result = variant.to_protein()
    assert len(result) == 2
    assert result[0].reference == "ENSP00000344798"
    assert result[0].start == 632
    assert result[0].end == 633
    assert result[0].strand == "+"
    assert result[0].ref_seq == "EL"


def test_to_transcript(variant):
    result = variant.to_rna()
    assert len(result) == 2
    assert result[0].reference == "ENST00000340058"
    assert result[0].start == 2074
    assert result[0].end == 2079
    assert result[0].strand == "+"
    assert result[0].ref_seq == "GAGCTG"


# class TestProteinDeletionGeneName:
@pytest.fixture
def variant():
    return ProteinDeletion(reference="RET", reference_type="gene", start=632, end=633, ref_seq="EL")


def test_to_cds(variant):
    result = variant.to_cds(gene_name=True)
    assert len(result) == 2
    assert result[0].reference == "RET"
    assert result[0].start == 1894
    assert result[0].end == 1899
    assert result[0].strand == "+"
    assert result[0].ref_seq == "GAGCTG"


def test_to_contig(variant):
    result = variant.to_contig(gene_name=True)
    assert len(result) == 1
    assert result[0].reference == "RET"
    assert result[0].start == 43609942
    assert result[0].end == 43609947
    assert result[0].strand == "+"
    assert result[0].ref_seq == "GAGCTG"


def test_to_gene(variant):
    result = variant.to_gene(gene_name=True)
    assert len(result) == 1
    assert result[0].reference == "RET"
    assert result[0].start == 43609942
    assert result[0].end == 43609947
    assert result[0].strand == "+"
    assert result[0].ref_seq == "GAGCTG"


def test_to_protein(variant):
    result = variant.to_protein(gene_name=True)
    assert len(result) == 2
    assert result[0].reference == "RET"
    assert result[0].start == 632
    assert result[0].end == 633
    assert result[0].strand == "+"
    assert result[0].ref_seq == "EL"


def test_to_transcript(variant):
    result = variant.to_rna(gene_name=True)
    assert len(result) == 2
    assert result[0].reference == "RET"
    assert result[0].start == 2074
    assert result[0].end == 2079
    assert result[0].strand == "+"
    assert result[0].ref_seq == "GAGCTG"


# class TestRnaDeletion:
@pytest.fixture
def variant():
    return RnaDeletion(reference="VHL", reference_type="gene", start=1318, end=1320, ref_seq="GAG")


def test_hgvs(variant):
    assert variant.hgvs() == "VHL:r.1318_1320del"


def test_to_cds(variant):
    result = variant.to_cds()
    assert len(result) == 1
    assert result[0].reference == "ENST00000256474"
    assert result[0].start == 478
    assert result[0].end == 480
    assert result[0].strand == "+"
    assert result[0].ref_seq == "GAG"


def test_to_contig(variant):
    result = variant.to_contig()
    assert len(result) == 1
    assert result[0].reference == "3"
    assert result[0].start == 10191485
    assert result[0].end == 10191487
    assert result[0].strand == "+"
    assert result[0].ref_seq == "GAG"


def test_to_gene(variant):
    result = variant.to_gene()
    assert len(result) == 1
    assert result[0].reference == "VHL"
    assert result[0].start == 10191485
    assert result[0].end == 10191487
    assert result[0].strand == "+"
    assert result[0].ref_seq == "GAG"


def test_to_protein(variant):
    result = variant.to_protein()
    assert len(result) == 1
    assert result[0].reference == "ENSP00000256474"
    assert result[0].start == 160
    assert result[0].end == 160
    assert result[0].strand == "+"
    assert result[0].ref_seq == "E"


def test_to_transcript(variant):
    result = variant.to_rna()
    assert len(result) == 1
    assert result[0].reference == "ENST00000256474"
    assert result[0].start == 1318
    assert result[0].end == 1320
    assert result[0].strand == "+"
    assert result[0].ref_seq == "GAG"


# class TestRnaDeletionGeneName:
@pytest.fixture
def variant():
    return RnaDeletion(reference="VHL", reference_type="gene", start=1318, end=1320, ref_seq="GAG")


def test_to_cds(variant):
    result = variant.to_cds(gene_name=True)
    assert len(result) == 1
    assert result[0].reference == "VHL"
    assert result[0].start == 478
    assert result[0].end == 480
    assert result[0].strand == "+"
    assert result[0].ref_seq == "GAG"


def test_to_contig(variant):
    result = variant.to_contig(gene_name=True)
    assert len(result) == 1
    assert result[0].reference == "VHL"
    assert result[0].start == 10191485
    assert result[0].end == 10191487
    assert result[0].strand == "+"
    assert result[0].ref_seq == "GAG"


def test_to_gene(variant):
    result = variant.to_gene(gene_name=True)
    assert len(result) == 1
    assert result[0].reference == "VHL"
    assert result[0].start == 10191485
    assert result[0].end == 10191487
    assert result[0].strand == "+"
    assert result[0].ref_seq == "GAG"


def test_to_protein(variant):
    result = variant.to_protein(gene_name=True)
    assert len(result) == 1
    assert result[0].reference == "VHL"
    assert result[0].start == 160
    assert result[0].end == 160
    assert result[0].strand == "+"
    assert result[0].ref_seq == "E"


def test_to_transcript(variant):
    result = variant.to_rna(gene_name=True)
    assert len(result) == 1
    assert result[0].reference == "VHL"
    assert result[0].start == 1318
    assert result[0].end == 1320
    assert result[0].strand == "+"
    assert result[0].ref_seq == "GAG"


# class TestProteinFrameShift:
@pytest.fixture
def variant():
    return ProteinFrameShift(
        reference="ENSP00000257430", reference_type="protein", start=1309, ref_seq="E"
    )


def test_hgvs(variant):
    assert variant.hgvs() == "ENSP00000257430:p.E1309fs"


def test_to_protein(variant):
    result = variant.to_protein()
    assert len(result) == 1
    assert result[0].reference == "ENSP00000257430"
    assert result[0].start == 1309
    assert result[0].end == 1309
    assert result[0].strand == "+"
    assert result[0].ref_seq == "E"
    assert result[0].alt_seq == ""


# class TestCdsDelins:
@pytest.fixture
def variant():
    return CdsDelins(
        reference="ENST00000269305",
        reference_type="transcript",
        start=878,
        end=880,
        ref_seq="GGG",
        alt_seq="TTT",
    )


def test_hgvs(variant):
    assert variant.hgvs() == "ENST00000269305:c.878_880delinsTTT"


def test_to_cds(variant):
    result = variant.to_cds()
    assert len(result) == 1
    assert result[0].reference == "ENST00000269305"
    assert result[0].start == 878
    assert result[0].end == 880
    assert result[0].strand == "-"
    assert result[0].ref_seq == "GGG"
    assert result[0].alt_seq == "TTT"


def test_to_contig(variant):
    result = variant.to_contig()
    assert len(result) == 1
    assert result[0].reference == "17"
    assert result[0].start == 7577058
    assert result[0].end == 7577060
    assert result[0].strand == "+"
    assert result[0].ref_seq == "CCC"
    assert result[0].alt_seq == "AAA"


def test_to_gene(variant):
    result = variant.to_gene()
    assert len(result) == 1
    assert result[0].reference == "TP53"
    assert result[0].start == 7577058
    assert result[0].end == 7577060
    assert result[0].strand == "-"
    assert result[0].ref_seq == "GGG"
    assert result[0].alt_seq == "TTT"


def test_to_protein(variant):
    result = variant.to_protein()
    assert len(result) == 1
    assert result[0].reference == "ENSP00000269305"
    assert result[0].start == 293
    assert result[0].end == 294
    assert result[0].strand == "-"
    assert result[0].ref_seq == "GE"
    assert result[0].alt_seq == "V*"


def test_to_transcript(variant):
    result = variant.to_rna()
    assert len(result) == 1
    assert result[0].reference == "ENST00000269305"
    assert result[0].start == 1068
    assert result[0].end == 1070
    assert result[0].strand == "-"
    assert result[0].ref_seq == "GGG"
    assert result[0].alt_seq == "TTT"


def test_values(variant):
    assert variant.values() == {
        "reference": "ENST00000269305",
        "reference_type": "transcript",
        "start": 878,
        "end": 880,
        "alt_seq": "TTT",
        "ref_seq": "GGG",
        "strand": "-",
        "start_offset": 0,
        "start_seq": "G",
        "end_offset": 0,
        "end_seq": "G",
        "user_repr": None,
        "graphkb_id": None,
    }


# class TestCdsDelinsError:
def test_to_protein_ref_mismatch_error(variant):
    return CdsDelins(
        reference="ENST00000269305",
        reference_type="transcript",
        start=878,
        end=880,
        ref_seq="NNN",
        alt_seq="TTT",
    )
    with pytest.raises(RefMismatchError):
        variant.to_protein()


# class TestCdsDelinsGeneName:
@pytest.fixture
def variant():
    return CdsDelins(
        reference="ENST00000269305",
        reference_type="transcript",
        start=878,
        end=880,
        ref_seq="GGG",
        alt_seq="TTT",
    )


def test_to_cds(variant):
    result = variant.to_cds(gene_name=True)
    assert len(result) == 1
    assert result[0].reference == "TP53"
    assert result[0].start == 878
    assert result[0].end == 880
    assert result[0].strand == "-"
    assert result[0].ref_seq == "GGG"
    assert result[0].alt_seq == "TTT"


def test_to_contig(variant):
    result = variant.to_contig(gene_name=True)
    assert len(result) == 1
    assert result[0].reference == "TP53"
    assert result[0].start == 7577058
    assert result[0].end == 7577060
    assert result[0].strand == "+"
    assert result[0].ref_seq == "CCC"
    assert result[0].alt_seq == "AAA"


def test_to_gene(variant):
    result = variant.to_gene(gene_name=True)
    assert len(result) == 1
    assert result[0].reference == "TP53"
    assert result[0].start == 7577058
    assert result[0].end == 7577060
    assert result[0].strand == "-"
    assert result[0].ref_seq == "GGG"
    assert result[0].alt_seq == "TTT"


def test_to_protein(variant):
    result = variant.to_protein(gene_name=True)
    assert len(result) == 1
    assert result[0].reference == "TP53"
    assert result[0].start == 293
    assert result[0].end == 294
    assert result[0].strand == "-"
    assert result[0].ref_seq == "GE"
    assert result[0].alt_seq == "V*"


def test_to_transcript(variant):
    result = variant.to_rna(gene_name=True)
    assert len(result) == 1
    assert result[0].reference == "TP53"
    assert result[0].start == 1068
    assert result[0].end == 1070
    assert result[0].strand == "-"
    assert result[0].ref_seq == "GGG"
    assert result[0].alt_seq == "TTT"


# class TestCdsDelinsOffsetSingle:
@pytest.fixture
def variant():
    return CdsDelins(
        reference="ENST00000323250",
        reference_type="transcript",
        start=1,
        start_offset=-86,
        ref_seq="G",
        alt_seq="CA",
    )


def test_hgvs(variant):
    assert variant.hgvs() == "ENST00000323250:c.-86delinsCA"


def test_to_contig(variant):
    result = variant.to_contig()
    assert result[0].reference == "18"
    assert result[0].start == 657657
    assert result[0].end == 657657
    assert result[0].strand == "+"
    assert result[0].ref_seq == "G"
    assert result[0].alt_seq == "CA"


# class TestCdsDelinsOffsetMultiple:
# NM_001354868.2:c.-86_-13delinsCGCCTGCCTCCGTCCCGCCGCGCCACTTCGCCTGCCTCCGTCCCGCCGCGCCACTTCGCCTGCCTCCGTCCCGCCGCGCCACTTCGCCTGCCTCCGTCCCG
@pytest.fixture
def variant():
    return CdsDelins(
        reference="ENST00000323250",
        reference_type="transcript",
        start=1,
        end=1,
        start_offset=-86,
        end_offset=-13,
        ref_seq="GGCCTGCCTCCGTCCCGCCGCGCCACTTGGCCTGCCTCCGTCCCGCCGCGCCACTTCGCCTGCCTCCGTCCCCC",
        alt_seq="CGCCTGCCTCCGTCCCGCCGCGCCACTTCGCCTGCCTCCGTCCCGCCGCGCCACTTCGCCTGCCTCCGTCCCGCCGCGCCACTTCGCCTGCCTCCGTCCCG",
    )


def test_hgvs(variant):
    assert variant.hgvs() == f"ENST00000323250:c.-86_-13delins{variant.alt_seq}"


def test_to_contig(variant):
    result = variant.to_contig()
    assert result[0].reference == "18"
    assert result[0].start == 657657
    assert result[0].end == 657730
    assert result[0].strand == "+"
    assert result[0].ref_seq == variant.ref_seq
    assert result[0].alt_seq == variant.alt_seq


# class TestCdsDelinsOffsetLengthOffsetAcrossZero:
@pytest.fixture
def variant():
    return CdsDelins(
        reference="ENST00000323250",
        reference_type="transcript",
        start=1,
        start_offset=-2,
        end=2,
        end_offset=0,
        alt_seq="N",
    )


def test_to_contig(variant):
    result = variant.to_contig()
    assert result[0].length == 4


# class TestCdsDelinsSynonymous:
@pytest.fixture
def variant():
    return CdsDelins(
        reference="ENST00000372098",
        reference_type="transcript",
        start=1039,
        end=1041,
        ref_seq="GAC",
        alt_seq="GAT",
    )


def test_to_protein(variant):
    result = variant.to_protein()
    assert len(result) == 1
    assert result[0].reference == "ENSP00000361170"
    assert result[0].start == 347
    assert result[0].strand == "-"
    assert result[0].ref_seq == "D"  # GAC == D
    assert result[0].alt_seq == "D"  # GAT == D
    assert result[0].hgvs() == "ENSP00000361170:p.D347="


# class TestDnaDelins:
@pytest.fixture
def variant():
    return DnaDelins(
        reference="TP53",
        reference_type="gene",
        start=7577058,
        end=7577060,
        ref_seq="GGG",
        alt_seq="TTT",
    )


def test_hgvs(variant):
    assert variant.hgvs() == "TP53:g.7577058_7577060delinsTTT"


def test_to_cds(variant):
    result = variant.to_cds()
    assert len(result) == 5
    assert result[0].reference == "ENST00000269305"
    assert result[0].start == 878
    assert result[0].end == 880
    assert result[0].strand == "-"
    assert result[0].ref_seq == "GGG"
    assert result[0].alt_seq == "TTT"


def test_to_contig(variant):
    result = variant.to_contig()
    assert len(result) == 1
    assert result[0].reference == "17"
    assert result[0].start == 7577058
    assert result[0].end == 7577060
    assert result[0].strand == "+"
    assert result[0].ref_seq == "CCC"
    assert result[0].alt_seq == "AAA"


def test_to_gene(variant):
    result = variant.to_gene()
    assert len(result) == 1
    assert result[0].reference == "TP53"
    assert result[0].start == 7577058
    assert result[0].end == 7577060
    assert result[0].strand == "-"
    assert result[0].ref_seq == "GGG"
    assert result[0].alt_seq == "TTT"


def test_to_protein(variant):
    result = variant.to_protein()
    assert len(result) == 5
    assert result[0].reference == "ENSP00000269305"
    assert result[0].start == 293
    assert result[0].end == 294
    assert result[0].strand == "-"
    assert result[0].ref_seq == "GE"
    assert result[0].alt_seq == "V*"


def test_to_transcript(variant):
    result = variant.to_rna()
    assert len(result) == 9
    assert result[0].reference == "ENST00000269305"
    assert result[0].start == 1068
    assert result[0].end == 1070
    assert result[0].strand == "-"
    assert result[0].ref_seq == "GGG"
    assert result[0].alt_seq == "TTT"


# class TestDnaDelinsGeneName:
@pytest.fixture
def variant():
    return DnaDelins(
        reference="17",
        reference_type="contig",
        start=7577058,
        end=7577060,
        ref_seq="CCC",
        alt_seq="AAA",
    )


def test_to_cds(variant):
    result = variant.to_cds(gene_name=True)
    assert len(result) == 5
    assert result[0].reference == "TP53"
    assert result[0].start == 878
    assert result[0].end == 880
    assert result[0].strand == "-"
    assert result[0].ref_seq == "GGG"
    assert result[0].alt_seq == "TTT"


def test_to_contig(variant):
    result = variant.to_contig(gene_name=True)
    assert len(result) == 1
    assert result[0].reference == "TP53"
    assert result[0].start == 7577058
    assert result[0].end == 7577060
    assert result[0].strand == "+"
    assert result[0].ref_seq == "CCC"
    assert result[0].alt_seq == "AAA"


def test_to_gene(variant):
    result = variant.to_gene(gene_name=True)
    assert len(result) == 1
    assert result[0].reference == "TP53"
    assert result[0].start == 7577058
    assert result[0].end == 7577060
    assert result[0].strand == "-"
    assert result[0].ref_seq == "GGG"
    assert result[0].alt_seq == "TTT"


def test_to_protein(variant):
    result = variant.to_protein(gene_name=True)
    assert len(result) == 1
    assert result[0].reference == "TP53"
    assert result[0].start == 293
    assert result[0].end == 294
    assert result[0].strand == "-"
    assert result[0].ref_seq == "GE"
    assert result[0].alt_seq == "V*"


def test_to_transcript(variant):
    result = variant.to_rna(gene_name=True)
    assert len(result) == 9
    assert result[0].reference == "TP53"
    assert result[0].start == 1068
    assert result[0].end == 1070
    assert result[0].strand == "-"
    assert result[0].ref_seq == "GGG"
    assert result[0].alt_seq == "TTT"


# class TestProteinDelins:
@pytest.fixture
def variant():
    return ProteinDelins(
        reference="KIT", reference_type="gene", start=417, end=419, ref_seq="TYD", alt_seq="I"
    )


def test_hgvs(variant):
    assert variant.hgvs() == "KIT:p.T417_D419delinsI"


def test_to_cds(variant):
    result = variant.to_cds()
    assert len(result) == 6
    assert result[0].reference == "ENST00000288135"
    assert result[0].start == 1249
    assert result[0].end == 1257
    assert result[0].strand == "+"
    assert result[0].ref_seq == "ACTTACGAC"
    assert result[0].alt_seq == "ATA"


def test_to_contig(variant):
    result = variant.to_contig()
    assert len(result) == 3
    assert result[0].reference == "4"
    assert result[0].start == 55589767
    assert result[0].end == 55589775
    assert result[0].strand == "+"
    assert result[0].ref_seq == "ACTTACGAC"
    assert result[0].alt_seq == "ATA"


def test_to_gene(variant):
    result = variant.to_gene()
    assert len(result) == 3
    assert result[0].reference == "KIT"
    assert result[0].start == 55589767
    assert result[0].end == 55589775
    assert result[0].strand == "+"
    assert result[0].ref_seq == "ACTTACGAC"
    assert result[0].alt_seq == "ATA"


def test_to_protein(variant):
    result = variant.to_protein()
    assert len(result) == 2
    assert result[0].reference == "ENSP00000288135"
    assert result[0].start == 417
    assert result[0].end == 419
    assert result[0].strand == "+"
    assert result[0].ref_seq == "TYD"
    assert result[0].alt_seq == "I"


def test_to_transcript(variant):
    result = variant.to_rna()
    assert len(result) == 6
    assert result[0].reference == "ENST00000288135"
    assert result[0].start == 1346
    assert result[0].end == 1354
    assert result[0].strand == "+"
    assert result[0].ref_seq == "ACTTACGAC"
    assert result[0].alt_seq == "ATA"


# class TestProteinDelinsGeneName:
@pytest.fixture
def variant():
    return ProteinDelins(
        reference="KIT", reference_type="gene", start=417, end=419, ref_seq="TYD", alt_seq="I"
    )


def test_to_cds(variant):
    result = variant.to_cds(gene_name=True)
    assert len(result) == 6
    assert result[0].reference == "KIT"
    assert result[0].start == 1249
    assert result[0].end == 1257
    assert result[0].strand == "+"
    assert result[0].ref_seq == "ACTTACGAC"
    assert result[0].alt_seq == "ATA"


def test_to_contig(variant):
    result = variant.to_contig(gene_name=True)
    assert len(result) == 3
    assert result[0].reference == "KIT"
    assert result[0].start == 55589767
    assert result[0].end == 55589775
    assert result[0].strand == "+"
    assert result[0].ref_seq == "ACTTACGAC"
    assert result[0].alt_seq == "ATA"


def test_to_gene(variant):
    result = variant.to_gene(gene_name=True)
    assert len(result) == 3
    assert result[0].reference == "KIT"
    assert result[0].start == 55589767
    assert result[0].end == 55589775
    assert result[0].strand == "+"
    assert result[0].ref_seq == "ACTTACGAC"
    assert result[0].alt_seq == "ATA"


def test_to_protein(variant):
    result = variant.to_protein(gene_name=True)
    assert len(result) == 2
    assert result[0].reference == "KIT"
    assert result[0].start == 417
    assert result[0].end == 419
    assert result[0].strand == "+"
    assert result[0].ref_seq == "TYD"
    assert result[0].alt_seq == "I"


def test_to_transcript(variant):
    result = variant.to_rna(gene_name=True)
    assert len(result) == 6
    assert result[0].reference == "KIT"
    assert result[0].start == 1346
    assert result[0].end == 1354
    assert result[0].strand == "+"
    assert result[0].ref_seq == "ACTTACGAC"
    assert result[0].alt_seq == "ATA"


# class TestRnaDelins:
"""Test mapping RNA delins to other types of delins."""


@pytest.fixture
def variant():
    return RnaDelins(
        reference="ENST00000288135",
        reference_type="transcript",
        start=976,
        end=977,
        ref_seq="TA",
        alt_seq="AG",
    )


def test_hgvs(variant):
    assert variant.hgvs() == "ENST00000288135:r.976_977delinsAG"


def test_to_cds(variant):
    result = variant.to_cds()
    assert len(result) == 1
    assert result[0].reference == "ENST00000288135"
    assert result[0].start == 879
    assert result[0].end == 880
    assert result[0].strand == "+"
    assert result[0].ref_seq == "TA"
    assert result[0].alt_seq == "AG"


def test_to_contig(variant):
    result = variant.to_contig()
    assert len(result) == 1
    assert result[0].reference == "4"
    assert result[0].start == 55570012
    assert result[0].end == 55570013
    assert result[0].strand == "+"
    assert result[0].ref_seq == "TA"
    assert result[0].alt_seq == "AG"


def test_to_gene(variant):
    result = variant.to_gene()
    assert len(result) == 1
    assert result[0].reference == "KIT"
    assert result[0].start == 55570012
    assert result[0].end == 55570013
    assert result[0].strand == "+"
    assert result[0].ref_seq == "TA"
    assert result[0].alt_seq == "AG"


def test_to_protein(variant):
    result = variant.to_protein()
    assert len(result) == 1
    assert result[0].reference == "ENSP00000288135"
    assert result[0].start == 293
    assert result[0].end == 294
    assert result[0].strand == "+"
    assert result[0].ref_seq == "NN"
    assert result[0].alt_seq == "KD"


def test_to_transcript(variant):
    result = variant.to_rna()
    assert len(result) == 1
    assert result[0].reference == "ENST00000288135"
    assert result[0].start == 976
    assert result[0].end == 977
    assert result[0].strand == "+"
    assert result[0].ref_seq == "TA"
    assert result[0].alt_seq == "AG"


# class TestRnaDelinsGeneName:
"""Test mapping RNA delins to other types of delins."""


@pytest.fixture
def variant():
    return RnaDelins(
        reference="ENST00000288135",
        reference_type="transcript",
        start=976,
        end=977,
        ref_seq="TA",
        alt_seq="AG",
    )


def test_to_cds(variant):
    result = variant.to_cds(gene_name=True)
    assert len(result) == 1
    assert result[0].reference == "KIT"
    assert result[0].start == 879
    assert result[0].end == 880
    assert result[0].strand == "+"
    assert result[0].ref_seq == "TA"
    assert result[0].alt_seq == "AG"


def test_to_contig(variant):
    result = variant.to_contig(gene_name=True)
    assert len(result) == 1
    assert result[0].reference == "KIT"
    assert result[0].start == 55570012
    assert result[0].end == 55570013
    assert result[0].strand == "+"
    assert result[0].ref_seq == "TA"
    assert result[0].alt_seq == "AG"


def test_to_gene(variant):
    result = variant.to_gene(gene_name=True)
    assert len(result) == 1
    assert result[0].reference == "KIT"
    assert result[0].start == 55570012
    assert result[0].end == 55570013
    assert result[0].strand == "+"
    assert result[0].ref_seq == "TA"
    assert result[0].alt_seq == "AG"


def test_to_protein(variant):
    result = variant.to_protein(gene_name=True)
    assert len(result) == 1
    assert result[0].reference == "KIT"
    assert result[0].start == 293
    assert result[0].end == 294
    assert result[0].strand == "+"
    assert result[0].ref_seq == "NN"
    assert result[0].alt_seq == "KD"


def test_to_transcript(variant):
    result = variant.to_rna(gene_name=True)
    assert len(result) == 1
    assert result[0].reference == "KIT"
    assert result[0].start == 976
    assert result[0].end == 977
    assert result[0].strand == "+"
    assert result[0].ref_seq == "TA"
    assert result[0].alt_seq == "AG"


# class TestRnaDelinsToSubstitution:
"""Test mapping RNA delins that collapse to substitutions."""


@pytest.fixture
def variant():
    """
    Since the 'A' in the 2nd position of alt allele matches the 'A' in the 2nd position of the
    reference allele, this is effectively a substitution because only the first base changes.
    """
    return RnaDelins(
        reference="ENST00000288135",
        reference_type="transcript",
        start=976,
        end=977,
        ref_seq="TA",
        alt_seq="AA",
    )


def test_hgvs(variant):
    assert variant.hgvs() == "ENST00000288135:r.976_977delinsAA"


def test_to_cds(variant):
    result = variant.to_cds()
    assert len(result) == 1
    assert isinstance(result[0], CdsSubstitution)
    assert result[0].reference == "ENST00000288135"
    assert result[0].start == 879
    assert result[0].strand == "+"
    assert result[0].ref_seq == "T"
    assert result[0].alt_seq == "A"


def test_to_contig(variant):
    result = variant.to_contig()
    assert len(result) == 1
    assert isinstance(result[0], DnaSubstitution)
    assert result[0].reference == "4"
    assert result[0].start == 55570012
    assert result[0].strand == "+"
    assert result[0].ref_seq == "T"
    assert result[0].alt_seq == "A"


def test_to_gene(variant):
    result = variant.to_gene()
    assert len(result) == 1
    assert isinstance(result[0], DnaSubstitution)
    assert result[0].reference == "KIT"
    assert result[0].start == 55570012
    assert result[0].strand == "+"
    assert result[0].ref_seq == "T"
    assert result[0].alt_seq == "A"


def test_to_protein(variant):
    result = variant.to_protein()
    assert len(result) == 1
    assert isinstance(result[0], ProteinSubstitution)
    assert result[0].reference == "ENSP00000288135"
    assert result[0].start == 293
    assert result[0].strand == "+"
    assert result[0].ref_seq == "N"
    assert result[0].alt_seq == "K"


def test_to_transcript(variant):
    result = variant.to_rna()
    assert len(result) == 1
    assert isinstance(result[0], RnaSubstitution)
    assert result[0].reference == "ENST00000288135"
    assert result[0].start == 976
    assert result[0].strand == "+"
    assert result[0].ref_seq == "T"
    assert result[0].alt_seq == "A"


# class TestCdsDuplication:
@pytest.fixture
def variant():
    return CdsDuplication(
        reference="ENST00000269305", reference_type="transcript", start=286, end=288, ref_seq="TCT"
    )


def test_hgvs(variant):
    assert variant.hgvs() == "ENST00000269305:c.286_288dup"


def test_to_cds(variant):
    result = variant.to_cds()
    assert len(result) == 1
    assert result[0].reference == "ENST00000269305"
    assert result[0].start == 286
    assert result[0].end == 288
    assert result[0].strand == "-"
    assert result[0].ref_seq == "TCT"


def test_to_contig(variant):
    result = variant.to_contig()
    assert len(result) == 1
    assert result[0].reference == "17"
    assert result[0].start == 7579399
    assert result[0].end == 7579401
    assert result[0].strand == "+"
    assert result[0].ref_seq == "AGA"


def test_to_gene(variant):
    result = variant.to_gene()
    assert len(result) == 1
    assert result[0].reference == "TP53"
    assert result[0].start == 7579399
    assert result[0].end == 7579401
    assert result[0].strand == "-"
    assert result[0].ref_seq == "TCT"


def test_to_protein(variant):
    result = variant.to_protein()
    assert len(result) == 1
    assert result[0].reference == "ENSP00000269305"
    assert result[0].start == 96
    assert result[0].end == 96
    assert result[0].strand == "-"
    assert result[0].ref_seq == "S"


def test_to_transcript(variant):
    result = variant.to_rna()
    assert len(result) == 1
    assert result[0].reference == "ENST00000269305"
    assert result[0].start == 476
    assert result[0].end == 478
    assert result[0].strand == "-"
    assert result[0].ref_seq == "TCT"


def test_values(variant):
    assert variant.values() == {
        "reference": "ENST00000269305",
        "reference_type": "transcript",
        "start": 286,
        "end": 288,
        "alt_seq": "TCTTCT",
        "ref_seq": "TCT",
        "strand": "-",
        "start_offset": 0,
        "start_seq": "T",
        "end_offset": 0,
        "end_seq": "T",
        "user_repr": None,
        "graphkb_id": None,
    }


# class TestCdsDuplicationGeneName:
@pytest.fixture
def variant():
    return CdsDuplication(
        reference="ENST00000269305", reference_type="transcript", start=286, end=288, ref_seq="TCT"
    )


def test_to_cds(variant):
    result = variant.to_cds(gene_name=True)
    assert len(result) == 1
    assert result[0].reference == "TP53"
    assert result[0].start == 286
    assert result[0].end == 288
    assert result[0].strand == "-"
    assert result[0].ref_seq == "TCT"


def test_to_contig(variant):
    result = variant.to_contig(gene_name=True)
    assert len(result) == 1
    assert result[0].reference == "TP53"
    assert result[0].start == 7579399
    assert result[0].end == 7579401
    assert result[0].strand == "+"
    assert result[0].ref_seq == "AGA"


def test_to_gene(variant):
    result = variant.to_gene(gene_name=True)
    assert len(result) == 1
    assert result[0].reference == "TP53"
    assert result[0].start == 7579399
    assert result[0].end == 7579401
    assert result[0].strand == "-"
    assert result[0].ref_seq == "TCT"


def test_to_protein(variant):
    result = variant.to_protein(gene_name=True)
    assert len(result) == 1
    assert result[0].reference == "TP53"
    assert result[0].start == 96
    assert result[0].end == 96
    assert result[0].strand == "-"
    assert result[0].ref_seq == "S"


def test_to_transcript(variant):
    result = variant.to_rna(gene_name=True)
    assert len(result) == 1
    assert result[0].reference == "TP53"
    assert result[0].start == 476
    assert result[0].end == 478
    assert result[0].strand == "-"
    assert result[0].ref_seq == "TCT"


# class TestCdsDuplicationNoRefSeq:
# VHL:c.165_166dup
@pytest.fixture
def variant():
    return CdsDuplication(reference="VHL", reference_type="gene", start=165, end=166)


def test_to_cds(variant):
    result = variant.to_cds()
    assert len(result) == 2
    assert result[0].reference == "ENST00000256474"
    assert result[0].start == 165
    assert result[0].end == 166
    assert result[0].strand == "+"
    assert result[0].ref_seq == "GG"


# class TestCdsDuplicationOffsetSingle:
@pytest.fixture
def variant():
    return CdsDuplication(
        reference="ENST00000323250",
        reference_type="transcript",
        start=1,
        start_offset=-86,
        ref_seq="G",
    )


def test_hgvs(variant):
    assert variant.hgvs() == "ENST00000323250:c.-86dup"


def test_to_contig(variant):
    result = variant.to_contig()
    assert result[0].reference == "18"
    assert result[0].start == 657657
    assert result[0].end == 657657
    assert result[0].strand == "+"
    assert result[0].ref_seq == "G"


# class TestCdsDuplicationOffsetMultiple:
@pytest.fixture
def variant():
    return CdsDuplication(
        reference="ENST00000323250",
        reference_type="transcript",
        start=1,
        start_offset=-86,
        end=1,
        end_offset=-84,
        ref_seq="GGC",
    )


def test_hgvs(variant):
    assert variant.hgvs() == "ENST00000323250:c.-86_-84dup"


def test_to_contig(variant):
    result = variant.to_contig()
    assert len(result) == 1
    assert result[0].reference == "18"
    assert result[0].start == 657657
    assert result[0].end == 657659
    assert result[0].strand == "+"
    assert result[0].ref_seq == "GGC"


# class TestDnaDuplication:
@pytest.fixture
def variant():
    return DnaDuplication(
        reference="KIT", reference_type="gene", start=55570013, end=55570015, ref_seq="AAT"
    )


def test_hgvs(variant):
    assert variant.hgvs() == "KIT:g.55570013_55570015dup"


def test_to_cds(variant):
    result = variant.to_cds()
    assert len(result) == 2
    assert result[0].reference == "ENST00000288135"
    assert result[0].start == 880
    assert result[0].end == 882
    assert result[0].strand == "+"
    assert result[0].ref_seq == "AAT"


def test_to_contig(variant):
    result = variant.to_contig()
    assert len(result) == 1
    assert result[0].reference == "4"
    assert result[0].start == 55570013
    assert result[0].end == 55570015
    assert result[0].strand == "+"
    assert result[0].ref_seq == "AAT"


def test_to_gene(variant):
    result = variant.to_gene()
    assert len(result) == 1
    assert result[0].reference == "KIT"
    assert result[0].start == 55570013
    assert result[0].end == 55570015
    assert result[0].strand == "+"
    assert result[0].ref_seq == "AAT"


def test_to_protein(variant):
    result = variant.to_protein()
    assert len(result) == 2
    assert result[0].reference == "ENSP00000288135"
    assert result[0].start == 294
    assert result[0].end == 294
    assert result[0].strand == "+"
    assert result[0].ref_seq == "N"


def test_to_transcript(variant):
    result = variant.to_rna()
    assert len(result) == 2
    assert result[0].reference == "ENST00000288135"
    assert result[0].start == 977
    assert result[0].end == 979
    assert result[0].strand == "+"
    assert result[0].ref_seq == "AAT"


# class TestDnaDuplicationGeneName:
@pytest.fixture
def variant():
    return DnaDuplication(
        reference="KIT", reference_type="gene", start=55570013, end=55570015, ref_seq="AAT"
    )


def test_to_cds(variant):
    result = variant.to_cds(gene_name=True)
    assert len(result) == 2
    assert result[0].reference == "KIT"
    assert result[0].start == 880
    assert result[0].end == 882
    assert result[0].strand == "+"
    assert result[0].ref_seq == "AAT"


def test_to_contig(variant):
    result = variant.to_contig(gene_name=True)
    assert len(result) == 1
    assert result[0].reference == "KIT"
    assert result[0].start == 55570013
    assert result[0].end == 55570015
    assert result[0].strand == "+"
    assert result[0].ref_seq == "AAT"


def test_to_gene(variant):
    result = variant.to_gene(gene_name=True)
    assert len(result) == 1
    assert result[0].reference == "KIT"
    assert result[0].start == 55570013
    assert result[0].end == 55570015
    assert result[0].strand == "+"
    assert result[0].ref_seq == "AAT"


def test_to_protein(variant):
    result = variant.to_protein(gene_name=True)
    assert len(result) == 1
    assert result[0].reference == "KIT"
    assert result[0].start == 294
    assert result[0].end == 294
    assert result[0].strand == "+"
    assert result[0].ref_seq == "N"


def test_to_transcript(variant):
    result = variant.to_rna(gene_name=True)
    assert len(result) == 2
    assert result[0].reference == "KIT"
    assert result[0].start == 977
    assert result[0].end == 979
    assert result[0].strand == "+"
    assert result[0].ref_seq == "AAT"


# class TestProteinDuplication:
@pytest.fixture
def variant():
    return ProteinDuplication(
        reference="KIT",
        reference_type="gene",
        start=501,
        end=502,
        start_seq="S",
        end_seq="A",
        ref_seq="SA",
    )


def test_hgvs(variant):
    assert variant.hgvs() == "KIT:p.S501_A502dup"


def test_to_cds(variant):
    result = variant.to_cds()
    assert len(result) == 2
    assert result[0].reference == "ENST00000288135"
    assert result[0].start == 1501
    assert result[0].end == 1506
    assert result[0].strand == "+"
    assert result[0].ref_seq == "TCTGCC"


def test_to_contig(variant):
    result = variant.to_contig()
    assert len(result) == 1
    assert result[0].reference == "4"
    assert result[0].start == 55592177
    assert result[0].end == 55592182
    assert result[0].strand == "+"
    assert result[0].ref_seq == "TCTGCC"


def test_to_gene(variant):
    result = variant.to_gene()
    assert len(result) == 1
    assert result[0].reference == "KIT"
    assert result[0].start == 55592177
    assert result[0].end == 55592182
    assert result[0].strand == "+"
    assert result[0].ref_seq == "TCTGCC"


def test_to_protein(variant):
    result = variant.to_protein()
    assert len(result) == 2
    assert result[0].reference == "ENSP00000288135"
    assert result[0].start == 501
    assert result[0].end == 502
    assert result[0].strand == "+"
    assert result[0].ref_seq == "SA"


def test_to_transcript(variant):
    result = variant.to_rna()
    assert len(result) == 2
    assert result[0].reference == "ENST00000288135"
    assert result[0].start == 1598
    assert result[0].end == 1603
    assert result[0].strand == "+"
    assert result[0].ref_seq == "TCTGCC"


# class TestProteinDuplicationGeneName:
@pytest.fixture
def variant():
    return ProteinDuplication(
        reference="KIT",
        reference_type="gene",
        start=501,
        end=502,
        start_seq="S",
        end_seq="A",
        ref_seq="SA",
    )


def test_to_cds(variant):
    result = variant.to_cds(gene_name=True)
    assert len(result) == 2
    assert result[0].reference == "KIT"
    assert result[0].start == 1501
    assert result[0].end == 1506
    assert result[0].strand == "+"
    assert result[0].ref_seq == "TCTGCC"


def test_to_contig(variant):
    result = variant.to_contig(gene_name=True)
    assert len(result) == 1
    assert result[0].reference == "KIT"
    assert result[0].start == 55592177
    assert result[0].end == 55592182
    assert result[0].strand == "+"
    assert result[0].ref_seq == "TCTGCC"


def test_to_gene(variant):
    result = variant.to_gene(gene_name=True)
    assert len(result) == 1
    assert result[0].reference == "KIT"
    assert result[0].start == 55592177
    assert result[0].end == 55592182
    assert result[0].strand == "+"
    assert result[0].ref_seq == "TCTGCC"


def test_to_protein(variant):
    result = variant.to_protein(gene_name=True)
    assert len(result) == 2
    assert result[0].reference == "KIT"
    assert result[0].start == 501
    assert result[0].end == 502
    assert result[0].strand == "+"
    assert result[0].ref_seq == "SA"


def test_to_transcript(variant):
    result = variant.to_rna(gene_name=True)
    assert len(result) == 2
    assert result[0].reference == "KIT"
    assert result[0].start == 1598
    assert result[0].end == 1603
    assert result[0].strand == "+"
    assert result[0].ref_seq == "TCTGCC"


# class TestRnaDuplication:
@pytest.fixture
def variant():
    return RnaDuplication(
        reference="ENST00000288135", reference_type="transcript", start=977, end=979, ref_seq="AAT"
    )


def test_hgvs(variant):
    assert variant.hgvs() == "ENST00000288135:r.977_979dup"


def test_to_cds(variant):
    result = variant.to_cds()
    assert len(result) == 1
    assert result[0].reference == "ENST00000288135"
    assert result[0].start == 880
    assert result[0].end == 882
    assert result[0].strand == "+"
    assert result[0].ref_seq == "AAT"


def test_to_contig(variant):
    result = variant.to_contig()
    assert len(result) == 1
    assert result[0].reference == "4"
    assert result[0].start == 55570013
    assert result[0].end == 55570015
    assert result[0].strand == "+"
    assert result[0].ref_seq == "AAT"


def test_to_gene(variant):
    result = variant.to_gene()
    assert len(result) == 1
    assert result[0].reference == "KIT"
    assert result[0].start == 55570013
    assert result[0].end == 55570015
    assert result[0].strand == "+"
    assert result[0].ref_seq == "AAT"


def test_to_protein(variant):
    result = variant.to_protein()
    assert len(result) == 1
    assert result[0].reference == "ENSP00000288135"
    assert result[0].start == 294
    assert result[0].end == 294
    assert result[0].strand == "+"
    assert result[0].ref_seq == "N"


def test_to_transcript(variant):
    result = variant.to_rna()
    assert len(result) == 1
    assert result[0].reference == "ENST00000288135"
    assert result[0].start == 977
    assert result[0].end == 979
    assert result[0].strand == "+"
    assert result[0].ref_seq == "AAT"


# class TestRnaDuplicationGeneName:
@pytest.fixture
def variant():
    return RnaDuplication(
        reference="ENST00000288135", reference_type="transcript", start=977, end=979, ref_seq="AAT"
    )


def test_to_cds(variant):
    result = variant.to_cds(gene_name=True)
    assert len(result) == 1
    assert result[0].reference == "KIT"
    assert result[0].start == 880
    assert result[0].end == 882
    assert result[0].strand == "+"
    assert result[0].ref_seq == "AAT"


def test_to_contig(variant):
    result = variant.to_contig(gene_name=True)
    assert len(result) == 1
    assert result[0].reference == "KIT"
    assert result[0].start == 55570013
    assert result[0].end == 55570015
    assert result[0].strand == "+"
    assert result[0].ref_seq == "AAT"


def test_to_gene(variant):
    result = variant.to_gene(gene_name=True)
    assert len(result) == 1
    assert result[0].reference == "KIT"
    assert result[0].start == 55570013
    assert result[0].end == 55570015
    assert result[0].strand == "+"
    assert result[0].ref_seq == "AAT"


def test_to_protein(variant):
    result = variant.to_protein(gene_name=True)
    assert len(result) == 1
    assert result[0].reference == "KIT"
    assert result[0].start == 294
    assert result[0].end == 294
    assert result[0].strand == "+"
    assert result[0].ref_seq == "N"


def test_to_transcript(variant):
    result = variant.to_rna(gene_name=True)
    assert len(result) == 1
    assert result[0].reference == "KIT"
    assert result[0].start == 977
    assert result[0].end == 979
    assert result[0].strand == "+"
    assert result[0].ref_seq == "AAT"


# class TestCdsFusion:
@pytest.fixture
def variant():
    return CdsFusion(
        break1_reference="ENST00000315869",
        break1_reference_type="transcript",
        break1_start=534,
        break2_reference="ENST00000271526",
        break2_reference_type="transcript",
        break2_start=469,
    )


def test_hgvs(variant):
    assert variant.hgvs() == "(ENST00000315869,ENST00000271526):fusion(c.534,c.469)"


def test_to_cds(variant):
    result = variant.to_cds()
    assert len(result) == 1
    assert result[0].break1_reference == "ENST00000315869"
    assert result[0].break1_start == 534
    assert result[0].break2_reference == "ENST00000271526"
    assert result[0].break2_start == 469


def test_to_exon(variant):
    result = variant.to_exon()
    assert len(result) == 1
    assert result[0].break1_reference == "TFE3"
    assert result[0].break1_start == 3
    assert result[0].break2_reference == "PRCC"
    assert result[0].break2_start == 2


def test_to_contig(variant):
    result = variant.to_contig()
    assert len(result) == 1
    assert result[0].break1_reference == "X"
    assert result[0].break1_start == 48896632
    assert result[0].break2_reference == "1"
    assert result[0].break2_start == 156752074


def test_to_gene(variant):
    result = variant.to_gene()
    assert len(result) == 1
    assert result[0].break1_reference == "TFE3"
    assert result[0].break1_start == 48896632
    assert result[0].break2_reference == "PRCC"
    assert result[0].break2_start == 156752074


def test_to_transcript(variant):
    result = variant.to_rna()
    assert len(result) == 1
    assert result[0].break1_reference == "ENST00000315869"
    assert result[0].break1_start == 794
    assert result[0].break2_reference == "ENST00000271526"
    assert result[0].break2_start == 741


# class TestCdsFusionGeneName:
@pytest.fixture
def variant():
    return CdsFusion(
        break1_reference="ENST00000315869",
        break1_reference_type="transcript",
        break1_start=534,
        break2_reference="ENST00000271526",
        break2_reference_type="transcript",
        break2_start=469,
    )


def test_to_cds(variant):
    result = variant.to_cds(gene_name=True)
    assert len(result) == 1
    assert result[0].break1_reference == "TFE3"
    assert result[0].break1_start == 534
    assert result[0].break2_reference == "PRCC"
    assert result[0].break2_start == 469


def test_to_exon(variant):
    result = variant.to_exon(gene_name=True)
    assert len(result) == 1
    assert result[0].break1_reference == "TFE3"
    assert result[0].break1_start == 3
    assert result[0].break2_reference == "PRCC"
    assert result[0].break2_start == 2


def test_to_contig(variant):
    result = variant.to_contig(gene_name=True)
    assert len(result) == 1
    assert result[0].break1_reference == "TFE3"
    assert result[0].break1_start == 48896632
    assert result[0].break2_reference == "PRCC"
    assert result[0].break2_start == 156752074


def test_to_gene(variant):
    result = variant.to_gene(gene_name=True)
    assert len(result) == 1
    assert result[0].break1_reference == "TFE3"
    assert result[0].break1_start == 48896632
    assert result[0].break2_reference == "PRCC"
    assert result[0].break2_start == 156752074


def test_to_transcript(variant):
    result = variant.to_rna(gene_name=True)
    assert len(result) == 1
    assert result[0].break1_reference == "TFE3"
    assert result[0].break1_start == 794
    assert result[0].break2_reference == "PRCC"
    assert result[0].break2_start == 741


# class TestDnaFusion:
@pytest.fixture
def variant():
    return DnaFusion(
        break1_reference="RYK",
        break1_reference_type="gene",
        break1_start=133921574,
        break2_reference="RABL3",
        break2_reference_type="gene",
        break2_start=120417420,
    )


def test_hgvs(variant):
    assert variant.hgvs() == "(RYK,RABL3):fusion(g.133921574,g.120417420)"


def test_to_cds(variant):
    result = variant.to_cds()
    assert len(result) == 2
    assert result[0].break1_reference == "ENST00000427044"
    assert result[0].break1_start == 212
    assert result[0].break2_reference == "ENST00000273375"
    assert result[0].break2_start == 384


def test_to_exon(variant):
    result = variant.to_exon()
    assert len(result) == 9
    assert result[0].break1_reference == "RYK"
    assert result[0].break1_start == 7
    assert result[0].break2_reference == "RABL3"
    assert result[0].break2_start == 4


def test_to_contig(variant):
    result = variant.to_contig()
    assert len(result) == 1
    assert result[0].break1_reference == "3"
    assert result[0].break1_start == 133921574
    assert result[0].break2_reference == "3"
    assert result[0].break2_start == 120417420


def test_to_gene(variant):
    result = variant.to_gene()
    assert len(result) == 1
    assert result[0].break1_reference == "RYK"
    assert result[0].break1_start == 133921574
    assert result[0].break2_reference == "RABL3"
    assert result[0].break2_start == 120417420


def test_to_transcript(variant):
    result = variant.to_rna()
    assert len(result) == 20
    assert result[0].break1_reference == "ENST00000296084"
    assert result[0].break1_start == 782
    assert result[0].break2_reference == "ENST00000273375"
    assert result[0].break2_start == 414


# class TestDnaFusionGeneName:
@pytest.fixture
def variant():
    return DnaFusion(
        break1_reference="RYK",
        break1_reference_type="gene",
        break1_start=133921574,
        break2_reference="RABL3",
        break2_reference_type="gene",
        break2_start=120417420,
    )


def test_to_cds(variant):
    result = variant.to_cds(gene_name=True)
    assert len(result) == 2
    assert result[0].break1_reference == "RYK"
    assert result[0].break1_start == 212
    assert result[0].break2_reference == "RABL3"
    assert result[0].break2_start == 384


def test_to_exon(variant):
    result = variant.to_exon(gene_name=True)
    assert len(result) == 9
    assert result[0].break1_reference == "RYK"
    assert result[0].break1_start == 7
    assert result[0].break2_reference == "RABL3"
    assert result[0].break2_start == 4


def test_to_contig(variant):
    result = variant.to_contig(gene_name=True)
    assert len(result) == 1
    assert result[0].break1_reference == "RYK"
    assert result[0].break1_start == 133921574
    assert result[0].break2_reference == "RABL3"
    assert result[0].break2_start == 120417420


def test_to_gene(variant):
    result = variant.to_gene(gene_name=True)
    assert len(result) == 1
    assert result[0].break1_reference == "RYK"
    assert result[0].break1_start == 133921574
    assert result[0].break2_reference == "RABL3"
    assert result[0].break2_start == 120417420


def test_to_transcript(variant):
    result = variant.to_rna(gene_name=True)
    assert len(result) == 20
    assert result[0].break1_reference == "RYK"
    assert result[0].break1_start == 782
    assert result[0].break2_reference == "RABL3"
    assert result[0].break2_start == 414


# class TestExonFusion:
@pytest.fixture
def variant():
    return ExonFusion(
        break1_reference="FGFR3",
        break1_reference_type="gene",
        break1_start=16,
        break2_reference="TACC3",
        break2_reference_type="gene",
        break2_start=8,
    )


def test_hgvs(variant):
    assert variant.hgvs() == "(FGFR3,TACC3):fusion(e.16,e.8)"


def test_to_cds(variant):
    result = variant.to_cds()
    assert len(result) == 18
    assert result[0].break1_reference == "ENST00000260795"
    assert result[0].break1_start == 2031
    assert result[0].break2_reference == "ENST00000313288"
    assert result[0].break2_start == 1645


def test_to_exon(variant):
    result = variant.to_exon()
    assert len(result) == 1
    assert result[0].break1_reference == "FGFR3"
    assert result[0].break1_start == 16
    assert result[0].break2_reference == "TACC3"
    assert result[0].break2_start == 8


def test_to_contig(variant):
    result = variant.to_contig()
    assert len(result) == 3
    assert result[0].break1_reference == "4"
    assert result[0].break1_start == 1808273
    assert result[0].break2_reference == "4"
    assert result[0].break2_start == 1737458


def test_to_gene(variant):
    result = variant.to_gene()
    assert len(result) == 3
    assert result[0].break1_reference == "FGFR3"
    assert result[0].break1_start == 1808273
    assert result[0].break2_reference == "TACC3"
    assert result[0].break2_start == 1737458


def test_to_transcript(variant):
    result = variant.to_rna()
    assert len(result) == 10
    assert result[0].break1_reference == "ENST00000260795"
    assert result[0].break1_start == 2271
    assert result[0].break2_reference == "ENST00000313288"
    assert result[0].break2_start == 1751


# class TestExonFusionGeneName:
@pytest.fixture
def variant():
    return ExonFusion(
        break1_reference="FGFR3",
        break1_reference_type="gene",
        break1_start=16,
        break2_reference="TACC3",
        break2_reference_type="gene",
        break2_start=8,
    )


def test_to_cds(variant):
    result = variant.to_cds(gene_name=True)
    assert len(result) == 18
    assert result[0].break1_reference == "FGFR3"
    assert result[0].break1_start == 2031
    assert result[0].break2_reference == "TACC3"
    assert result[0].break2_start == 1645


def test_to_exon(variant):
    result = variant.to_exon(gene_name=True)
    assert len(result) == 1
    assert result[0].break1_reference == "FGFR3"
    assert result[0].break1_start == 16
    assert result[0].break2_reference == "TACC3"
    assert result[0].break2_start == 8


def test_to_contig(variant):
    result = variant.to_contig(gene_name=True)
    assert len(result) == 3
    assert result[0].break1_reference == "FGFR3"
    assert result[0].break1_start == 1808273
    assert result[0].break2_reference == "TACC3"
    assert result[0].break2_start == 1737458


def test_to_gene(variant):
    result = variant.to_gene(gene_name=True)
    assert len(result) == 3
    assert result[0].break1_reference == "FGFR3"
    assert result[0].break1_start == 1808273
    assert result[0].break2_reference == "TACC3"
    assert result[0].break2_start == 1737458


def test_to_transcript(variant):
    result = variant.to_rna(gene_name=True)
    assert len(result) == 10
    assert result[0].break1_reference == "FGFR3"
    assert result[0].break1_start == 2271
    assert result[0].break2_reference == "TACC3"
    assert result[0].break2_start == 1751


# class TestRnaFusion:
@pytest.fixture
def variant():
    return RnaFusion(
        break1_reference="TFE3",
        break1_reference_type="gene",
        break1_start=794,
        break2_reference="PRCC",
        break2_reference_type="gene",
        break2_start=741,
    )


def test_hgvs(variant):
    assert variant.hgvs() == "(TFE3,PRCC):fusion(r.794,r.741)"


def test_to_cds(variant):
    result = variant.to_cds()
    assert len(result) == 2
    assert result[0].break1_reference == "ENST00000315869"
    assert result[0].break1_start == 534
    assert result[0].break2_reference == "ENST00000271526"
    assert result[0].break2_start == 469


def test_to_exon(variant):
    result = variant.to_exon()
    assert len(result) == 15
    assert result[0].break1_reference == "TFE3"
    assert result[0].break1_start == 3
    assert result[0].break2_reference == "PRCC"
    assert result[0].break2_start == 6


def test_to_contig(variant):
    result = variant.to_contig()
    assert len(result) == 20
    assert result[0].break1_reference == "X"
    assert result[0].break1_start == 48891281
    assert result[0].break2_reference == "1"
    assert result[0].break2_start == 156752074


def test_to_gene(variant):
    result = variant.to_gene()
    assert len(result) == 20
    assert result[0].break1_reference == "TFE3"
    assert result[0].break1_start == 48891281
    assert result[0].break2_reference == "PRCC"
    assert result[0].break2_start == 156752074


def test_to_transcript(variant):
    result = variant.to_rna()
    assert len(result) == 20
    assert result[0].break1_reference == "ENST00000315869"
    assert result[0].break1_start == 794
    assert result[0].break2_reference == "ENST00000271526"
    assert result[0].break2_start == 741


# class TestRnaFusionGeneName:
@pytest.fixture
def variant():
    return RnaFusion(
        break1_reference="TFE3",
        break1_reference_type="gene",
        break1_start=794,
        break2_reference="PRCC",
        break2_reference_type="gene",
        break2_start=741,
    )


def test_to_cds(variant):
    result = variant.to_cds(gene_name=True)
    assert len(result) == 2
    assert result[0].break1_reference == "TFE3"
    assert result[0].break1_start == 534
    assert result[0].break2_reference == "PRCC"
    assert result[0].break2_start == 469


def test_to_exon(variant):
    result = variant.to_exon(gene_name=True)
    assert len(result) == 15
    assert result[0].break1_reference == "TFE3"
    assert result[0].break1_start == 3
    assert result[0].break2_reference == "PRCC"
    assert result[0].break2_start == 6


def test_to_contig(variant):
    result = variant.to_contig(gene_name=True)
    assert len(result) == 20
    assert result[0].break1_reference == "TFE3"
    assert result[0].break1_start == 48891281
    assert result[0].break2_reference == "PRCC"
    assert result[0].break2_start == 156752074


def test_to_gene(variant):
    result = variant.to_gene(gene_name=True)
    assert len(result) == 20
    assert result[0].break1_reference == "TFE3"
    assert result[0].break1_start == 48891281
    assert result[0].break2_reference == "PRCC"
    assert result[0].break2_start == 156752074


def test_to_transcript(variant):
    result = variant.to_rna(gene_name=True)
    assert len(result) == 20
    assert result[0].break1_reference == "TFE3"
    assert result[0].break1_start == 794
    assert result[0].break2_reference == "PRCC"
    assert result[0].break2_start == 741


# class TestCdsInsertion:
@pytest.fixture
def variant():
    return CdsInsertion(
        reference="ENST00000288135",
        reference_type="transcript",
        start=1672,
        end=1673,
        alt_seq="TTC",
    )


def test_hgvs(variant):
    assert variant.hgvs() == "ENST00000288135:c.1672_1673insTTC"


def test_to_cds(variant):
    result = variant.to_cds()
    assert len(result) == 1
    assert result[0].reference == "ENST00000288135"
    assert result[0].start == 1672
    assert result[0].end == 1673
    assert result[0].strand == "+"
    assert result[0].alt_seq == "TTC"


def test_to_contig(variant):
    result = variant.to_contig()
    assert len(result) == 1
    assert result[0].reference == "4"
    assert result[0].start == 55593606
    assert result[0].end == 55593607
    assert result[0].strand == "+"
    assert result[0].alt_seq == "TTC"


def test_to_gene(variant):
    result = variant.to_gene()
    assert len(result) == 1
    assert result[0].reference == "KIT"
    assert result[0].start == 55593606
    assert result[0].end == 55593607
    assert result[0].strand == "+"
    assert result[0].alt_seq == "TTC"


def test_to_protein(variant):
    result = variant.to_protein()
    assert len(result) == 1
    assert result[0].reference == "ENSP00000288135"
    assert result[0].start == 558
    assert result[0].end == 558
    assert result[0].strand == "+"
    assert result[0].alt_seq == "F"


def test_to_transcript(variant):
    result = variant.to_rna()
    assert len(result) == 1
    assert result[0].reference == "ENST00000288135"
    assert result[0].start == 1769
    assert result[0].end == 1770
    assert result[0].strand == "+"
    assert result[0].alt_seq == "TTC"


def test_values(variant):
    assert variant.values() == {
        "reference": "ENST00000288135",
        "reference_type": "transcript",
        "start": 1672,
        "end": 1673,
        "alt_seq": "TTC",
        "ref_seq": "NN",
        "strand": "+",
        "start_offset": 0,
        "start_seq": "N",
        "end_offset": 0,
        "end_seq": "N",
        "user_repr": None,
        "graphkb_id": None,
    }


# class TestCdsInsertionGeneName:
@pytest.fixture
def variant():
    return CdsInsertion(
        reference="ENST00000288135",
        reference_type="transcript",
        start=1672,
        end=1673,
        alt_seq="TTC",
    )


def test_to_cds(variant):
    result = variant.to_cds(gene_name=True)
    assert len(result) == 1
    assert result[0].reference == "KIT"
    assert result[0].start == 1672
    assert result[0].end == 1673
    assert result[0].strand == "+"
    assert result[0].alt_seq == "TTC"


def test_to_contig(variant):
    result = variant.to_contig(gene_name=True)
    assert len(result) == 1
    assert result[0].reference == "KIT"
    assert result[0].start == 55593606
    assert result[0].end == 55593607
    assert result[0].strand == "+"
    assert result[0].alt_seq == "TTC"


def test_to_gene(variant):
    result = variant.to_gene(gene_name=True)
    assert len(result) == 1
    assert result[0].reference == "KIT"
    assert result[0].start == 55593606
    assert result[0].end == 55593607
    assert result[0].strand == "+"
    assert result[0].alt_seq == "TTC"


def test_to_protein(variant):
    result = variant.to_protein(gene_name=True)
    assert len(result) == 1
    assert result[0].reference == "KIT"
    assert result[0].start == 558
    assert result[0].end == 558
    assert result[0].strand == "+"
    assert result[0].alt_seq == "F"


def test_to_transcript(variant):
    result = variant.to_rna(gene_name=True)
    assert len(result) == 1
    assert result[0].reference == "KIT"
    assert result[0].start == 1769
    assert result[0].end == 1770
    assert result[0].strand == "+"
    assert result[0].alt_seq == "TTC"


# class TestCdsInsertionOffsetSingle:
@pytest.fixture
def variant():
    return CdsInsertion(
        reference="ENST00000288135",
        reference_type="transcript",
        start=1,
        start_offset=-3,
        alt_seq="TTC",
    )


def test_hgvs(variant):
    assert variant.hgvs() == "ENST00000288135:c.-3_-2insTTC"


def test_to_contig(variant):
    result = variant.to_contig()
    assert len(result) == 1
    assert result[0].reference == "4"
    assert result[0].start == 55524179
    assert result[0].end == 55524180
    assert result[0].strand == "+"
    assert result[0].alt_seq == "TTC"


# class TestCdsInsertionOffsetMultiple:
@pytest.fixture
def variant():
    return CdsInsertion(
        reference="ENST00000288135",
        reference_type="transcript",
        start=1,
        start_offset=-3,
        end=2,
        end_offset=-3,
        alt_seq="TTC",
    )


def test_hgvs(variant):
    assert variant.hgvs() == "ENST00000288135:c.-3_-2insTTC"


def test_to_contig(variant):
    result = variant.to_contig()
    assert len(result) == 1
    assert result[0].reference == "4"
    assert result[0].start == 55524179
    assert result[0].end == 55524180
    assert result[0].strand == "+"
    assert result[0].alt_seq == "TTC"


# class TestCdsInsertionOffsetMultiple2:
@pytest.fixture
def variant():
    return CdsInsertion(
        reference="ENST00000288135",
        reference_type="transcript",
        start=1,
        start_offset=-1,
        end=1,
        end_offset=0,
        alt_seq="TTC",
    )


def test_hgvs(variant):
    assert variant.hgvs() == "ENST00000288135:c.-1_1insTTC"


def test_to_contig(variant):
    result = variant.to_contig()
    assert len(result) == 1
    assert result[0].reference == "4"
    assert result[0].start == 55524181
    assert result[0].end == 55524182
    assert result[0].strand == "+"
    assert result[0].alt_seq == "TTC"


# class TestCdsInsertionToDuplication:
@pytest.fixture
def variant():
    return CdsInsertion(
        reference="ENST00000288135", reference_type="transcript", start=1672, end=1673, alt_seq="A"
    )


def test_to_cds(variant):
    result = variant.to_cds()
    assert len(result) == 1
    assert type(result[0]) == CdsDuplication
    assert result[0].reference == "ENST00000288135"
    assert result[0].start == 1672
    assert result[0].end == 1673
    assert result[0].strand == "+"
    assert result[0].ref_seq == "AA"


# class TestDnaInsertion:
@pytest.fixture
def variant():
    return DnaInsertion(
        reference="17", reference_type="contig", start=7579470, end=7579471, alt_seq="TTT"
    )


def test_hgvs(variant):
    assert variant.hgvs() == "17:g.7579470_7579471insTTT"


def test_to_cds(variant):
    result = variant.to_cds()
    assert len(result) == 6
    assert result[0].reference == "ENST00000269305"
    assert result[0].start == 216
    assert result[0].end == 217
    assert result[0].strand == "-"
    assert result[0].alt_seq == "AAA"


def test_to_contig(variant):
    result = variant.to_contig()
    assert len(result) == 1
    assert result[0].reference == "17"
    assert result[0].start == 7579470
    assert result[0].end == 7579471
    assert result[0].strand == "+"
    assert result[0].alt_seq == "TTT"


def test_to_gene(variant):
    result = variant.to_gene()
    assert len(result) == 1
    assert result[0].reference == "TP53"
    assert result[0].start == 7579470
    assert result[0].end == 7579471
    assert result[0].strand == "-"
    assert result[0].alt_seq == "AAA"


def test_to_protein(variant):
    result = variant.to_protein()
    assert len(result) == 6
    assert result[0].reference == "ENSP00000269305"
    assert result[0].start == 72
    assert result[0].end == 73
    assert result[0].strand == "-"
    assert result[0].alt_seq == "K"  # AAA


def test_to_transcript(variant):
    result = variant.to_rna()
    assert len(result) == 9
    assert result[0].reference == "ENST00000269305"
    assert result[0].start == 406
    assert result[0].end == 407
    assert result[0].strand == "-"
    assert result[0].alt_seq == "AAA"


# class TestDnaInsertionGeneName:
@pytest.fixture
def variant():
    return DnaInsertion(
        reference="17", reference_type="contig", start=7579470, end=7579471, alt_seq="TTT"
    )


def test_to_cds(variant):
    result = variant.to_cds(gene_name=True)
    assert len(result) == 6
    assert result[0].reference == "TP53"
    assert result[0].start == 216
    assert result[0].end == 217
    assert result[0].strand == "-"
    assert result[0].alt_seq == "AAA"


def test_to_contig(variant):
    result = variant.to_contig(gene_name=True)
    assert len(result) == 1
    assert result[0].reference == "TP53"
    assert result[0].start == 7579470
    assert result[0].end == 7579471
    assert result[0].strand == "+"
    assert result[0].alt_seq == "TTT"


def test_to_gene(variant):
    result = variant.to_gene(gene_name=True)
    assert len(result) == 1
    assert result[0].reference == "TP53"
    assert result[0].start == 7579470
    assert result[0].end == 7579471
    assert result[0].strand == "-"
    assert result[0].alt_seq == "AAA"


def test_to_protein(variant):
    result = variant.to_protein(gene_name=True)
    assert len(result) == 1
    assert result[0].reference == "TP53"
    assert result[0].start == 72
    assert result[0].end == 73
    assert result[0].strand == "-"
    assert result[0].alt_seq == "K"  # AAA


def test_to_transcript(variant):
    result = variant.to_rna(gene_name=True)
    assert len(result) == 9
    assert result[0].reference == "TP53"
    assert result[0].start == 406
    assert result[0].end == 407
    assert result[0].strand == "-"
    assert result[0].alt_seq == "AAA"


# class TestProteinInsertion:
@pytest.fixture
def variant():
    return ProteinInsertion(
        reference="KRAS",
        reference_type="gene",
        start=11,
        end=12,
        start_seq="A",
        end_seq="G",
        alt_seq="T",
    )


def test_hgvs(variant):
    assert variant.hgvs() == "KRAS:p.A11_G12insT"


def test_to_cds(variant):
    result = variant.to_cds()
    assert len(result) == 16
    assert result[0].reference == "ENST00000256078"
    assert result[0].start == 33
    assert result[0].end == 34
    assert result[0].strand == "-"
    assert result[0].alt_seq == "ACA"


def test_to_contig(variant):
    result = variant.to_contig()
    assert len(result) == 4
    assert result[0].reference == "12"
    assert result[0].start == 25398285
    assert result[0].end == 25398286
    assert result[0].strand == "+"
    assert result[0].alt_seq == "TGT"


def test_to_gene(variant):
    result = variant.to_gene()
    assert len(result) == 4
    assert result[0].reference == "KRAS"
    assert result[0].start == 25398285
    assert result[0].end == 25398286
    assert result[0].strand == "-"
    assert result[0].alt_seq == "ACA"


def test_to_protein(variant):
    result = variant.to_protein()
    assert len(result) == 4
    assert result[0].reference == "ENSP00000256078"
    assert result[0].start == 11
    assert result[0].end == 12
    assert result[0].strand == "-"
    assert result[0].alt_seq == "T"


def test_to_transcript(variant):
    result = variant.to_rna()
    assert len(result) == 16
    assert result[0].reference == "ENST00000256078"
    assert result[0].start == 97
    assert result[0].end == 98
    assert result[0].strand == "-"
    assert result[0].alt_seq == "ACA"


# class TestProteinInsertionGeneName:
@pytest.fixture
def variant():
    return ProteinInsertion(
        reference="KRAS",
        reference_type="gene",
        start=11,
        end=12,
        start_seq="A",
        end_seq="G",
        alt_seq="T",
    )


def test_to_cds(variant):
    result = variant.to_cds(gene_name=True)
    assert len(result) == 16
    assert result[0].reference == "KRAS"
    assert result[0].start == 33
    assert result[0].end == 34
    assert result[0].strand == "-"
    assert result[0].alt_seq == "ACA"


def test_to_contig(variant):
    result = variant.to_contig(gene_name=True)
    assert len(result) == 4
    assert result[0].reference == "KRAS"
    assert result[0].start == 25398285
    assert result[0].end == 25398286
    assert result[0].strand == "+"
    assert result[0].alt_seq == "TGT"


def test_to_gene(variant):
    result = variant.to_gene(gene_name=True)
    assert len(result) == 4
    assert result[0].reference == "KRAS"
    assert result[0].start == 25398285
    assert result[0].end == 25398286
    assert result[0].strand == "-"
    assert result[0].alt_seq == "ACA"


def test_to_protein(variant):
    result = variant.to_protein(gene_name=True)
    assert len(result) == 4
    assert result[0].reference == "KRAS"
    assert result[0].start == 11
    assert result[0].end == 12
    assert result[0].strand == "-"
    assert result[0].alt_seq == "T"


def test_to_transcript(variant):
    result = variant.to_rna(gene_name=True)
    assert len(result) == 16
    assert result[0].reference == "KRAS"
    assert result[0].start == 97
    assert result[0].end == 98
    assert result[0].strand == "-"
    assert result[0].alt_seq == "ACA"


# class TestProteinInsertionRefSeq:
# KIT:p.Y503_F504insAY #157:5434
@pytest.fixture
def variant():
    return ProteinInsertion(
        reference="KIT", reference_type="gene", start=503, end=504, alt_seq="AY"
    )


def test_to_transcript(variant):
    result = variant.to_rna()
    assert len(result) == 16
    assert result[0].reference == "ENST00000288135"
    assert result[0].is_insertion
    assert result[0].start == 1606
    assert result[0].end == 1607
    assert result[0].strand == "+"
    assert result[0].ref_seq == "TT"
    assert result[0].alt_seq == "GCATAC"


# class TestRnaInsertion:
@pytest.fixture
def variant():
    return RnaInsertion(
        reference="ENST00000288135",
        reference_type="transcript",
        start=1769,
        end=1770,
        alt_seq="CAT",
    )


def test_hgvs(variant):
    assert variant.hgvs() == "ENST00000288135:r.1769_1770insCAT"


def test_to_cds(variant):
    result = variant.to_cds()
    assert len(result) == 1
    assert result[0].reference == "ENST00000288135"
    assert result[0].start == 1672
    assert result[0].end == 1673
    assert result[0].strand == "+"
    assert result[0].alt_seq == "CAT"


def test_to_contig(variant):
    result = variant.to_contig()
    assert len(result) == 1
    assert result[0].reference == "4"
    assert result[0].start == 55593606
    assert result[0].end == 55593607
    assert result[0].strand == "+"
    assert result[0].alt_seq == "CAT"


def test_to_gene(variant):
    result = variant.to_gene()
    assert len(result) == 1
    assert result[0].reference == "KIT"
    assert result[0].start == 55593606
    assert result[0].end == 55593607
    assert result[0].strand == "+"
    assert result[0].alt_seq == "CAT"


def test_to_protein(variant):
    result = variant.to_protein()
    assert len(result) == 1
    assert result[0].reference == "ENSP00000288135"
    assert result[0].start == 558
    assert result[0].end == 558
    assert result[0].strand == "+"
    assert result[0].alt_seq == "H"


def test_to_transcript(variant):
    result = variant.to_rna()
    assert len(result) == 1
    assert result[0].reference == "ENST00000288135"
    assert result[0].start == 1769
    assert result[0].end == 1770
    assert result[0].strand == "+"
    assert result[0].alt_seq == "CAT"


# class TestRnaInsertionGeneName:
@pytest.fixture
def variant():
    return RnaInsertion(
        reference="ENST00000288135",
        reference_type="transcript",
        start=1769,
        end=1770,
        alt_seq="CAT",
    )


def test_to_cds(variant):
    result = variant.to_cds(gene_name=True)
    assert len(result) == 1
    assert result[0].reference == "KIT"
    assert result[0].start == 1672
    assert result[0].end == 1673
    assert result[0].strand == "+"
    assert result[0].alt_seq == "CAT"


def test_to_contig(variant):
    result = variant.to_contig(gene_name=True)
    assert len(result) == 1
    assert result[0].reference == "KIT"
    assert result[0].start == 55593606
    assert result[0].end == 55593607
    assert result[0].strand == "+"
    assert result[0].alt_seq == "CAT"


def test_to_gene(variant):
    result = variant.to_gene(gene_name=True)
    assert len(result) == 1
    assert result[0].reference == "KIT"
    assert result[0].start == 55593606
    assert result[0].end == 55593607
    assert result[0].strand == "+"
    assert result[0].alt_seq == "CAT"


def test_to_protein(variant):
    result = variant.to_protein(gene_name=True)
    assert len(result) == 1
    assert result[0].reference == "KIT"
    assert result[0].start == 558
    assert result[0].end == 558
    assert result[0].strand == "+"
    assert result[0].alt_seq == "H"


def test_to_transcript(variant):
    result = variant.to_rna(gene_name=True)
    assert len(result) == 1
    assert result[0].reference == "KIT"
    assert result[0].start == 1769
    assert result[0].end == 1770
    assert result[0].strand == "+"
    assert result[0].alt_seq == "CAT"


# class TestCdsSubstitution:
@pytest.fixture
def variant():
    return CdsSubstitution(
        reference="ENST00000256078", reference_type="transcript", start=38, ref_seq="G", alt_seq="A"
    )


def test_hgvs(variant):
    assert variant.hgvs() == "ENST00000256078:c.38G>A"


def test_to_cds(variant):
    result = variant.to_cds()
    assert len(result) == 1
    assert result[0].reference == "ENST00000256078"
    assert result[0].start == 38
    assert result[0].strand == "-"
    assert result[0].ref_seq == "G"
    assert result[0].alt_seq == "A"


def test_to_contig(variant):
    result = variant.to_contig()
    assert len(result) == 1
    assert result[0].reference == "12"
    assert result[0].start == 25398281
    assert result[0].strand == "+"
    assert result[0].ref_seq == "C"
    assert result[0].alt_seq == "T"


def test_to_gene(variant):
    result = variant.to_gene()
    assert len(result) == 1
    assert result[0].reference == "KRAS"
    assert result[0].start == 25398281
    assert result[0].strand == "-"
    assert result[0].ref_seq == "G"
    assert result[0].alt_seq == "A"


def test_to_protein(variant):
    result = variant.to_protein()
    assert len(result) == 1
    assert result[0].reference == "ENSP00000256078"
    assert result[0].start == 13
    assert result[0].strand == "-"
    assert result[0].ref_seq == "G"  # GGC
    assert result[0].alt_seq == "D"  # GAC


def test_to_transcript(variant):
    result = variant.to_rna()
    assert len(result) == 1
    assert result[0].reference == "ENST00000256078"
    assert result[0].start == 102
    assert result[0].strand == "-"
    assert result[0].ref_seq == "G"
    assert result[0].alt_seq == "A"


def test_values(variant):
    assert variant.values() == {
        "reference": "ENST00000256078",
        "reference_type": "transcript",
        "start": 38,
        "end": 38,
        "alt_seq": "A",
        "ref_seq": "G",
        "strand": "-",
        "start_offset": 0,
        "start_seq": "G",
        "end_offset": 0,
        "end_seq": "G",
        "user_repr": None,
        "graphkb_id": None,
    }


# class TestCdsSubstitutionError:
def test_to_protein_ref_mismatch_error(variant):
    return CdsSubstitution(
        reference="ENST00000256078", reference_type="transcript", start=38, ref_seq="N", alt_seq="A"
    )
    with pytest.raises(RefMismatchError):
        variant.to_protein()


# class TestCdsSubstitutionGeneName:
@pytest.fixture
def variant():
    return CdsSubstitution(
        reference="ENST00000256078", reference_type="transcript", start=38, ref_seq="G", alt_seq="A"
    )


def test_to_cds(variant):
    result = variant.to_cds(gene_name=True)
    assert len(result) == 1
    assert result[0].reference == "KRAS"
    assert result[0].start == 38
    assert result[0].strand == "-"
    assert result[0].ref_seq == "G"
    assert result[0].alt_seq == "A"
    assert result[0].hgvs() == "KRAS:c.38G>A"


def test_to_contig(variant):
    result = variant.to_contig(gene_name=True)
    assert len(result) == 1
    assert result[0].reference == "KRAS"
    assert result[0].start == 25398281
    assert result[0].strand == "+"
    assert result[0].ref_seq == "C"
    assert result[0].alt_seq == "T"
    assert result[0].hgvs() == "KRAS:g.25398281C>T"


def test_to_gene(variant):
    result = variant.to_gene(gene_name=True)
    assert len(result) == 1
    assert result[0].reference == "KRAS"
    assert result[0].start == 25398281
    assert result[0].strand == "-"
    assert result[0].ref_seq == "G"
    assert result[0].alt_seq == "A"
    assert result[0].hgvs() == "KRAS:g.25398281G>A"


def test_to_protein(variant):
    result = variant.to_protein(gene_name=True)
    assert len(result) == 1
    assert result[0].reference == "KRAS"
    assert result[0].start == 13
    assert result[0].strand == "-"
    assert result[0].ref_seq == "G"  # GGC
    assert result[0].alt_seq == "D"  # GAC
    assert result[0].hgvs() == "KRAS:p.G13D"


def test_to_transcript(variant):
    result = variant.to_rna(gene_name=True)
    assert len(result) == 1
    assert result[0].reference == "KRAS"
    assert result[0].start == 102
    assert result[0].strand == "-"
    assert result[0].ref_seq == "G"
    assert result[0].alt_seq == "A"
    assert result[0].hgvs() == "KRAS:r.102G>A"


# class TestCdsSubstitutionWithOffsetCheckCorrectSequence:
# NM_130799:c.654+3A>G (#159:3164)
@pytest.fixture
def variant():
    return CdsSubstitution(
        reference="ENST00000312049",
        reference_type="transcript",
        start=654,
        ref_seq="A",
        alt_seq="G",
        start_offset=3,
    )


def test_hgvs(variant):
    assert variant.hgvs() == "ENST00000312049:c.654+3A>G"


def test_to_cds(variant):
    result = variant.to_contig()
    assert len(result) == 1
    assert result[0].reference == "11"
    assert result[0].start == 64575360
    assert result[0].strand == "+"
    assert result[0].ref_seq == "T"
    assert result[0].alt_seq == "C"


# class TestCdsSubstitutionOffset:
# example: TERT:c.-124C>T
@pytest.fixture
def variant():
    return CdsSubstitution(
        reference="TERT",
        reference_type="gene",
        start=1,
        ref_seq="C",
        alt_seq="T",
        start_offset=-124,
    )


def test_hgvs(variant):
    assert variant.hgvs() == "TERT:c.-124C>T"


def test_to_contig(variant):
    result = variant.to_contig()
    assert len(result) == 1
    assert result[0].reference == "5"
    assert result[0].start == 1295228
    assert result[0].end == 1295228
    assert result[0].ref_seq == "G"
    assert result[0].alt_seq == "A"


# class TestCdsSubstitutionSynonymous:
@pytest.fixture
def variant():
    return CdsSubstitution(
        reference="ENST00000372098",
        reference_type="transcript",
        start=1041,
        ref_seq="C",
        alt_seq="T",
    )


def test_to_protein(variant):
    result = variant.to_protein()
    assert len(result) == 1
    assert result[0].reference == "ENSP00000361170"
    assert result[0].start == 347
    assert result[0].strand == "-"
    assert result[0].ref_seq == "D"  # GAC == D
    assert result[0].alt_seq == "D"  # GAT == D
    assert result[0].hgvs() == "ENSP00000361170:p.D347="


# class TestDnaSubstitution:
@pytest.fixture
def variant():
    return DnaSubstitution(
        reference="12", reference_type="contig", start=25380283, ref_seq="C", alt_seq="A"
    )


def test_hgvs(variant):
    assert variant.hgvs() == "12:g.25380283C>A"


def test_to_cds(variant):
    result = variant.to_cds()
    assert len(result) == 2
    assert result[0].reference == "ENST00000256078"
    assert result[0].start == 175
    assert result[0].strand == "-"
    assert result[0].ref_seq == "G"
    assert result[0].alt_seq == "T"


def test_to_contig(variant):
    result = variant.to_contig()
    assert len(result) == 1
    assert result[0].reference == "12"
    assert result[0].start == 25380283
    assert result[0].strand == "+"
    assert result[0].ref_seq == "C"
    assert result[0].alt_seq == "A"


def test_to_gene(variant):
    result = variant.to_gene()
    assert len(result) == 1
    assert result[0].reference == "KRAS"
    assert result[0].start == 25380283
    assert result[0].strand == "-"
    assert result[0].ref_seq == "G"
    assert result[0].alt_seq == "T"


def test_to_protein(variant):
    result = variant.to_protein()
    assert len(result) == 2
    assert result[0].reference == "ENSP00000256078"
    assert result[0].start == 59
    assert result[0].strand == "-"
    assert result[0].ref_seq == "A"
    assert result[0].alt_seq == "S"


def test_to_transcript(variant):
    result = variant.to_rna()
    assert len(result) == 2
    assert result[0].reference == "ENST00000256078"
    assert result[0].start == 239
    assert result[0].strand == "-"
    assert result[0].ref_seq == "G"
    assert result[0].alt_seq == "T"


# class TestDnaSubstitutionGeneName:
@pytest.fixture
def variant():
    return DnaSubstitution(
        reference="12", reference_type="contig", start=25380283, ref_seq="C", alt_seq="A"
    )


def test_to_cds(variant):
    result = variant.to_cds(gene_name=True)
    assert len(result) == 2
    assert result[0].reference == "KRAS"
    assert result[0].start == 175
    assert result[0].strand == "-"
    assert result[0].ref_seq == "G"
    assert result[0].alt_seq == "T"


def test_to_contig(variant):
    result = variant.to_contig(gene_name=True)
    assert len(result) == 1
    assert result[0].reference == "KRAS"
    assert result[0].start == 25380283
    assert result[0].strand == "+"
    assert result[0].ref_seq == "C"
    assert result[0].alt_seq == "A"


def test_to_gene(variant):
    result = variant.to_gene(gene_name=True)
    assert len(result) == 1
    assert result[0].reference == "KRAS"
    assert result[0].start == 25380283
    assert result[0].strand == "-"
    assert result[0].ref_seq == "G"
    assert result[0].alt_seq == "T"


def test_to_protein(variant):
    result = variant.to_protein(gene_name=True)
    assert len(result) == 1
    assert result[0].reference == "KRAS"
    assert result[0].start == 59
    assert result[0].strand == "-"
    assert result[0].ref_seq == "A"
    assert result[0].alt_seq == "S"


def test_to_transcript(variant):
    result = variant.to_rna(gene_name=True)
    assert len(result) == 2
    assert result[0].reference == "KRAS"
    assert result[0].start == 239
    assert result[0].strand == "-"
    assert result[0].ref_seq == "G"
    assert result[0].alt_seq == "T"


# class TestProteinSubstitution:
@pytest.fixture
def variant():
    return ProteinSubstitution(
        reference="KRAS", reference_type="gene", start=13, ref_seq="G", alt_seq="V"
    )


def test_hgvs(variant):
    assert variant.hgvs() == "KRAS:p.G13V"


def test_to_cds(variant):
    result = variant.to_cds()
    assert len(result) == 16
    assert result[0].reference == "ENST00000256078"
    assert result[0].start == 38
    assert result[0].end == 39
    assert result[0].strand == "-"
    assert result[0].ref_seq == "GC"  # GGC
    assert result[0].alt_seq == "TA"  # GTA


def test_to_contig(variant):
    result = variant.to_contig()
    assert len(result) == 4
    assert result[0].reference == "12"
    assert result[0].start == 25398280
    assert result[0].end == 25398281
    assert result[0].strand == "+"
    assert result[0].ref_seq == "GC"  # GCC
    assert result[0].alt_seq == "TA"  # TAC


def test_to_gene(variant):
    result = variant.to_gene()
    assert len(result) == 4
    assert result[0].reference == "KRAS"
    assert result[0].start == 25398281
    assert result[0].end == 25398282
    assert result[0].strand == "-"
    assert result[0].ref_seq == "GC"
    assert result[0].alt_seq == "TA"


def test_to_protein(variant):
    result = variant.to_protein()
    assert len(result) == 4
    assert result[0].reference == "ENSP00000256078"
    assert result[0].start == 13
    assert result[0].strand == "-"
    assert result[0].ref_seq == "G"
    assert result[0].alt_seq == "V"


def test_to_transcript(variant):
    result = variant.to_rna()
    assert len(result) == 16
    assert result[0].reference == "ENST00000256078"
    assert result[0].start == 102
    assert result[0].end == 103
    assert result[0].strand == "-"
    assert result[0].ref_seq == "GC"
    assert result[0].alt_seq == "TA"


# class TestProteinSubstitutionGeneName:
@pytest.fixture
def variant():
    return ProteinSubstitution(
        reference="KRAS", reference_type="gene", start=13, ref_seq="G", alt_seq="V"
    )


def test_to_cds(variant):
    result = variant.to_cds(gene_name=True)
    assert len(result) == 16
    assert result[0].reference == "KRAS"
    assert result[0].start == 38
    assert result[0].end == 39
    assert result[0].strand == "-"
    assert result[0].ref_seq == "GC"  # GGC
    assert result[0].alt_seq == "TA"  # GTA


def test_to_contig(variant):
    result = variant.to_contig(gene_name=True)
    assert len(result) == 4
    assert result[0].reference == "KRAS"
    assert result[0].start == 25398280
    assert result[0].end == 25398281
    assert result[0].strand == "+"
    assert result[0].ref_seq == "GC"  # GCC
    assert result[0].alt_seq == "TA"  # TAC


def test_to_gene(variant):
    result = variant.to_gene(gene_name=True)
    assert len(result) == 4
    assert result[0].reference == "KRAS"
    assert result[0].start == 25398281
    assert result[0].end == 25398282
    assert result[0].strand == "-"
    assert result[0].ref_seq == "GC"
    assert result[0].alt_seq == "TA"


def test_to_protein(variant):
    result = variant.to_protein(gene_name=True)
    assert len(result) == 4
    assert result[0].reference == "KRAS"
    assert result[0].start == 13
    assert result[0].strand == "-"
    assert result[0].ref_seq == "G"
    assert result[0].alt_seq == "V"


def test_to_transcript(variant):
    result = variant.to_rna(gene_name=True)
    assert len(result) == 16
    assert result[0].reference == "KRAS"
    assert result[0].start == 102
    assert result[0].end == 103
    assert result[0].strand == "-"
    assert result[0].ref_seq == "GC"
    assert result[0].alt_seq == "TA"


# class TestProteinSubstitutionMissingAlt:
# example: KRAS:p.G13 #157:2942
def test_missing_alt_seq(variant):
    return ProteinSubstitution(reference="KRAS", reference_type="gene", start=13, ref_seq="G")
    assert variant.alt_seq == "X"


# class TestProteinSubstitutionRefMismatch:
# SDEV-2077 - ASXL2:p.R591X was incorrectly loaded as ASXL2:p.L591V
@pytest.fixture
def variant():
    return ProteinSubstitution(
        reference="ASXL2", reference_type="gene", start=591, ref_seq="R", alt_seq="X"
    )


def test_to_protein(variant):
    result = variant.to_protein()
    assert len(result) == 1
    assert result[0].reference == "ENSP00000391447"
    assert result[0].start == 591
    assert result[0].strand == "-"
    assert result[0].ref_seq == "R"
    assert result[0].alt_seq == "X"


# class TestRnaSubstitution:
@pytest.fixture
def variant():
    return RnaSubstitution(
        reference="ENST00000256078", reference_type="transcript", start=96, ref_seq="C", alt_seq="A"
    )


def test_hgvs(variant):
    assert variant.hgvs() == "ENST00000256078:r.96C>A"


def test_to_cds(variant):
    result = variant.to_cds()
    assert len(result) == 1
    assert result[0].reference == "ENST00000256078"
    assert result[0].start == 32
    assert result[0].strand == "-"
    assert result[0].ref_seq == "C"
    assert result[0].alt_seq == "A"


def test_to_contig(variant):
    result = variant.to_contig()
    assert len(result) == 1
    assert result[0].reference == "12"
    assert result[0].start == 25398287
    assert result[0].strand == "+"
    assert result[0].ref_seq == "G"
    assert result[0].alt_seq == "T"


def test_to_gene(variant):
    result = variant.to_gene()
    assert len(result) == 1
    assert result[0].reference == "KRAS"
    assert result[0].start == 25398287
    assert result[0].strand == "-"
    assert result[0].ref_seq == "C"
    assert result[0].alt_seq == "A"


def test_to_protein(variant):
    result = variant.to_protein()
    assert len(result) == 1
    assert result[0].reference == "ENSP00000256078"
    assert result[0].start == 11
    assert result[0].strand == "-"
    assert result[0].ref_seq == "A"
    assert result[0].alt_seq == "D"


def test_to_transcript(variant):
    result = variant.to_rna()
    assert len(result) == 1
    assert result[0].reference == "ENST00000256078"
    assert result[0].start == 96
    assert result[0].strand == "-"
    assert result[0].ref_seq == "C"
    assert result[0].alt_seq == "A"


# class TestRnaSubstitutionGeneName:
@pytest.fixture
def variant():
    return RnaSubstitution(
        reference="ENST00000256078", reference_type="transcript", start=96, ref_seq="C", alt_seq="A"
    )


def test_to_cds(variant):
    result = variant.to_cds(gene_name=True)
    assert len(result) == 1
    assert result[0].reference == "KRAS"
    assert result[0].start == 32
    assert result[0].strand == "-"
    assert result[0].ref_seq == "C"
    assert result[0].alt_seq == "A"


def test_to_contig(variant):
    result = variant.to_contig(gene_name=True)
    assert len(result) == 1
    assert result[0].reference == "KRAS"
    assert result[0].start == 25398287
    assert result[0].strand == "+"
    assert result[0].ref_seq == "G"
    assert result[0].alt_seq == "T"


def test_to_gene(variant):
    result = variant.to_gene(gene_name=True)
    assert len(result) == 1
    assert result[0].reference == "KRAS"
    assert result[0].start == 25398287
    assert result[0].strand == "-"
    assert result[0].ref_seq == "C"
    assert result[0].alt_seq == "A"


def test_to_protein(variant):
    result = variant.to_protein(gene_name=True)
    assert len(result) == 1
    assert result[0].reference == "KRAS"
    assert result[0].start == 11
    assert result[0].strand == "-"
    assert result[0].ref_seq == "A"
    assert result[0].alt_seq == "D"


def test_to_transcript(variant):
    result = variant.to_rna(gene_name=True)
    assert len(result) == 1
    assert result[0].reference == "KRAS"
    assert result[0].start == 96
    assert result[0].strand == "-"
    assert result[0].ref_seq == "C"
    assert result[0].alt_seq == "A"


def test_variant_repr_hgvs():
    variant = RnaSubstitution(
        reference="ENST00000256078", reference_type="transcript", start=96, ref_seq="C", alt_seq="A"
    )
    assert repr(variant) == "ENST00000256078:r.96C>A"


def test_variant_repr_error_user_repr_fallback():
    def raise_error():
        raise Exception

    variant = RnaSubstitution(
        reference="ENST00000256078",
        reference_type="transcript",
        start=96,
        ref_seq="C",
        alt_seq="A",
        user_repr="spam",
    )
    variant.hgvs = raise_error
    assert repr(variant) == "spam"


def test_variant_repr_error_classname_fallback():
    def raise_error():
        raise Exception

    variant = RnaSubstitution(
        reference="ENST00000256078",
        reference_type="transcript",
        start=96,
        ref_seq="C",
        alt_seq="A",
        user_repr=None,
    )
    variant.hgvs = raise_error
    assert repr(variant) == "RnaSubstitution"


def test_get_strand():
    variant = CdsSubstitution(
        reference="KIT", reference_type="gene", start=38, ref_seq="G", alt_seq="A", strand=None
    )
    assert variant.strand == "+"

    variant = CdsSubstitution(
        reference="TERT", reference_type="gene", start=38, ref_seq="G", alt_seq="A", strand=None
    )
    assert variant.strand == "-"

    variant = DnaSubstitution(
        reference="12", reference_type="contig", start=25380283, ref_seq="C", alt_seq="A"
    )
    assert variant.strand == "+"
