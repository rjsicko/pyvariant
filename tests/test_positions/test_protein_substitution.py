import pytest

from pyvariant.constants import SUBSTITUTION
from pyvariant.positions import (
    CdnaDelins,
    CdnaSubstitution,
    DnaDelins,
    DnaSubstitution,
    ProteinSubstitution,
    RnaDelins,
    RnaSubstitution,
)


@pytest.fixture()
def variant(ensembl69):
    return ProteinSubstitution(
        _core=ensembl69,
        contig_id="12",
        start=13,
        start_offset=0,
        end=13,
        end_offset=0,
        strand="-",
        gene_id="ENSG00000133703",
        gene_name="KRAS",
        transcript_id="ENST00000256078",
        transcript_name="KRAS-004",
        protein_id="ENSP00000256078",
        refseq="G",
        altseq="V",
    )


def test_str(variant):
    assert str(variant) == "ENSP00000256078:p.G13V"


def test_to_string_gene_name(variant):
    assert variant.to_string(reference="gene_name") == "KRAS:p.G13V"


def test_type(variant):
    assert not variant.is_cdna
    assert not variant.is_dna
    assert not variant.is_exon
    assert variant.is_protein
    assert not variant.is_rna


def test_variant_type(variant):
    assert variant.variant_type == SUBSTITUTION
    assert not variant.is_deletion
    assert not variant.is_delins
    assert not variant.is_duplication
    assert not variant.is_frameshift
    assert not variant.is_fusion
    assert not variant.is_insertion
    assert variant.is_substitution


def test_to_cdna(ensembl69, variant):
    expected_substitution = CdnaSubstitution(
        _core=ensembl69,
        contig_id="12",
        start=38,
        start_offset=0,
        end=38,
        end_offset=0,
        strand="-",
        gene_id="ENSG00000133703",
        gene_name="KRAS",
        transcript_id="ENST00000256078",
        transcript_name="KRAS-004",
        protein_id="ENSP00000256078",
        refseq="G",  # GGC
        altseq="T",  # GTC
    )
    expected_delins = CdnaDelins(
        _core=ensembl69,
        contig_id="12",
        start=38,
        start_offset=0,
        end=39,
        end_offset=0,
        strand="-",
        gene_id="ENSG00000133703",
        gene_name="KRAS",
        transcript_id="ENST00000256078",
        transcript_name="KRAS-004",
        protein_id="ENSP00000256078",
        refseq="GC",  # GGC
        altseq="TA",  # GTA
    )
    result = variant.to_cdna()
    assert expected_substitution in result
    assert expected_delins in result
    assert len(result) == 4


def test_to_dna(ensembl69, variant):
    expected_substitution = DnaSubstitution(
        _core=ensembl69,
        contig_id="12",
        start=25398281,
        start_offset=0,
        end=25398281,
        end_offset=0,
        strand="-",
        refseq="G",  # GGC
        altseq="T",  # GTC
    )
    expected_delins = DnaDelins(
        _core=ensembl69,
        contig_id="12",
        start=25398280,
        start_offset=0,
        end=25398281,
        end_offset=0,
        strand="-",
        refseq="GC",  # GGC
        altseq="TA",  # GTA
    )
    result = variant.to_dna()
    assert expected_substitution in result
    assert expected_delins in result
    assert len(result) == 4


def test_to_protein(ensembl69, variant):
    expected = [
        ProteinSubstitution(
            _core=ensembl69,
            contig_id="12",
            start=13,
            start_offset=0,
            end=13,
            end_offset=0,
            strand="-",
            gene_id="ENSG00000133703",
            gene_name="KRAS",
            transcript_id="ENST00000256078",
            transcript_name="KRAS-004",
            protein_id="ENSP00000256078",
            refseq="G",
            altseq="V",
        )
    ]
    assert variant.to_protein() == expected


def test_to_rna(ensembl69, variant):
    expected_substitution = RnaSubstitution(
        _core=ensembl69,
        contig_id="12",
        start=102,
        start_offset=0,
        end=102,
        end_offset=0,
        strand="-",
        gene_id="ENSG00000133703",
        gene_name="KRAS",
        transcript_id="ENST00000256078",
        transcript_name="KRAS-004",
        refseq="G",  # GGC
        altseq="T",  # GTC
    )
    expected_delins = RnaDelins(
        _core=ensembl69,
        contig_id="12",
        start=102,
        start_offset=0,
        end=103,
        end_offset=0,
        strand="-",
        gene_id="ENSG00000133703",
        gene_name="KRAS",
        transcript_id="ENST00000256078",
        transcript_name="KRAS-004",
        refseq="GC",  # GGC
        altseq="TA",  # GTA
    )
    result = variant.to_rna()
    assert expected_substitution in result
    assert expected_delins in result
    assert len(result) == 4
