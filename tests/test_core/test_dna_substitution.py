import pytest

from variant_map.constants import SUBTITUTION
from variant_map.core import CdnaSubstitution, DnaSubstitution, ProteinSubstitution, RnaSubstitution


@pytest.fixture()
def variant(ensembl69):
    return DnaSubstitution(
        _data=ensembl69,
        contig_id="12",
        start=25380283,
        start_offset=0,
        end=25380283,
        end_offset=0,
        strand="-",
        refseq="G",
        altseq="T",
    )


def test_str(variant):
    assert str(variant) == "12:g.25380283G>T"


def test_type(variant):
    assert not variant.is_cdna
    assert variant.is_dna
    assert not variant.is_exon
    assert not variant.is_protein
    assert not variant.is_rna


def test_variant_type(variant):
    assert variant.type == SUBTITUTION
    assert not variant.is_deletion
    assert not variant.is_delins
    assert not variant.is_duplication
    assert not variant.is_frameshift
    assert not variant.is_fusion
    assert not variant.is_insertion
    assert variant.is_substitution


def test_to_cdna(ensembl69, variant):
    expected = [
        CdnaSubstitution(
            _data=ensembl69,
            contig_id="12",
            start=175,
            start_offset=0,
            end=175,
            end_offset=0,
            strand="-",
            gene_id="ENSG00000133703",
            gene_name="KRAS",
            transcript_id="ENST00000256078",
            transcript_name="KRAS-004",
            protein_id="ENSP00000256078",
            refseq="G",
            altseq="T",
        ),
        CdnaSubstitution(
            _data=ensembl69,
            contig_id="12",
            start=175,
            start_offset=0,
            end=175,
            end_offset=0,
            strand="-",
            gene_id="ENSG00000133703",
            gene_name="KRAS",
            transcript_id="ENST00000311936",
            transcript_name="KRAS-001",
            protein_id="ENSP00000308495",
            refseq="G",
            altseq="T",
        ),
    ]
    assert variant.to_cdna() == expected


def test_to_dna(ensembl69, variant):
    expected = [
        DnaSubstitution(
            _data=ensembl69,
            contig_id="12",
            start=25380283,
            start_offset=0,
            end=25380283,
            end_offset=0,
            strand="-",
            refseq="G",
            altseq="T",
        )
    ]
    assert variant.to_dna() == expected


def test_to_protein(ensembl69, variant):
    expected = [
        ProteinSubstitution(
            _data=ensembl69,
            contig_id="12",
            start=59,
            start_offset=0,
            end=59,
            end_offset=0,
            strand="-",
            gene_id="ENSG00000133703",
            gene_name="KRAS",
            transcript_id="ENST00000256078",
            transcript_name="KRAS-004",
            protein_id="ENSP00000256078",
            refseq="A",
            altseq="S",
        ),
        ProteinSubstitution(
            _data=ensembl69,
            contig_id="12",
            start=59,
            start_offset=0,
            end=59,
            end_offset=0,
            strand="-",
            gene_id="ENSG00000133703",
            gene_name="KRAS",
            transcript_id="ENST00000311936",
            transcript_name="KRAS-001",
            protein_id="ENSP00000308495",
            refseq="A",
            altseq="S",
        ),
    ]
    assert variant.to_protein() == expected


def test_to_rna(ensembl69, variant):
    expected = [
        RnaSubstitution(
            _data=ensembl69,
            contig_id="12",
            start=239,
            start_offset=0,
            end=239,
            end_offset=0,
            strand="-",
            gene_id="ENSG00000133703",
            gene_name="KRAS",
            transcript_id="ENST00000256078",
            transcript_name="KRAS-004",
            refseq="G",
            altseq="T",
        ),
        RnaSubstitution(
            _data=ensembl69,
            contig_id="12",
            start=367,
            start_offset=0,
            end=367,
            end_offset=0,
            strand="-",
            gene_id="ENSG00000133703",
            gene_name="KRAS",
            transcript_id="ENST00000311936",
            transcript_name="KRAS-001",
            refseq="G",
            altseq="T",
        ),
    ]
    assert variant.to_rna() == expected
