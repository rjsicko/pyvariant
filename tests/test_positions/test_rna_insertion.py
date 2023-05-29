import pytest

from pyvariant.constants import INSERTION
from pyvariant.positions import CdnaInsertion, DnaInsertion, ProteinInsertion, RnaInsertion


@pytest.fixture()
def variant(ensembl69):
    return RnaInsertion(
        _core=ensembl69,
        contig_id="4",
        start=1771,
        start_offset=0,
        end=1772,
        end_offset=0,
        strand="+",
        gene_id="ENSG00000157404",
        gene_name="KIT",
        transcript_id="ENST00000288135",
        transcript_name="KIT-001",
        refseq="GG",
        altseq="GCATG",
    )


def test_str(variant):
    assert str(variant) == "ENST00000288135:r.1771_1772insCAT"


def test_to_string_gene_name(variant):
    assert variant.to_string(reference="gene_name") == "KIT:r.1771_1772insCAT"


def test_type(variant):
    assert not variant.is_cdna
    assert not variant.is_dna
    assert not variant.is_exon
    assert not variant.is_protein
    assert variant.is_rna


def test_variant_type(variant):
    assert variant.variant_type == INSERTION
    assert not variant.is_deletion
    assert not variant.is_delins
    assert not variant.is_duplication
    assert not variant.is_frameshift
    assert not variant.is_fusion
    assert variant.is_insertion
    assert not variant.is_substitution


def test_to_cdna(ensembl69, variant):
    expected = [
        CdnaInsertion(
            _core=ensembl69,
            contig_id="4",
            start=1674,
            start_offset=0,
            end=1675,
            end_offset=0,
            strand="+",
            gene_id="ENSG00000157404",
            gene_name="KIT",
            transcript_id="ENST00000288135",
            transcript_name="KIT-001",
            protein_id="ENSP00000288135",
            refseq="GG",
            altseq="GCATG",
        )
    ]
    assert variant.to_cdna() == expected


def test_to_dna(ensembl69, variant):
    expected = [
        DnaInsertion(
            _core=ensembl69,
            contig_id="4",
            start=55593608,
            start_offset=0,
            end=55593609,
            end_offset=0,
            strand="+",
            refseq="GG",
            altseq="GCATG",
        )
    ]
    assert variant.to_dna() == expected


def test_to_protein(ensembl69, variant):
    expected = [
        ProteinInsertion(
            _core=ensembl69,
            contig_id="4",
            start=558,
            start_offset=0,
            end=559,
            end_offset=0,
            strand="+",
            gene_id="ENSG00000157404",
            gene_name="KIT",
            transcript_id="ENST00000288135",
            transcript_name="KIT-001",
            protein_id="ENSP00000288135",
            refseq="KV",
            altseq="KHV",
        )
    ]
    assert variant.to_protein() == expected


def test_to_rna(ensembl69, variant):
    expected = [
        RnaInsertion(
            _core=ensembl69,
            contig_id="4",
            start=1771,
            start_offset=0,
            end=1772,
            end_offset=0,
            strand="+",
            gene_id="ENSG00000157404",
            gene_name="KIT",
            transcript_id="ENST00000288135",
            transcript_name="KIT-001",
            refseq="GG",
            altseq="GCATG",
        )
    ]
    assert variant.to_rna() == expected
