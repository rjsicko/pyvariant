import pytest

from variant_map.positions import (
    CdnaPosition,
    DnaPosition,
    ExonPosition,
    ProteinPosition,
    RnaPosition,
)


def test_is_cdna():
    from_position = RnaPosition(
        contig_id="5",
        start=100,
        start_offset=-1,
        end=100,
        end_offset=-1,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
    )
    assert from_position.is_cdna is False


def test_is_dna():
    from_position = RnaPosition(
        contig_id="5",
        start=100,
        start_offset=-1,
        end=100,
        end_offset=-1,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
    )
    assert from_position.is_dna is False


def test_is_exon():
    from_position = RnaPosition(
        contig_id="5",
        start=100,
        start_offset=-1,
        end=100,
        end_offset=-1,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
    )
    assert from_position.is_exon is False


def test_is_protein():
    from_position = RnaPosition(
        contig_id="5",
        start=100,
        start_offset=-1,
        end=100,
        end_offset=-1,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
    )
    assert from_position.is_protein is False


def test_is_rna():
    from_position = RnaPosition(
        contig_id="5",
        start=100,
        start_offset=-1,
        end=100,
        end_offset=-1,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
    )
    assert from_position.is_rna is True


def test_is_on_negative_strand():
    from_position = RnaPosition(
        contig_id="5",
        start=100,
        start_offset=-1,
        end=100,
        end_offset=-1,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
    )
    assert from_position.on_negative_strand is True


def test_is_on_positive_strand():
    from_position = RnaPosition(
        contig_id="5",
        start=100,
        start_offset=-1,
        end=100,
        end_offset=-1,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
    )
    assert from_position.on_positive_strand is False


def test_sequence(ensembl100):
    from_position = RnaPosition(
        contig_id="5",
        start=1652,
        start_offset=0,
        end=1654,
        end_offset=0,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
    )
    assert ensembl100.sequence(from_position) == "GGG"


@pytest.mark.skip(reason="how to handle?")
def test_sequence_offset(ensembl100):
    # Position crosses exon boundary so getting DNA sequence returns intron sequences as well
    # TODO: How to handle sequences for offset variants?
    from_position = RnaPosition(
        contig_id="5",
        start=1651,
        start_offset=1,
        end=1653,
        end_offset=1,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
    )
    assert ensembl100.sequence(from_position) == "GGG"


def test_str_mutli_position():
    from_position = RnaPosition(
        contig_id="5",
        start=100,
        start_offset=-1,
        end=101,
        end_offset=-1,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
    )
    assert str(from_position) == "ENST00000310581:r.100-1_101-1"


def test_str_one_position():
    from_position = RnaPosition(
        contig_id="5",
        start=100,
        start_offset=-1,
        end=100,
        end_offset=-1,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
    )
    assert str(from_position) == "ENST00000310581:r.100-1"


def test_to_cdna_rna_end_minus_strand(ensembl100):
    from_position = RnaPosition(
        contig_id="5",
        start=4037,
        start_offset=0,
        end=4039,
        end_offset=0,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
    )
    to_positions = []
    assert ensembl100.to_cdna(from_position) == to_positions


def test_to_cdna_rna_end_plus_strand(ensembl100):
    from_position = RnaPosition(
        contig_id="13",
        start=11952,
        start_offset=0,
        end=11954,
        end_offset=0,
        strand="+",
        gene_id="ENSG00000139618",
        gene_name="BRCA2",
        transcript_id="ENST00000380152",
        transcript_name="BRCA2-201",
    )
    to_positions = []
    assert ensembl100.to_cdna(from_position) == to_positions


def test_to_cdna_rna_exon_boundary_minus_strand(ensembl100):
    from_position = RnaPosition(
        contig_id="5",
        start=1652,
        start_offset=0,
        end=1654,
        end_offset=0,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
    )
    to_positions = [
        CdnaPosition(
            contig_id="5",
            start=1573,
            start_offset=0,
            end=1575,
            end_offset=0,
            strand="-",
            gene_id="ENSG00000164362",
            gene_name="TERT",
            transcript_id="ENST00000310581",
            transcript_name="TERT-201",
            protein_id="ENSP00000309572",
        )
    ]
    assert ensembl100.to_cdna(from_position) == to_positions


def test_to_cdna_rna_exon_boundary_plus_strand(ensembl100):
    from_position = RnaPosition(
        contig_id="13",
        start=300,
        start_offset=0,
        end=302,
        end_offset=0,
        strand="+",
        gene_id="ENSG00000139618",
        gene_name="BRCA2",
        transcript_id="ENST00000380152",
        transcript_name="BRCA2-201",
    )
    to_positions = [
        CdnaPosition(
            contig_id="13",
            start=67,
            start_offset=0,
            end=69,
            end_offset=0,
            strand="+",
            gene_id="ENSG00000139618",
            gene_name="BRCA2",
            transcript_id="ENST00000380152",
            transcript_name="BRCA2-201",
            protein_id="ENSP00000369497",
        )
    ]
    assert ensembl100.to_cdna(from_position) == to_positions


def test_to_cdna_rna_negative_offset_normalized_minus_strand(ensembl100):
    from_position = RnaPosition(
        contig_id="5",
        start=100,
        start_offset=-1,
        end=100,
        end_offset=-1,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
    )
    to_positions = [
        CdnaPosition(
            contig_id="5",
            start=20,
            start_offset=0,
            end=20,
            end_offset=0,
            strand="-",
            gene_id="ENSG00000164362",
            gene_name="TERT",
            transcript_id="ENST00000310581",
            transcript_name="TERT-201",
            protein_id="ENSP00000309572",
        )
    ]
    assert ensembl100.to_cdna(from_position) == to_positions


def test_to_cdna_rna_positive_offset_minus_strand(ensembl100):
    from_position = RnaPosition(
        contig_id="5",
        start=1652,
        start_offset=0,
        end=1654,
        end_offset=0,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
    )
    to_positions = [
        CdnaPosition(
            contig_id="5",
            start=1573,
            start_offset=0,
            end=1575,
            end_offset=0,
            strand="-",
            gene_id="ENSG00000164362",
            gene_name="TERT",
            transcript_id="ENST00000310581",
            transcript_name="TERT-201",
            protein_id="ENSP00000309572",
        )
    ]
    assert ensembl100.to_cdna(from_position) == to_positions


def test_to_cdna_rna_start_minus_strand(ensembl100):
    from_position = RnaPosition(
        contig_id="5",
        start=1,
        start_offset=0,
        end=3,
        end_offset=0,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
    )
    to_positions = []
    assert ensembl100.to_cdna(from_position) == to_positions


def test_to_cdna_rna_start_plus_strand(ensembl100):
    from_position = RnaPosition(
        contig_id="13",
        start=1,
        start_offset=0,
        end=3,
        end_offset=0,
        strand="+",
        gene_id="ENSG00000139618",
        gene_name="BRCA2",
        transcript_id="ENST00000380152",
        transcript_name="BRCA2-201",
    )
    to_positions = []
    assert ensembl100.to_cdna(from_position) == to_positions


def test_to_dna_rna_end_minus_strand(ensembl100):
    from_position = RnaPosition(
        contig_id="5",
        start=4037,
        start_offset=0,
        end=4039,
        end_offset=0,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
    )
    to_positions = [
        DnaPosition(
            contig_id="5", start=1253167, start_offset=0, end=1253169, end_offset=0, strand="-"
        )
    ]
    assert ensembl100.to_dna(from_position) == to_positions


def test_to_dna_rna_end_plus_strand(ensembl100):
    from_position = RnaPosition(
        contig_id="13",
        start=11952,
        start_offset=0,
        end=11954,
        end_offset=0,
        strand="+",
        gene_id="ENSG00000139618",
        gene_name="BRCA2",
        transcript_id="ENST00000380152",
        transcript_name="BRCA2-201",
    )
    to_positions = [
        DnaPosition(
            contig_id="13", start=32400232, start_offset=0, end=32400234, end_offset=0, strand="+"
        )
    ]
    assert ensembl100.to_dna(from_position) == to_positions


def test_to_dna_rna_exon_boundary_minus_strand(ensembl100):
    from_position = RnaPosition(
        contig_id="5",
        start=1652,
        start_offset=0,
        end=1654,
        end_offset=0,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
    )
    to_positions = [
        DnaPosition(
            contig_id="5", start=1282623, start_offset=0, end=1293313, end_offset=0, strand="-"
        )
    ]
    assert ensembl100.to_dna(from_position) == to_positions


def test_to_dna_rna_exon_boundary_plus_strand(ensembl100):
    from_position = RnaPosition(
        contig_id="13",
        start=300,
        start_offset=0,
        end=302,
        end_offset=0,
        strand="+",
        gene_id="ENSG00000139618",
        gene_name="BRCA2",
        transcript_id="ENST00000380152",
        transcript_name="BRCA2-201",
    )
    to_positions = [
        DnaPosition(
            contig_id="13", start=32316527, start_offset=0, end=32319078, end_offset=0, strand="+"
        )
    ]
    assert ensembl100.to_dna(from_position) == to_positions


def test_to_dna_rna_negative_offset_normalized_minus_strand(ensembl100):
    from_position = RnaPosition(
        contig_id="5",
        start=100,
        start_offset=-1,
        end=100,
        end_offset=-1,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
    )
    to_positions = [
        DnaPosition(
            contig_id="5", start=1294970, start_offset=0, end=1294970, end_offset=0, strand="-"
        )
    ]
    assert ensembl100.to_dna(from_position) == to_positions


def test_to_dna_rna_positive_offset_minus_strand(ensembl100):
    from_position = RnaPosition(
        contig_id="5",
        start=1652,
        start_offset=0,
        end=1654,
        end_offset=0,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
    )
    to_positions = [
        DnaPosition(
            contig_id="5", start=1282623, start_offset=0, end=1293313, end_offset=0, strand="-"
        )
    ]
    assert ensembl100.to_dna(from_position) == to_positions


def test_to_dna_rna_start_minus_strand(ensembl100):
    from_position = RnaPosition(
        contig_id="5",
        start=1,
        start_offset=0,
        end=3,
        end_offset=0,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
    )
    to_positions = [
        DnaPosition(
            contig_id="5", start=1295066, start_offset=0, end=1295068, end_offset=0, strand="-"
        )
    ]
    assert ensembl100.to_dna(from_position) == to_positions


def test_to_dna_rna_start_plus_strand(ensembl100):
    from_position = RnaPosition(
        contig_id="13",
        start=1,
        start_offset=0,
        end=3,
        end_offset=0,
        strand="+",
        gene_id="ENSG00000139618",
        gene_name="BRCA2",
        transcript_id="ENST00000380152",
        transcript_name="BRCA2-201",
    )
    to_positions = [
        DnaPosition(
            contig_id="13", start=32315474, start_offset=0, end=32315476, end_offset=0, strand="+"
        )
    ]
    assert ensembl100.to_dna(from_position) == to_positions


def test_to_exon_rna_end_minus_strand(ensembl100):
    from_position = RnaPosition(
        contig_id="5",
        start=4037,
        start_offset=0,
        end=4039,
        end_offset=0,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
    )
    to_positions = [
        ExonPosition(
            contig_id="5",
            start=16,
            start_offset=0,
            end=16,
            end_offset=0,
            strand="-",
            gene_id="ENSG00000164362",
            gene_name="TERT",
            transcript_id="ENST00000310581",
            transcript_name="TERT-201",
            exon_id="ENSE00001863787",
        )
    ]
    assert ensembl100.to_exon(from_position) == to_positions


def test_to_exon_rna_end_plus_strand(ensembl100):
    from_position = RnaPosition(
        contig_id="13",
        start=11952,
        start_offset=0,
        end=11954,
        end_offset=0,
        strand="+",
        gene_id="ENSG00000139618",
        gene_name="BRCA2",
        transcript_id="ENST00000380152",
        transcript_name="BRCA2-201",
    )
    to_positions = [
        ExonPosition(
            contig_id="13",
            start=27,
            start_offset=0,
            end=27,
            end_offset=0,
            strand="+",
            gene_id="ENSG00000139618",
            gene_name="BRCA2",
            transcript_id="ENST00000380152",
            transcript_name="BRCA2-201",
            exon_id="ENSE00003717596",
        )
    ]
    assert ensembl100.to_exon(from_position) == to_positions


def test_to_exon_rna_exon_boundary_minus_strand(ensembl100):
    from_position = RnaPosition(
        contig_id="5",
        start=1652,
        start_offset=0,
        end=1654,
        end_offset=0,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
    )
    to_positions = [
        ExonPosition(
            contig_id="5",
            start=2,
            start_offset=0,
            end=3,
            end_offset=0,
            strand="-",
            gene_id="ENSG00000164362",
            gene_name="TERT",
            transcript_id="ENST00000310581",
            transcript_name="TERT-201",
            exon_id="ENSE00001197112",
        )
    ]
    assert ensembl100.to_exon(from_position) == to_positions


def test_to_exon_rna_exon_boundary_plus_strand(ensembl100):
    from_position = RnaPosition(
        contig_id="13",
        start=300,
        start_offset=0,
        end=302,
        end_offset=0,
        strand="+",
        gene_id="ENSG00000139618",
        gene_name="BRCA2",
        transcript_id="ENST00000380152",
        transcript_name="BRCA2-201",
    )
    to_positions = [
        ExonPosition(
            contig_id="13",
            start=2,
            start_offset=0,
            end=3,
            end_offset=0,
            strand="+",
            gene_id="ENSG00000139618",
            gene_name="BRCA2",
            transcript_id="ENST00000380152",
            transcript_name="BRCA2-201",
            exon_id="ENSE00001484009",
        )
    ]
    assert ensembl100.to_exon(from_position) == to_positions


def test_to_exon_rna_negative_offset_normalized_minus_strand(ensembl100):
    from_position = RnaPosition(
        contig_id="5",
        start=100,
        start_offset=-1,
        end=100,
        end_offset=-1,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
    )
    to_positions = [
        ExonPosition(
            contig_id="5",
            start=1,
            start_offset=0,
            end=1,
            end_offset=0,
            strand="-",
            gene_id="ENSG00000164362",
            gene_name="TERT",
            transcript_id="ENST00000310581",
            transcript_name="TERT-201",
            exon_id="ENSE00003896691",
        )
    ]
    assert ensembl100.to_exon(from_position) == to_positions


def test_to_exon_rna_positive_offset_minus_strand(ensembl100):
    from_position = RnaPosition(
        contig_id="5",
        start=1652,
        start_offset=0,
        end=1654,
        end_offset=0,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
    )
    to_positions = [
        ExonPosition(
            contig_id="5",
            start=2,
            start_offset=0,
            end=3,
            end_offset=0,
            strand="-",
            gene_id="ENSG00000164362",
            gene_name="TERT",
            transcript_id="ENST00000310581",
            transcript_name="TERT-201",
            exon_id="ENSE00001197112",
        )
    ]
    assert ensembl100.to_exon(from_position) == to_positions


def test_to_exon_rna_start_minus_strand(ensembl100):
    from_position = RnaPosition(
        contig_id="5",
        start=1,
        start_offset=0,
        end=3,
        end_offset=0,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
    )
    to_positions = [
        ExonPosition(
            contig_id="5",
            start=1,
            start_offset=0,
            end=1,
            end_offset=0,
            strand="-",
            gene_id="ENSG00000164362",
            gene_name="TERT",
            transcript_id="ENST00000310581",
            transcript_name="TERT-201",
            exon_id="ENSE00003896691",
        )
    ]
    assert ensembl100.to_exon(from_position) == to_positions


def test_to_exon_rna_start_plus_strand(ensembl100):
    from_position = RnaPosition(
        contig_id="13",
        start=1,
        start_offset=0,
        end=3,
        end_offset=0,
        strand="+",
        gene_id="ENSG00000139618",
        gene_name="BRCA2",
        transcript_id="ENST00000380152",
        transcript_name="BRCA2-201",
    )
    to_positions = [
        ExonPosition(
            contig_id="13",
            start=1,
            start_offset=0,
            end=1,
            end_offset=0,
            strand="+",
            gene_id="ENSG00000139618",
            gene_name="BRCA2",
            transcript_id="ENST00000380152",
            transcript_name="BRCA2-201",
            exon_id="ENSE00001184784",
        )
    ]
    assert ensembl100.to_exon(from_position) == to_positions


def test_to_protein_rna_end_minus_strand(ensembl100):
    from_position = RnaPosition(
        contig_id="5",
        start=4037,
        start_offset=0,
        end=4039,
        end_offset=0,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
    )
    to_positions = []
    assert ensembl100.to_protein(from_position) == to_positions


def test_to_protein_rna_end_plus_strand(ensembl100):
    from_position = RnaPosition(
        contig_id="13",
        start=11952,
        start_offset=0,
        end=11954,
        end_offset=0,
        strand="+",
        gene_id="ENSG00000139618",
        gene_name="BRCA2",
        transcript_id="ENST00000380152",
        transcript_name="BRCA2-201",
    )
    to_positions = []
    assert ensembl100.to_protein(from_position) == to_positions


def test_to_protein_rna_exon_boundary_minus_strand(ensembl100):
    from_position = RnaPosition(
        contig_id="5",
        start=1652,
        start_offset=0,
        end=1654,
        end_offset=0,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
    )
    to_positions = [
        ProteinPosition(
            contig_id="5",
            start=525,
            start_offset=0,
            end=525,
            end_offset=0,
            strand="-",
            gene_id="ENSG00000164362",
            gene_name="TERT",
            transcript_id="ENST00000310581",
            transcript_name="TERT-201",
            protein_id="ENSP00000309572",
        )
    ]
    assert ensembl100.to_protein(from_position) == to_positions


def test_to_protein_rna_exon_boundary_plus_strand(ensembl100):
    from_position = RnaPosition(
        contig_id="13",
        start=300,
        start_offset=0,
        end=302,
        end_offset=0,
        strand="+",
        gene_id="ENSG00000139618",
        gene_name="BRCA2",
        transcript_id="ENST00000380152",
        transcript_name="BRCA2-201",
    )
    to_positions = [
        ProteinPosition(
            contig_id="13",
            start=23,
            start_offset=0,
            end=23,
            end_offset=0,
            strand="+",
            gene_id="ENSG00000139618",
            gene_name="BRCA2",
            transcript_id="ENST00000380152",
            transcript_name="BRCA2-201",
            protein_id="ENSP00000369497",
        )
    ]
    assert ensembl100.to_protein(from_position) == to_positions


def test_to_protein_rna_negative_offset_normalized_minus_strand(ensembl100):
    from_position = RnaPosition(
        contig_id="5",
        start=100,
        start_offset=-1,
        end=100,
        end_offset=-1,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
    )
    to_positions = [
        ProteinPosition(
            contig_id="5",
            start=7,
            start_offset=0,
            end=7,
            end_offset=0,
            strand="-",
            gene_id="ENSG00000164362",
            gene_name="TERT",
            transcript_id="ENST00000310581",
            transcript_name="TERT-201",
            protein_id="ENSP00000309572",
        )
    ]
    assert ensembl100.to_protein(from_position) == to_positions


def test_to_protein_rna_positive_offset_minus_strand(ensembl100):
    from_position = RnaPosition(
        contig_id="5",
        start=1652,
        start_offset=0,
        end=1654,
        end_offset=0,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
    )
    to_positions = [
        ProteinPosition(
            contig_id="5",
            start=525,
            start_offset=0,
            end=525,
            end_offset=0,
            strand="-",
            gene_id="ENSG00000164362",
            gene_name="TERT",
            transcript_id="ENST00000310581",
            transcript_name="TERT-201",
            protein_id="ENSP00000309572",
        )
    ]
    assert ensembl100.to_protein(from_position) == to_positions


def test_to_protein_rna_start_minus_strand(ensembl100):
    from_position = RnaPosition(
        contig_id="5",
        start=1,
        start_offset=0,
        end=3,
        end_offset=0,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
    )
    to_positions = []
    assert ensembl100.to_protein(from_position) == to_positions


def test_to_protein_rna_start_plus_strand(ensembl100):
    from_position = RnaPosition(
        contig_id="13",
        start=1,
        start_offset=0,
        end=3,
        end_offset=0,
        strand="+",
        gene_id="ENSG00000139618",
        gene_name="BRCA2",
        transcript_id="ENST00000380152",
        transcript_name="BRCA2-201",
    )
    to_positions = []
    assert ensembl100.to_protein(from_position) == to_positions


def test_to_rna_rna_end_minus_strand(ensembl100):
    from_position = RnaPosition(
        contig_id="5",
        start=4037,
        start_offset=0,
        end=4039,
        end_offset=0,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
    )
    to_positions = [
        RnaPosition(
            contig_id="5",
            start=4037,
            start_offset=0,
            end=4039,
            end_offset=0,
            strand="-",
            gene_id="ENSG00000164362",
            gene_name="TERT",
            transcript_id="ENST00000310581",
            transcript_name="TERT-201",
        )
    ]
    assert ensembl100.to_rna(from_position) == to_positions


def test_to_rna_rna_end_plus_strand(ensembl100):
    from_position = RnaPosition(
        contig_id="13",
        start=11952,
        start_offset=0,
        end=11954,
        end_offset=0,
        strand="+",
        gene_id="ENSG00000139618",
        gene_name="BRCA2",
        transcript_id="ENST00000380152",
        transcript_name="BRCA2-201",
    )
    to_positions = [
        RnaPosition(
            contig_id="13",
            start=11952,
            start_offset=0,
            end=11954,
            end_offset=0,
            strand="+",
            gene_id="ENSG00000139618",
            gene_name="BRCA2",
            transcript_id="ENST00000380152",
            transcript_name="BRCA2-201",
        )
    ]
    assert ensembl100.to_rna(from_position) == to_positions


def test_to_rna_rna_exon_boundary_minus_strand(ensembl100):
    from_position = RnaPosition(
        contig_id="5",
        start=1652,
        start_offset=0,
        end=1654,
        end_offset=0,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
    )
    to_positions = [
        RnaPosition(
            contig_id="5",
            start=1652,
            start_offset=0,
            end=1654,
            end_offset=0,
            strand="-",
            gene_id="ENSG00000164362",
            gene_name="TERT",
            transcript_id="ENST00000310581",
            transcript_name="TERT-201",
        )
    ]
    assert ensembl100.to_rna(from_position) == to_positions


def test_to_rna_rna_exon_boundary_plus_strand(ensembl100):
    from_position = RnaPosition(
        contig_id="13",
        start=300,
        start_offset=0,
        end=302,
        end_offset=0,
        strand="+",
        gene_id="ENSG00000139618",
        gene_name="BRCA2",
        transcript_id="ENST00000380152",
        transcript_name="BRCA2-201",
    )
    to_positions = [
        RnaPosition(
            contig_id="13",
            start=300,
            start_offset=0,
            end=302,
            end_offset=0,
            strand="+",
            gene_id="ENSG00000139618",
            gene_name="BRCA2",
            transcript_id="ENST00000380152",
            transcript_name="BRCA2-201",
        )
    ]
    assert ensembl100.to_rna(from_position) == to_positions


def test_to_rna_rna_negative_offset_normalized_minus_strand(ensembl100):
    from_position = RnaPosition(
        contig_id="5",
        start=100,
        start_offset=-1,
        end=100,
        end_offset=-1,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
    )
    to_positions = [
        RnaPosition(
            contig_id="5",
            start=99,
            start_offset=0,
            end=99,
            end_offset=0,
            strand="-",
            gene_id="ENSG00000164362",
            gene_name="TERT",
            transcript_id="ENST00000310581",
            transcript_name="TERT-201",
        )
    ]
    assert ensembl100.to_rna(from_position) == to_positions


def test_to_rna_rna_positive_offset_minus_strand(ensembl100):
    from_position = RnaPosition(
        contig_id="5",
        start=1652,
        start_offset=0,
        end=1654,
        end_offset=0,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
    )
    to_positions = [
        RnaPosition(
            contig_id="5",
            start=1652,
            start_offset=0,
            end=1654,
            end_offset=0,
            strand="-",
            gene_id="ENSG00000164362",
            gene_name="TERT",
            transcript_id="ENST00000310581",
            transcript_name="TERT-201",
        )
    ]
    assert ensembl100.to_rna(from_position) == to_positions


def test_to_rna_rna_start_minus_strand(ensembl100):
    from_position = RnaPosition(
        contig_id="5",
        start=1,
        start_offset=0,
        end=3,
        end_offset=0,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
    )
    to_positions = [
        RnaPosition(
            contig_id="5",
            start=1,
            start_offset=0,
            end=3,
            end_offset=0,
            strand="-",
            gene_id="ENSG00000164362",
            gene_name="TERT",
            transcript_id="ENST00000310581",
            transcript_name="TERT-201",
        )
    ]
    assert ensembl100.to_rna(from_position) == to_positions


def test_to_rna_rna_start_plus_strand(ensembl100):
    from_position = RnaPosition(
        contig_id="13",
        start=1,
        start_offset=0,
        end=3,
        end_offset=0,
        strand="+",
        gene_id="ENSG00000139618",
        gene_name="BRCA2",
        transcript_id="ENST00000380152",
        transcript_name="BRCA2-201",
    )
    to_positions = [
        RnaPosition(
            contig_id="13",
            start=1,
            start_offset=0,
            end=3,
            end_offset=0,
            strand="+",
            gene_id="ENSG00000139618",
            gene_name="BRCA2",
            transcript_id="ENST00000380152",
            transcript_name="BRCA2-201",
        )
    ]
    assert ensembl100.to_rna(from_position) == to_positions
