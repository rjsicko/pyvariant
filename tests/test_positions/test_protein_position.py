import pytest

from variant_map.positions import (
    CdnaPosition,
    DnaPosition,
    ExonPosition,
    ProteinPosition,
    RnaPosition,
)


def test_is_cdna():
    from_position = ProteinPosition(
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
    assert from_position.is_cdna is False


def test_is_dna():
    from_position = ProteinPosition(
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
    assert from_position.is_dna is False


def test_is_exon():
    from_position = ProteinPosition(
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
    assert from_position.is_exon is False


def test_is_protein():
    from_position = ProteinPosition(
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
    assert from_position.is_protein is True


def test_is_rna():
    from_position = ProteinPosition(
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
    assert from_position.is_rna is False


def test_is_on_negative_strand():
    from_position = ProteinPosition(
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
    assert from_position.on_negative_strand is True


def test_is_on_positive_strand():
    from_position = ProteinPosition(
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
    assert from_position.on_positive_strand is False


def test_sequence(ensembl100):
    from_position = ProteinPosition(
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
    assert ensembl100.sequence(from_position) == "G"


def test_sequence_offset(ensembl100):
    from_position = ProteinPosition(
        contig_id="5",
        start=525,
        start_offset=1,
        end=525,
        end_offset=1,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
        protein_id="ENSP00000309572",
    )
    with pytest.raises(ValueError):
        ensembl100.sequence(from_position)


def test_str_mutli_position():
    from_position = ProteinPosition(
        contig_id="5",
        start=525,
        start_offset=0,
        end=526,
        end_offset=0,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
        protein_id="ENSP00000309572",
    )
    assert str(from_position) == "ENSP00000309572:p.525_526"


def test_str_one_position():
    from_position = ProteinPosition(
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
    assert str(from_position) == "ENSP00000309572:p.525"


def test_to_cdna_protein_end_minus_strand(ensembl100):
    from_position = ProteinPosition(
        contig_id="5",
        start=1132,
        start_offset=0,
        end=1132,
        end_offset=0,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
        protein_id="ENSP00000309572",
    )
    to_positions = [
        CdnaPosition(
            contig_id="5",
            start=3394,
            start_offset=0,
            end=3396,
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


def test_to_cdna_protein_end_plus_strand(ensembl100):
    from_position = ProteinPosition(
        contig_id="13",
        start=3418,
        start_offset=0,
        end=3418,
        end_offset=0,
        strand="+",
        gene_id="ENSG00000139618",
        gene_name="BRCA2",
        transcript_id="ENST00000380152",
        transcript_name="BRCA2-201",
        protein_id="ENSP00000369497",
    )
    to_positions = [
        CdnaPosition(
            contig_id="13",
            start=10252,
            start_offset=0,
            end=10254,
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


def test_to_cdna_protein_exon_boundary_minus_strand(ensembl100):
    from_position = ProteinPosition(
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


def test_to_cdna_protein_exon_boundary_plus_strand(ensembl100):
    from_position = ProteinPosition(
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


def test_to_cdna_protein_start_minus_strand(ensembl100):
    from_position = ProteinPosition(
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
        protein_id="ENSP00000309572",
    )
    to_positions = [
        CdnaPosition(
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
            protein_id="ENSP00000309572",
        )
    ]
    assert ensembl100.to_cdna(from_position) == to_positions


def test_to_cdna_protein_start_plus_strand(ensembl100):
    from_position = ProteinPosition(
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
        protein_id="ENSP00000369497",
    )
    to_positions = [
        CdnaPosition(
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
            protein_id="ENSP00000369497",
        )
    ]
    assert ensembl100.to_cdna(from_position) == to_positions


def test_to_dna_protein_end_minus_strand(ensembl100):
    from_position = ProteinPosition(
        contig_id="5",
        start=1132,
        start_offset=0,
        end=1132,
        end_offset=0,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
        protein_id="ENSP00000309572",
    )
    to_positions = [
        DnaPosition(
            contig_id="5", start=1253731, start_offset=0, end=1253733, end_offset=0, strand="-"
        )
    ]
    assert ensembl100.to_dna(from_position) == to_positions


def test_to_dna_protein_end_plus_strand(ensembl100):
    from_position = ProteinPosition(
        contig_id="13",
        start=3418,
        start_offset=0,
        end=3418,
        end_offset=0,
        strand="+",
        gene_id="ENSG00000139618",
        gene_name="BRCA2",
        transcript_id="ENST00000380152",
        transcript_name="BRCA2-201",
        protein_id="ENSP00000369497",
    )
    to_positions = [
        DnaPosition(
            contig_id="13", start=32398765, start_offset=0, end=32398767, end_offset=0, strand="+"
        )
    ]
    assert ensembl100.to_dna(from_position) == to_positions


def test_to_dna_protein_exon_boundary_minus_strand(ensembl100):
    from_position = ProteinPosition(
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
    to_positions = [
        DnaPosition(
            contig_id="5", start=1282623, start_offset=0, end=1293313, end_offset=0, strand="-"
        )
    ]
    assert ensembl100.to_dna(from_position) == to_positions


def test_to_dna_protein_exon_boundary_plus_strand(ensembl100):
    from_position = ProteinPosition(
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
    to_positions = [
        DnaPosition(
            contig_id="13", start=32316527, start_offset=0, end=32319078, end_offset=0, strand="+"
        )
    ]
    assert ensembl100.to_dna(from_position) == to_positions


def test_to_dna_protein_start_minus_strand(ensembl100):
    from_position = ProteinPosition(
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
        protein_id="ENSP00000309572",
    )
    to_positions = [
        DnaPosition(
            contig_id="5", start=1294987, start_offset=0, end=1294989, end_offset=0, strand="-"
        )
    ]
    assert ensembl100.to_dna(from_position) == to_positions


def test_to_dna_protein_start_plus_strand(ensembl100):
    from_position = ProteinPosition(
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
        protein_id="ENSP00000369497",
    )
    to_positions = [
        DnaPosition(
            contig_id="13", start=32316461, start_offset=0, end=32316463, end_offset=0, strand="+"
        )
    ]
    assert ensembl100.to_dna(from_position) == to_positions


def test_to_exon_protein_end_minus_strand(ensembl100):
    from_position = ProteinPosition(
        contig_id="5",
        start=1132,
        start_offset=0,
        end=1132,
        end_offset=0,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
        protein_id="ENSP00000309572",
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


def test_to_exon_protein_end_plus_strand(ensembl100):
    from_position = ProteinPosition(
        contig_id="13",
        start=3418,
        start_offset=0,
        end=3418,
        end_offset=0,
        strand="+",
        gene_id="ENSG00000139618",
        gene_name="BRCA2",
        transcript_id="ENST00000380152",
        transcript_name="BRCA2-201",
        protein_id="ENSP00000369497",
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


def test_to_exon_protein_exon_boundary_minus_strand(ensembl100):
    from_position = ProteinPosition(
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


def test_to_exon_protein_exon_boundary_plus_strand(ensembl100):
    from_position = ProteinPosition(
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


def test_to_exon_protein_start_minus_strand(ensembl100):
    from_position = ProteinPosition(
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
        protein_id="ENSP00000309572",
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


def test_to_exon_protein_start_plus_strand(ensembl100):
    from_position = ProteinPosition(
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
        protein_id="ENSP00000369497",
    )
    to_positions = [
        ExonPosition(
            contig_id="13",
            start=2,
            start_offset=0,
            end=2,
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


def test_to_protein_protein_end_minus_strand(ensembl100):
    from_position = ProteinPosition(
        contig_id="5",
        start=1132,
        start_offset=0,
        end=1132,
        end_offset=0,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
        protein_id="ENSP00000309572",
    )
    to_positions = [
        ProteinPosition(
            contig_id="5",
            start=1132,
            start_offset=0,
            end=1132,
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


def test_to_protein_protein_end_plus_strand(ensembl100):
    from_position = ProteinPosition(
        contig_id="13",
        start=3418,
        start_offset=0,
        end=3418,
        end_offset=0,
        strand="+",
        gene_id="ENSG00000139618",
        gene_name="BRCA2",
        transcript_id="ENST00000380152",
        transcript_name="BRCA2-201",
        protein_id="ENSP00000369497",
    )
    to_positions = [
        ProteinPosition(
            contig_id="13",
            start=3418,
            start_offset=0,
            end=3418,
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


def test_to_protein_protein_exon_boundary_minus_strand(ensembl100):
    from_position = ProteinPosition(
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


def test_to_protein_protein_exon_boundary_plus_strand(ensembl100):
    from_position = ProteinPosition(
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


def test_to_protein_protein_start_minus_strand(ensembl100):
    from_position = ProteinPosition(
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
        protein_id="ENSP00000309572",
    )
    to_positions = [
        ProteinPosition(
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
            protein_id="ENSP00000309572",
        )
    ]
    assert ensembl100.to_protein(from_position) == to_positions


def test_to_protein_protein_start_plus_strand(ensembl100):
    from_position = ProteinPosition(
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
        protein_id="ENSP00000369497",
    )
    to_positions = [
        ProteinPosition(
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
            protein_id="ENSP00000369497",
        )
    ]
    assert ensembl100.to_protein(from_position) == to_positions


def test_to_rna_protein_end_minus_strand(ensembl100):
    from_position = ProteinPosition(
        contig_id="5",
        start=1132,
        start_offset=0,
        end=1132,
        end_offset=0,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
        protein_id="ENSP00000309572",
    )
    to_positions = [
        RnaPosition(
            contig_id="5",
            start=3473,
            start_offset=0,
            end=3475,
            end_offset=0,
            strand="-",
            gene_id="ENSG00000164362",
            gene_name="TERT",
            transcript_id="ENST00000310581",
            transcript_name="TERT-201",
        )
    ]
    assert ensembl100.to_rna(from_position) == to_positions


def test_to_rna_protein_end_plus_strand(ensembl100):
    from_position = ProteinPosition(
        contig_id="13",
        start=3418,
        start_offset=0,
        end=3418,
        end_offset=0,
        strand="+",
        gene_id="ENSG00000139618",
        gene_name="BRCA2",
        transcript_id="ENST00000380152",
        transcript_name="BRCA2-201",
        protein_id="ENSP00000369497",
    )
    to_positions = [
        RnaPosition(
            contig_id="13",
            start=10485,
            start_offset=0,
            end=10487,
            end_offset=0,
            strand="+",
            gene_id="ENSG00000139618",
            gene_name="BRCA2",
            transcript_id="ENST00000380152",
            transcript_name="BRCA2-201",
        )
    ]
    assert ensembl100.to_rna(from_position) == to_positions


def test_to_rna_protein_exon_boundary_minus_strand(ensembl100):
    from_position = ProteinPosition(
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


def test_to_rna_protein_exon_boundary_plus_strand(ensembl100):
    from_position = ProteinPosition(
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


def test_to_rna_protein_start_minus_strand(ensembl100):
    from_position = ProteinPosition(
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
        protein_id="ENSP00000309572",
    )
    to_positions = [
        RnaPosition(
            contig_id="5",
            start=80,
            start_offset=0,
            end=82,
            end_offset=0,
            strand="-",
            gene_id="ENSG00000164362",
            gene_name="TERT",
            transcript_id="ENST00000310581",
            transcript_name="TERT-201",
        )
    ]
    assert ensembl100.to_rna(from_position) == to_positions


def test_to_rna_protein_start_plus_strand(ensembl100):
    from_position = ProteinPosition(
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
        protein_id="ENSP00000369497",
    )
    to_positions = [
        RnaPosition(
            contig_id="13",
            start=234,
            start_offset=0,
            end=236,
            end_offset=0,
            strand="+",
            gene_id="ENSG00000139618",
            gene_name="BRCA2",
            transcript_id="ENST00000380152",
            transcript_name="BRCA2-201",
        )
    ]
    assert ensembl100.to_rna(from_position) == to_positions


def test_to_cdna_offset_error(ensembl100):
    from_position = ProteinPosition(
        contig_id="13",
        start=1,
        start_offset=1,
        end=1,
        end_offset=1,
        strand="+",
        gene_id="ENSG00000139618",
        gene_name="BRCA2",
        transcript_id="ENST00000380152",
        transcript_name="BRCA2-201",
        protein_id="ENSP00000369497",
    )
    with pytest.raises(AssertionError):
        ensembl100.to_cdna(from_position)


def test_to_dna_offset_error(ensembl100):
    from_position = ProteinPosition(
        contig_id="13",
        start=1,
        start_offset=1,
        end=1,
        end_offset=1,
        strand="+",
        gene_id="ENSG00000139618",
        gene_name="BRCA2",
        transcript_id="ENST00000380152",
        transcript_name="BRCA2-201",
        protein_id="ENSP00000369497",
    )
    with pytest.raises(AssertionError):
        ensembl100.to_dna(from_position)


def test_to_exon_offset_error(ensembl100):
    from_position = ProteinPosition(
        contig_id="13",
        start=1,
        start_offset=1,
        end=1,
        end_offset=1,
        strand="+",
        gene_id="ENSG00000139618",
        gene_name="BRCA2",
        transcript_id="ENST00000380152",
        transcript_name="BRCA2-201",
        protein_id="ENSP00000369497",
    )
    with pytest.raises(AssertionError):
        ensembl100.to_exon(from_position)


def test_to_protein_offset_error(ensembl100):
    from_position = ProteinPosition(
        contig_id="13",
        start=1,
        start_offset=1,
        end=1,
        end_offset=1,
        strand="+",
        gene_id="ENSG00000139618",
        gene_name="BRCA2",
        transcript_id="ENST00000380152",
        transcript_name="BRCA2-201",
        protein_id="ENSP00000369497",
    )
    with pytest.raises(AssertionError):
        ensembl100.to_protein(from_position)


def test_to_rna_offset_error(ensembl100):
    from_position = ProteinPosition(
        contig_id="13",
        start=1,
        start_offset=1,
        end=1,
        end_offset=1,
        strand="+",
        gene_id="ENSG00000139618",
        gene_name="BRCA2",
        transcript_id="ENST00000380152",
        transcript_name="BRCA2-201",
        protein_id="ENSP00000369497",
    )
    with pytest.raises(AssertionError):
        ensembl100.to_rna(from_position)
