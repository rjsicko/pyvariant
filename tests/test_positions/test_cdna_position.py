import pytest

from variant_map.positions import (
    CdnaPosition,
    DnaPosition,
    ExonPosition,
    ProteinPosition,
    RnaPosition,
)


# TODO: use fixtures for repeated positions
def test_is_cdna():
    from_position = CdnaPosition(
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
    assert from_position.is_cdna is True


def test_is_dna():
    from_position = CdnaPosition(
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
    assert from_position.is_dna is False


def test_is_exon():
    from_position = CdnaPosition(
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
    assert from_position.is_exon is False


def test_is_protein():
    from_position = CdnaPosition(
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
    assert from_position.is_protein is False


def test_is_rna():
    from_position = CdnaPosition(
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
    assert from_position.is_rna is False


def test_is_on_negative_strand():
    from_position = CdnaPosition(
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
    assert from_position.on_negative_strand is True


def test_is_on_positive_strand():
    from_position = CdnaPosition(
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
    assert from_position.on_positive_strand is False


def test_sequence(ensembl100):
    from_position = CdnaPosition(
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
    assert ensembl100.sequence(from_position) == "GGG"


@pytest.mark.skip(reason="how to handle?")
def test_sequence_offset(ensembl100):
    # Position crosses exon boundary so getting DNA sequence returns intron sequences as well
    # TODO: How to handle sequences for offset variants?
    from_position = CdnaPosition(
        contig_id="5",
        start=1572,
        start_offset=1,
        end=1574,
        end_offset=1,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
        protein_id="ENSP00000309572",
    )
    assert ensembl100.sequence(from_position) == "GGG"


def test_str_mutli_position():
    from_position = CdnaPosition(
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
    assert str(from_position) == "ENST00000310581:c.1573_1575"


def test_str_one_position():
    from_position = CdnaPosition(
        contig_id="5",
        start=1573,
        start_offset=0,
        end=1573,
        end_offset=0,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
        protein_id="ENSP00000309572",
    )
    assert str(from_position) == "ENST00000310581:c.1573"


def test_to_cdna_cdna_end_minus_strand(ensembl100):
    from_position = CdnaPosition(
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


def test_to_cdna_cdna_end_plus_strand(ensembl100):
    from_position = CdnaPosition(
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


def test_to_cdna_cdna_exon_boundary_minus_strand(ensembl100):
    from_position = CdnaPosition(
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


def test_to_cdna_cdna_exon_boundary_plus_strand(ensembl100):
    from_position = CdnaPosition(
        contig_id="13",
        start=1573,
        start_offset=0,
        end=1575,
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
            start=1573,
            start_offset=0,
            end=1575,
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


def test_to_cdna_cdna_negative_offset_into_intron_plus_strand(ensembl100):
    from_position = CdnaPosition(
        contig_id="4",
        start=68,
        start_offset=-1,
        end=68,
        end_offset=-1,
        strand="+",
        gene_id="ENSG00000157404",
        gene_name="KIT",
        transcript_id="ENST00000288135",
        transcript_name="KIT-201",
        protein_id="ENSP00000288135",
    )
    to_positions = [
        CdnaPosition(
            contig_id="4",
            start=68,
            start_offset=-1,
            end=68,
            end_offset=-1,
            strand="+",
            gene_id="ENSG00000157404",
            gene_name="KIT",
            transcript_id="ENST00000288135",
            transcript_name="KIT-201",
            protein_id="ENSP00000288135",
        )
    ]
    assert ensembl100.to_cdna(from_position) == to_positions


def test_to_cdna_cdna_negative_offset_normalized_minus_strand(ensembl100):
    from_position = CdnaPosition(
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
        protein_id="ENSP00000309572",
    )
    to_positions = [
        CdnaPosition(
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
            protein_id="ENSP00000309572",
        )
    ]
    assert ensembl100.to_cdna(from_position) == to_positions


def test_to_cdna_cdna_positive_offset_into_intron_plus_strand(ensembl100):
    from_position = CdnaPosition(
        contig_id="4",
        start=67,
        start_offset=1,
        end=67,
        end_offset=1,
        strand="+",
        gene_id="ENSG00000157404",
        gene_name="KIT",
        transcript_id="ENST00000288135",
        transcript_name="KIT-201",
        protein_id="ENSP00000288135",
    )
    to_positions = [
        CdnaPosition(
            contig_id="4",
            start=67,
            start_offset=1,
            end=67,
            end_offset=1,
            strand="+",
            gene_id="ENSG00000157404",
            gene_name="KIT",
            transcript_id="ENST00000288135",
            transcript_name="KIT-201",
            protein_id="ENSP00000288135",
        )
    ]
    assert ensembl100.to_cdna(from_position) == to_positions


def test_to_cdna_cdna_start_minus_strand(ensembl100):
    from_position = CdnaPosition(
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


def test_to_cdna_cdna_start_plus_strand(ensembl100):
    from_position = CdnaPosition(
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


def test_to_dna_cdna_end_minus_strand(ensembl100):
    from_position = CdnaPosition(
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
    to_positions = [
        DnaPosition(
            contig_id="5", start=1253731, start_offset=0, end=1253733, end_offset=0, strand="-"
        )
    ]
    assert ensembl100.to_dna(from_position) == to_positions


def test_to_dna_cdna_end_plus_strand(ensembl100):
    from_position = CdnaPosition(
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
    to_positions = [
        DnaPosition(
            contig_id="13", start=32398765, start_offset=0, end=32398767, end_offset=0, strand="+"
        )
    ]
    assert ensembl100.to_dna(from_position) == to_positions


def test_to_dna_cdna_exon_boundary_minus_strand(ensembl100):
    from_position = CdnaPosition(
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
    to_positions = [
        DnaPosition(
            contig_id="5", start=1282623, start_offset=0, end=1293313, end_offset=0, strand="-"
        )
    ]
    assert ensembl100.to_dna(from_position) == to_positions


def test_to_dna_cdna_exon_boundary_plus_strand(ensembl100):
    from_position = CdnaPosition(
        contig_id="13",
        start=1573,
        start_offset=0,
        end=1575,
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
            contig_id="13", start=32333051, start_offset=0, end=32333053, end_offset=0, strand="+"
        )
    ]
    assert ensembl100.to_dna(from_position) == to_positions


def test_to_dna_cdna_negative_offset_into_intron_plus_strand(ensembl100):
    from_position = CdnaPosition(
        contig_id="4",
        start=68,
        start_offset=-1,
        end=68,
        end_offset=-1,
        strand="+",
        gene_id="ENSG00000157404",
        gene_name="KIT",
        transcript_id="ENST00000288135",
        transcript_name="KIT-201",
        protein_id="ENSP00000288135",
    )
    to_positions = [
        DnaPosition(
            contig_id="4", start=54695511, start_offset=0, end=54695511, end_offset=0, strand="+"
        )
    ]
    assert ensembl100.to_dna(from_position) == to_positions


def test_to_dna_cdna_negative_offset_normalized_minus_strand(ensembl100):
    from_position = CdnaPosition(
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
        protein_id="ENSP00000309572",
    )
    to_positions = [
        DnaPosition(
            contig_id="5", start=1294891, start_offset=0, end=1294891, end_offset=0, strand="-"
        )
    ]
    assert ensembl100.to_dna(from_position) == to_positions


def test_to_dna_cdna_positive_offset_into_intron_plus_strand(ensembl100):
    from_position = CdnaPosition(
        contig_id="4",
        start=67,
        start_offset=1,
        end=67,
        end_offset=1,
        strand="+",
        gene_id="ENSG00000157404",
        gene_name="KIT",
        transcript_id="ENST00000288135",
        transcript_name="KIT-201",
        protein_id="ENSP00000288135",
    )
    to_positions = [
        DnaPosition(
            contig_id="4", start=54658082, start_offset=0, end=54658082, end_offset=0, strand="+"
        )
    ]
    assert ensembl100.to_dna(from_position) == to_positions


def test_to_dna_cdna_start_minus_strand(ensembl100):
    from_position = CdnaPosition(
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
    to_positions = [
        DnaPosition(
            contig_id="5", start=1294987, start_offset=0, end=1294989, end_offset=0, strand="-"
        )
    ]
    assert ensembl100.to_dna(from_position) == to_positions


def test_to_dna_cdna_start_plus_strand(ensembl100):
    from_position = CdnaPosition(
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
    to_positions = [
        DnaPosition(
            contig_id="13", start=32316461, start_offset=0, end=32316463, end_offset=0, strand="+"
        )
    ]
    assert ensembl100.to_dna(from_position) == to_positions


def test_to_exon_cdna_end_minus_strand(ensembl100):
    from_position = CdnaPosition(
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


def test_to_exon_cdna_end_plus_strand(ensembl100):
    from_position = CdnaPosition(
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


def test_to_exon_cdna_exon_boundary_minus_strand(ensembl100):
    from_position = CdnaPosition(
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


def test_to_exon_cdna_exon_boundary_plus_strand(ensembl100):
    from_position = CdnaPosition(
        contig_id="13",
        start=1573,
        start_offset=0,
        end=1575,
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
            start=10,
            start_offset=0,
            end=10,
            end_offset=0,
            strand="+",
            gene_id="ENSG00000139618",
            gene_name="BRCA2",
            transcript_id="ENST00000380152",
            transcript_name="BRCA2-201",
            exon_id="ENSE00000939167",
        )
    ]
    assert ensembl100.to_exon(from_position) == to_positions


def test_to_exon_cdna_negative_offset_into_intron_plus_strand(ensembl100):
    from_position = CdnaPosition(
        contig_id="4",
        start=68,
        start_offset=-1,
        end=68,
        end_offset=-1,
        strand="+",
        gene_id="ENSG00000157404",
        gene_name="KIT",
        transcript_id="ENST00000288135",
        transcript_name="KIT-201",
        protein_id="ENSP00000288135",
    )
    to_positions = []
    assert ensembl100.to_exon(from_position) == to_positions


def test_to_exon_cdna_negative_offset_normalized_minus_strand(ensembl100):
    from_position = CdnaPosition(
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


def test_to_exon_cdna_positive_offset_into_intron_plus_strand(ensembl100):
    from_position = CdnaPosition(
        contig_id="4",
        start=67,
        start_offset=1,
        end=67,
        end_offset=1,
        strand="+",
        gene_id="ENSG00000157404",
        gene_name="KIT",
        transcript_id="ENST00000288135",
        transcript_name="KIT-201",
        protein_id="ENSP00000288135",
    )
    to_positions = []
    assert ensembl100.to_exon(from_position) == to_positions


def test_to_exon_cdna_start_minus_strand(ensembl100):
    from_position = CdnaPosition(
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


def test_to_exon_cdna_start_plus_strand(ensembl100):
    from_position = CdnaPosition(
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


def test_to_protein_cdna_end_minus_strand(ensembl100):
    from_position = CdnaPosition(
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


def test_to_protein_cdna_end_plus_strand(ensembl100):
    from_position = CdnaPosition(
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


def test_to_protein_cdna_exon_boundary_minus_strand(ensembl100):
    from_position = CdnaPosition(
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


def test_to_protein_cdna_exon_boundary_plus_strand(ensembl100):
    from_position = CdnaPosition(
        contig_id="13",
        start=1573,
        start_offset=0,
        end=1575,
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
            start=525,
            start_offset=0,
            end=525,
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


def test_to_protein_cdna_negative_offset_into_intron_plus_strand(ensembl100):
    from_position = CdnaPosition(
        contig_id="4",
        start=68,
        start_offset=-1,
        end=68,
        end_offset=-1,
        strand="+",
        gene_id="ENSG00000157404",
        gene_name="KIT",
        transcript_id="ENST00000288135",
        transcript_name="KIT-201",
        protein_id="ENSP00000288135",
    )
    to_positions = []
    assert ensembl100.to_protein(from_position) == to_positions


def test_to_protein_cdna_negative_offset_normalized_minus_strand(ensembl100):
    from_position = CdnaPosition(
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
        protein_id="ENSP00000309572",
    )
    to_positions = [
        ProteinPosition(
            contig_id="5",
            start=33,
            start_offset=0,
            end=33,
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


def test_to_protein_cdna_positive_offset_into_intron_plus_strand(ensembl100):
    from_position = CdnaPosition(
        contig_id="4",
        start=67,
        start_offset=1,
        end=67,
        end_offset=1,
        strand="+",
        gene_id="ENSG00000157404",
        gene_name="KIT",
        transcript_id="ENST00000288135",
        transcript_name="KIT-201",
        protein_id="ENSP00000288135",
    )
    to_positions = []
    assert ensembl100.to_protein(from_position) == to_positions


def test_to_protein_cdna_start_minus_strand(ensembl100):
    from_position = CdnaPosition(
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


def test_to_protein_cdna_start_plus_strand(ensembl100):
    from_position = CdnaPosition(
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


def test_to_rna_cdna_end_minus_strand(ensembl100):
    from_position = CdnaPosition(
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


def test_to_rna_cdna_end_plus_strand(ensembl100):
    from_position = CdnaPosition(
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


def test_to_rna_cdna_exon_boundary_minus_strand(ensembl100):
    from_position = CdnaPosition(
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


def test_to_rna_cdna_exon_boundary_plus_strand(ensembl100):
    from_position = CdnaPosition(
        contig_id="13",
        start=1573,
        start_offset=0,
        end=1575,
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
            start=1806,
            start_offset=0,
            end=1808,
            end_offset=0,
            strand="+",
            gene_id="ENSG00000139618",
            gene_name="BRCA2",
            transcript_id="ENST00000380152",
            transcript_name="BRCA2-201",
        )
    ]
    assert ensembl100.to_rna(from_position) == to_positions


def test_to_rna_cdna_negative_offset_into_intron_plus_strand(ensembl100):
    from_position = CdnaPosition(
        contig_id="4",
        start=68,
        start_offset=-1,
        end=68,
        end_offset=-1,
        strand="+",
        gene_id="ENSG00000157404",
        gene_name="KIT",
        transcript_id="ENST00000288135",
        transcript_name="KIT-201",
        protein_id="ENSP00000288135",
    )
    to_positions = [
        RnaPosition(
            contig_id="4",
            start=126,
            start_offset=-1,
            end=126,
            end_offset=-1,
            strand="+",
            gene_id="ENSG00000157404",
            gene_name="KIT",
            transcript_id="ENST00000288135",
            transcript_name="KIT-201",
        )
    ]
    assert ensembl100.to_rna(from_position) == to_positions


def test_to_rna_cdna_negative_offset_normalized_minus_strand(ensembl100):
    from_position = CdnaPosition(
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
        protein_id="ENSP00000309572",
    )
    to_positions = [
        RnaPosition(
            contig_id="5",
            start=178,
            start_offset=0,
            end=178,
            end_offset=0,
            strand="-",
            gene_id="ENSG00000164362",
            gene_name="TERT",
            transcript_id="ENST00000310581",
            transcript_name="TERT-201",
        )
    ]
    assert ensembl100.to_rna(from_position) == to_positions


def test_to_rna_cdna_positive_offset_into_intron_plus_strand(ensembl100):
    from_position = CdnaPosition(
        contig_id="4",
        start=67,
        start_offset=1,
        end=67,
        end_offset=1,
        strand="+",
        gene_id="ENSG00000157404",
        gene_name="KIT",
        transcript_id="ENST00000288135",
        transcript_name="KIT-201",
        protein_id="ENSP00000288135",
    )
    to_positions = [
        RnaPosition(
            contig_id="4",
            start=125,
            start_offset=1,
            end=125,
            end_offset=1,
            strand="+",
            gene_id="ENSG00000157404",
            gene_name="KIT",
            transcript_id="ENST00000288135",
            transcript_name="KIT-201",
        )
    ]
    assert ensembl100.to_rna(from_position) == to_positions


def test_to_rna_cdna_start_minus_strand(ensembl100):
    from_position = CdnaPosition(
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


def test_to_rna_cdna_start_plus_strand(ensembl100):
    from_position = CdnaPosition(
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
