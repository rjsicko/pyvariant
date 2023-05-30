import pytest

from pyvariant.variants import CdnaPosition, DnaPosition, ExonPosition, ProteinPosition, RnaPosition


def test_is_cdna(ensembl100):
    from_position = ExonPosition(
        _core=ensembl100,
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
    assert from_position.is_cdna is False


def test_is_dna(ensembl100):
    from_position = ExonPosition(
        _core=ensembl100,
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
    assert from_position.is_dna is False


def test_is_exon(ensembl100):
    from_position = ExonPosition(
        _core=ensembl100,
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
    assert from_position.is_exon is True


def test_is_protein(ensembl100):
    from_position = ExonPosition(
        _core=ensembl100,
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
    assert from_position.is_protein is False


def test_is_rna(ensembl100):
    from_position = ExonPosition(
        _core=ensembl100,
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
    assert from_position.is_rna is False


def test_is_on_negative_strand(ensembl100):
    from_position = ExonPosition(
        _core=ensembl100,
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
    assert from_position.on_negative_strand is True


def test_is_on_positive_strand(ensembl100):
    from_position = ExonPosition(
        _core=ensembl100,
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
    assert from_position.on_positive_strand is False


def test_str_mutli_position(ensembl100):
    from_position = ExonPosition(
        _core=ensembl100,
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
    assert str(from_position) == "ENST00000310581:e.2_3"


def test_str_one_position(ensembl100):
    from_position = ExonPosition(
        _core=ensembl100,
        contig_id="5",
        start=2,
        start_offset=0,
        end=2,
        end_offset=0,
        strand="-",
        gene_id="ENSG00000164362",
        gene_name="TERT",
        transcript_id="ENST00000310581",
        transcript_name="TERT-201",
        exon_id="ENSE00001197112",
    )
    assert str(from_position) == "ENST00000310581:e.2"


def test_to_cdna_exon_end_minus_strand(ensembl100):
    from_position = ExonPosition(
        _core=ensembl100,
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
    to_positions = [
        CdnaPosition(
            _core=ensembl100,
            contig_id="5",
            start=3296,
            start_offset=0,
            end=3396,
            end_offset=0,
            strand="-",
            gene_id="ENSG00000164362",
            gene_name="TERT",
            transcript_id="ENST00000310581",
            transcript_name="TERT-201",
            protein_id="ENSP00000309572",
        ),
        CdnaPosition(
            _core=ensembl100,
            contig_id="5",
            start=3397,
            start_offset=0,
            end=3399,
            end_offset=0,
            strand="-",
            gene_id="ENSG00000164362",
            gene_name="TERT",
            transcript_id="ENST00000310581",
            transcript_name="TERT-201",
            protein_id="ENSP00000309572",
        ),
    ]
    assert ensembl100.to_cdna(from_position) == to_positions


def test_to_cdna_exon_end_plus_strand(ensembl100):
    from_position = ExonPosition(
        _core=ensembl100,
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
        exon_id="ENSE00001863787",
    )
    to_positions = [
        CdnaPosition(
            _core=ensembl100,
            contig_id="13",
            start=10255,
            start_offset=0,
            end=10257,
            end_offset=0,
            strand="+",
            gene_id="ENSG00000139618",
            gene_name="BRCA2",
            transcript_id="ENST00000380152",
            transcript_name="BRCA2-201",
            protein_id="ENSP00000369497",
        ),
        CdnaPosition(
            _core=ensembl100,
            contig_id="13",
            start=9649,
            start_offset=0,
            end=10254,
            end_offset=0,
            strand="+",
            gene_id="ENSG00000139618",
            gene_name="BRCA2",
            transcript_id="ENST00000380152",
            transcript_name="BRCA2-201",
            protein_id="ENSP00000369497",
        ),
    ]
    assert ensembl100.to_cdna(from_position) == to_positions


def test_to_cdna_exon_exon_boundary_minus_strand(ensembl100):
    from_position = ExonPosition(
        _core=ensembl100,
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
    to_positions = [
        CdnaPosition(
            _core=ensembl100,
            contig_id="5",
            start=220,
            start_offset=0,
            end=1769,
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


def test_to_cdna_exon_exon_boundary_plus_strand(ensembl100):
    from_position = ExonPosition(
        _core=ensembl100,
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
    to_positions = [
        CdnaPosition(
            _core=ensembl100,
            contig_id="13",
            start=1,
            start_offset=0,
            end=316,
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


def test_to_cdna_exon_overlapping_genes_different_strands_minus_strand(ensembl100):
    from_position = ExonPosition(
        _core=ensembl100,
        contig_id="6",
        start=1,
        start_offset=0,
        end=1,
        end_offset=0,
        strand="-",
        gene_id="ENSG00000204531",
        gene_name="POU5F1",
        transcript_id="ENST00000471529",
        transcript_name="POU5F1-204",
        exon_id="ENSE00002568331",
    )
    to_positions = []
    assert ensembl100.to_cdna(from_position) == to_positions


def test_to_cdna_exon_overlapping_genes_different_strands_plus_strand(ensembl100):
    from_position = ExonPosition(
        _core=ensembl100,
        contig_id="6",
        start=2,
        start_offset=0,
        end=2,
        end_offset=0,
        strand="+",
        gene_id="ENSG00000137310",
        gene_name="TCF19",
        transcript_id="ENST00000542218",
        transcript_name="TCF19-204",
        exon_id="ENSE00002283659",
    )
    to_positions = [
        CdnaPosition(
            _core=ensembl100,
            contig_id="6",
            start=558,
            start_offset=0,
            end=643,
            end_offset=0,
            strand="+",
            gene_id="ENSG00000137310",
            gene_name="TCF19",
            transcript_id="ENST00000542218",
            transcript_name="TCF19-204",
            protein_id="ENSP00000439397",
        )
    ]
    assert ensembl100.to_cdna(from_position) == to_positions


def test_to_cdna_exon_start_minus_strand(ensembl100):
    from_position = ExonPosition(
        _core=ensembl100,
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
    to_positions = [
        CdnaPosition(
            _core=ensembl100,
            contig_id="5",
            start=1,
            start_offset=0,
            end=219,
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


def test_to_cdna_exon_start_plus_strand(ensembl100):
    from_position = ExonPosition(
        _core=ensembl100,
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
        exon_id="ENSE00003896691",
    )
    to_positions = []
    assert ensembl100.to_cdna(from_position) == to_positions


def test_to_dna_exon_end_minus_strand(ensembl100):
    from_position = ExonPosition(
        _core=ensembl100,
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
    to_positions = [
        DnaPosition(
            _core=ensembl100,
            contig_id="5",
            start=1253167,
            start_offset=0,
            end=1253831,
            end_offset=0,
            strand="-",
        )
    ]
    assert ensembl100.to_dna(from_position) == to_positions


def test_to_dna_exon_end_plus_strand(ensembl100):
    from_position = ExonPosition(
        _core=ensembl100,
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
        exon_id="ENSE00001863787",
    )
    to_positions = [
        DnaPosition(
            _core=ensembl100,
            contig_id="13",
            start=32398162,
            start_offset=0,
            end=32400266,
            end_offset=0,
            strand="+",
        )
    ]
    assert ensembl100.to_dna(from_position) == to_positions


def test_to_dna_exon_exon_boundary_minus_strand(ensembl100):
    from_position = ExonPosition(
        _core=ensembl100,
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
    to_positions = [
        DnaPosition(
            _core=ensembl100,
            contig_id="5",
            start=1282429,
            start_offset=0,
            end=1294666,
            end_offset=0,
            strand="-",
        )
    ]
    assert ensembl100.to_dna(from_position) == to_positions


def test_to_dna_exon_exon_boundary_plus_strand(ensembl100):
    from_position = ExonPosition(
        _core=ensembl100,
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
    to_positions = [
        DnaPosition(
            _core=ensembl100,
            contig_id="13",
            start=32316422,
            start_offset=0,
            end=32319325,
            end_offset=0,
            strand="+",
        )
    ]
    assert ensembl100.to_dna(from_position) == to_positions


def test_to_dna_exon_overlapping_genes_different_strands_minus_strand(ensembl100):
    from_position = ExonPosition(
        _core=ensembl100,
        contig_id="6",
        start=1,
        start_offset=0,
        end=1,
        end_offset=0,
        strand="-",
        gene_id="ENSG00000204531",
        gene_name="POU5F1",
        transcript_id="ENST00000471529",
        transcript_name="POU5F1-204",
        exon_id="ENSE00002568331",
    )
    to_positions = [
        DnaPosition(
            _core=ensembl100,
            contig_id="6",
            start=31165927,
            start_offset=0,
            end=31166838,
            end_offset=0,
            strand="-",
        )
    ]
    assert ensembl100.to_dna(from_position) == to_positions


def test_to_dna_exon_overlapping_genes_different_strands_plus_strand(ensembl100):
    from_position = ExonPosition(
        _core=ensembl100,
        contig_id="6",
        start=2,
        start_offset=0,
        end=2,
        end_offset=0,
        strand="+",
        gene_id="ENSG00000137310",
        gene_name="TCF19",
        transcript_id="ENST00000542218",
        transcript_name="TCF19-204",
        exon_id="ENSE00002283659",
    )
    to_positions = [
        DnaPosition(
            _core=ensembl100,
            contig_id="6",
            start=31166280,
            start_offset=0,
            end=31166365,
            end_offset=0,
            strand="+",
        )
    ]
    assert ensembl100.to_dna(from_position) == to_positions


def test_to_dna_exon_start_minus_strand(ensembl100):
    from_position = ExonPosition(
        _core=ensembl100,
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
    to_positions = [
        DnaPosition(
            _core=ensembl100,
            contig_id="5",
            start=1294771,
            start_offset=0,
            end=1295068,
            end_offset=0,
            strand="-",
        )
    ]
    assert ensembl100.to_dna(from_position) == to_positions


def test_to_dna_exon_start_plus_strand(ensembl100):
    from_position = ExonPosition(
        _core=ensembl100,
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
        exon_id="ENSE00003896691",
    )
    to_positions = [
        DnaPosition(
            _core=ensembl100,
            contig_id="13",
            start=32315474,
            start_offset=0,
            end=32315667,
            end_offset=0,
            strand="+",
        )
    ]
    assert ensembl100.to_dna(from_position) == to_positions


def test_to_exon_exon_end_minus_strand(ensembl100):
    from_position = ExonPosition(
        _core=ensembl100,
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
    to_positions = [
        ExonPosition(
            _core=ensembl100,
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


def test_to_exon_exon_end_plus_strand(ensembl100):
    from_position = ExonPosition(
        _core=ensembl100,
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
        exon_id="ENSE00001863787",
    )
    to_positions = [
        ExonPosition(
            _core=ensembl100,
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


def test_to_exon_exon_exon_boundary_minus_strand(ensembl100):
    from_position = ExonPosition(
        _core=ensembl100,
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
    to_positions = [
        ExonPosition(
            _core=ensembl100,
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


def test_to_exon_exon_exon_boundary_plus_strand(ensembl100):
    from_position = ExonPosition(
        _core=ensembl100,
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
    to_positions = [
        ExonPosition(
            _core=ensembl100,
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


def test_to_exon_exon_overlapping_genes_different_strands_minus_strand(ensembl100):
    from_position = ExonPosition(
        _core=ensembl100,
        contig_id="6",
        start=1,
        start_offset=0,
        end=1,
        end_offset=0,
        strand="-",
        gene_id="ENSG00000204531",
        gene_name="POU5F1",
        transcript_id="ENST00000471529",
        transcript_name="POU5F1-204",
        exon_id="ENSE00002568331",
    )
    to_positions = [
        ExonPosition(
            _core=ensembl100,
            contig_id="6",
            start=1,
            start_offset=0,
            end=1,
            end_offset=0,
            strand="-",
            gene_id="ENSG00000204531",
            gene_name="POU5F1",
            transcript_id="ENST00000471529",
            transcript_name="POU5F1-204",
            exon_id="ENSE00002568331",
        )
    ]
    assert ensembl100.to_exon(from_position) == to_positions


def test_to_exon_exon_overlapping_genes_different_strands_plus_strand(ensembl100):
    from_position = ExonPosition(
        _core=ensembl100,
        contig_id="6",
        start=2,
        start_offset=0,
        end=2,
        end_offset=0,
        strand="+",
        gene_id="ENSG00000137310",
        gene_name="TCF19",
        transcript_id="ENST00000542218",
        transcript_name="TCF19-204",
        exon_id="ENSE00002283659",
    )
    to_positions = [
        ExonPosition(
            _core=ensembl100,
            contig_id="6",
            start=2,
            start_offset=0,
            end=2,
            end_offset=0,
            strand="+",
            gene_id="ENSG00000137310",
            gene_name="TCF19",
            transcript_id="ENST00000542218",
            transcript_name="TCF19-204",
            exon_id="ENSE00002283659",
        )
    ]
    assert ensembl100.to_exon(from_position) == to_positions


def test_to_exon_exon_start_minus_strand(ensembl100):
    from_position = ExonPosition(
        _core=ensembl100,
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
    to_positions = [
        ExonPosition(
            _core=ensembl100,
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


def test_to_exon_exon_start_plus_strand(ensembl100):
    from_position = ExonPosition(
        _core=ensembl100,
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
        exon_id="ENSE00003896691",
    )
    to_positions = [
        ExonPosition(
            _core=ensembl100,
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


def test_to_protein_exon_end_minus_strand(ensembl100):
    from_position = ExonPosition(
        _core=ensembl100,
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
    to_positions = [
        ProteinPosition(
            _core=ensembl100,
            contig_id="5",
            start=1099,
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


def test_to_protein_exon_end_plus_strand(ensembl100):
    from_position = ExonPosition(
        _core=ensembl100,
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
        exon_id="ENSE00001863787",
    )
    to_positions = [
        ProteinPosition(
            _core=ensembl100,
            contig_id="13",
            start=3217,
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


def test_to_protein_exon_exon_boundary_minus_strand(ensembl100):
    from_position = ExonPosition(
        _core=ensembl100,
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
    to_positions = [
        ProteinPosition(
            _core=ensembl100,
            contig_id="5",
            start=74,
            start_offset=0,
            end=590,
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


def test_to_protein_exon_exon_boundary_plus_strand(ensembl100):
    from_position = ExonPosition(
        _core=ensembl100,
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
    to_positions = [
        ProteinPosition(
            _core=ensembl100,
            contig_id="13",
            start=1,
            start_offset=0,
            end=106,
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


def test_to_protein_exon_overlapping_genes_different_strands_minus_strand(ensembl100):
    from_position = ExonPosition(
        _core=ensembl100,
        contig_id="6",
        start=1,
        start_offset=0,
        end=1,
        end_offset=0,
        strand="-",
        gene_id="ENSG00000204531",
        gene_name="POU5F1",
        transcript_id="ENST00000471529",
        transcript_name="POU5F1-204",
        exon_id="ENSE00002568331",
    )
    to_positions = []
    assert ensembl100.to_protein(from_position) == to_positions


def test_to_protein_exon_overlapping_genes_different_strands_plus_strand(ensembl100):
    from_position = ExonPosition(
        _core=ensembl100,
        contig_id="6",
        start=2,
        start_offset=0,
        end=2,
        end_offset=0,
        strand="+",
        gene_id="ENSG00000137310",
        gene_name="TCF19",
        transcript_id="ENST00000542218",
        transcript_name="TCF19-204",
        exon_id="ENSE00002283659",
    )
    to_positions = [
        ProteinPosition(
            _core=ensembl100,
            contig_id="6",
            start=186,
            start_offset=0,
            end=215,
            end_offset=0,
            strand="+",
            gene_id="ENSG00000137310",
            gene_name="TCF19",
            transcript_id="ENST00000542218",
            transcript_name="TCF19-204",
            protein_id="ENSP00000439397",
        )
    ]
    assert ensembl100.to_protein(from_position) == to_positions


def test_to_protein_exon_start_minus_strand(ensembl100):
    from_position = ExonPosition(
        _core=ensembl100,
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
    to_positions = [
        ProteinPosition(
            _core=ensembl100,
            contig_id="5",
            start=1,
            start_offset=0,
            end=73,
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


def test_to_protein_exon_start_plus_strand(ensembl100):
    from_position = ExonPosition(
        _core=ensembl100,
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
        exon_id="ENSE00003896691",
    )
    to_positions = []
    assert ensembl100.to_protein(from_position) == to_positions


def test_to_rna_exon_end_minus_strand(ensembl100):
    from_position = ExonPosition(
        _core=ensembl100,
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
    to_positions = [
        RnaPosition(
            _core=ensembl100,
            contig_id="5",
            start=3375,
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


def test_to_rna_exon_end_plus_strand(ensembl100):
    from_position = ExonPosition(
        _core=ensembl100,
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
        exon_id="ENSE00001863787",
    )
    to_positions = [
        RnaPosition(
            _core=ensembl100,
            contig_id="13",
            start=9882,
            start_offset=0,
            end=11986,
            end_offset=0,
            strand="+",
            gene_id="ENSG00000139618",
            gene_name="BRCA2",
            transcript_id="ENST00000380152",
            transcript_name="BRCA2-201",
        )
    ]
    assert ensembl100.to_rna(from_position) == to_positions


def test_to_rna_exon_exon_boundary_minus_strand(ensembl100):
    from_position = ExonPosition(
        _core=ensembl100,
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
    to_positions = [
        RnaPosition(
            _core=ensembl100,
            contig_id="5",
            start=299,
            start_offset=0,
            end=1848,
            end_offset=0,
            strand="-",
            gene_id="ENSG00000164362",
            gene_name="TERT",
            transcript_id="ENST00000310581",
            transcript_name="TERT-201",
        )
    ]
    assert ensembl100.to_rna(from_position) == to_positions


def test_to_rna_exon_exon_boundary_plus_strand(ensembl100):
    from_position = ExonPosition(
        _core=ensembl100,
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
    to_positions = [
        RnaPosition(
            _core=ensembl100,
            contig_id="13",
            start=195,
            start_offset=0,
            end=549,
            end_offset=0,
            strand="+",
            gene_id="ENSG00000139618",
            gene_name="BRCA2",
            transcript_id="ENST00000380152",
            transcript_name="BRCA2-201",
        )
    ]
    assert ensembl100.to_rna(from_position) == to_positions


def test_to_rna_exon_overlapping_genes_different_strands_minus_strand(ensembl100):
    from_position = ExonPosition(
        _core=ensembl100,
        contig_id="6",
        start=1,
        start_offset=0,
        end=1,
        end_offset=0,
        strand="-",
        gene_id="ENSG00000204531",
        gene_name="POU5F1",
        transcript_id="ENST00000471529",
        transcript_name="POU5F1-204",
        exon_id="ENSE00002568331",
    )
    to_positions = [
        RnaPosition(
            _core=ensembl100,
            contig_id="6",
            start=1,
            start_offset=0,
            end=912,
            end_offset=0,
            strand="-",
            gene_id="ENSG00000204531",
            gene_name="POU5F1",
            transcript_id="ENST00000471529",
            transcript_name="POU5F1-204",
        )
    ]
    assert ensembl100.to_rna(from_position) == to_positions


def test_to_rna_exon_overlapping_genes_different_strands_plus_strand(ensembl100):
    from_position = ExonPosition(
        _core=ensembl100,
        contig_id="6",
        start=2,
        start_offset=0,
        end=2,
        end_offset=0,
        strand="+",
        gene_id="ENSG00000137310",
        gene_name="TCF19",
        transcript_id="ENST00000542218",
        transcript_name="TCF19-204",
        exon_id="ENSE00002283659",
    )
    to_positions = [
        RnaPosition(
            _core=ensembl100,
            contig_id="6",
            start=558,
            start_offset=0,
            end=643,
            end_offset=0,
            strand="+",
            gene_id="ENSG00000137310",
            gene_name="TCF19",
            transcript_id="ENST00000542218",
            transcript_name="TCF19-204",
        )
    ]
    assert ensembl100.to_rna(from_position) == to_positions


def test_to_rna_exon_start_minus_strand(ensembl100):
    from_position = ExonPosition(
        _core=ensembl100,
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
    to_positions = [
        RnaPosition(
            _core=ensembl100,
            contig_id="5",
            start=1,
            start_offset=0,
            end=298,
            end_offset=0,
            strand="-",
            gene_id="ENSG00000164362",
            gene_name="TERT",
            transcript_id="ENST00000310581",
            transcript_name="TERT-201",
        )
    ]
    assert ensembl100.to_rna(from_position) == to_positions


def test_to_rna_exon_start_plus_strand(ensembl100):
    from_position = ExonPosition(
        _core=ensembl100,
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
        exon_id="ENSE00003896691",
    )
    to_positions = [
        RnaPosition(
            _core=ensembl100,
            contig_id="13",
            start=1,
            start_offset=0,
            end=194,
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
    from_position = ExonPosition(
        _core=ensembl100,
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
        exon_id="ENSE00003896691",
    )
    with pytest.raises(AssertionError):
        ensembl100.to_cdna(from_position)


def test_to_dna_offset_error(ensembl100):
    from_position = ExonPosition(
        _core=ensembl100,
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
        exon_id="ENSE00003896691",
    )
    with pytest.raises(AssertionError):
        ensembl100.to_dna(from_position)


def test_to_exon_offset_error(ensembl100):
    from_position = ExonPosition(
        _core=ensembl100,
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
        exon_id="ENSE00003896691",
    )
    with pytest.raises(AssertionError):
        ensembl100.to_exon(from_position)


def test_to_protein_offset_error(ensembl100):
    from_position = ExonPosition(
        _core=ensembl100,
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
        exon_id="ENSE00003896691",
    )
    with pytest.raises(AssertionError):
        ensembl100.to_protein(from_position)


def test_to_rna_offset_error(ensembl100):
    from_position = ExonPosition(
        _core=ensembl100,
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
        exon_id="ENSE00003896691",
    )
    with pytest.raises(AssertionError):
        ensembl100.to_rna(from_position)


def test_sequence(ensembl100):
    from_position = ExonPosition(
        _core=ensembl100,
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
        exon_id="ENSE00003896691",
    )
    with pytest.raises(NotImplementedError):
        ensembl100.sequence(from_position)
