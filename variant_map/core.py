"""Core logic for handling annotations and mapping positions/variants between different types."""
from __future__ import annotations

from itertools import product
from typing import Any, Callable, Dict, List, Optional, Tuple, Type, Union, cast

import pandas as pd
from gtfparse import read_gtf
from pyfaidx import Fasta

from .constants import (
    CDNA,
    CDS,
    CONTIG_ID,
    DNA,
    EXON,
    EXON_ID,
    FUSION,
    GENE_ID,
    GENE_NAME,
    PROTEIN,
    PROTEIN_ID,
    RNA,
    STOP_CODON,
    TRANSCRIPT_ID,
    TRANSCRIPT_NAME,
)
from .files import read_fasta, tsv_to_dict, txt_to_list
from .parser import parse
from .positions import (
    CdnaDeletion,
    CdnaDelins,
    CdnaDuplication,
    CdnaFusion,
    CdnaInsertion,
    CdnaPosition,
    CdnaSubstitution,
    DnaDeletion,
    DnaDelins,
    DnaDuplication,
    DnaFusion,
    DnaInsertion,
    DnaPosition,
    DnaSubstitution,
    ExonFusion,
    ExonPosition,
    ProteinDeletion,
    ProteinDelins,
    ProteinDuplication,
    ProteinFusion,
    ProteinInsertion,
    ProteinPosition,
    ProteinSubstitution,
    RnaDeletion,
    RnaDelins,
    RnaDuplication,
    RnaFusion,
    RnaInsertion,
    RnaPosition,
    RnaSubstitution,
    _CdnaSmallVariant,
    _DnaSmallVariant,
    _ExonSmallVariant,
    _Fusion,
    _ProteinSmallVariant,
    _RnaSmallVariant,
)
from .tables import AMINO_ACID_TABLE
from .utils import (
    calc_cdna_to_protein,
    collapse_seq_change,
    expand_nt,
    is_deletion,
    is_delins,
    is_duplication,
    is_frameshift,
    is_insertion,
    is_substitution,
    reverse_complement,
    reverse_translate,
    split_by_codon,
    strip_version,
)


class Core:
    """Core class that handles converting between position types and retrieving information on
    biological features.
    """

    def __init__(
        self,
        gtf: str,
        cds: List[str],
        dna: List[str],
        peptide: List[str],
        rna: List[str],
        canonical_transcript: Union[List[str], str] = [],
        contig_alias: Union[str, Dict] = {},
        exon_alias: Union[str, Dict] = {},
        gene_alias: Union[str, Dict] = {},
        protein_alias: Union[str, Dict] = {},
        transcript_alias: Union[str, Dict] = {},
    ):
        """_summary_

        Args:
            gtf (str): Path to a GTF files with feature annotations
            cds (List[str]): List of paths to a FASTA files of CDS sequences
            dna (List[str]): List of paths to a FASTA files of DNA sequences
            peptide (List[str]): List of paths to a FASTA files of peptide sequences
            rna (List[str]): List of paths to a FASTA files of RNA sequences
            canonical_transcript (Union[List[str], str], optional): List of canonical transcript IDs, or the path to a text file
            contig_alias (Union[str, Dict], optional): Dictionary mapping contig aliases to their normalized ID, or a path to a text file
            exon_alias (Union[str, Dict], optional): Dictionary mapping exon aliases to their normalized ID, or a path to a text file
            gene_alias (Union[str, Dict], optional): Dictionary mapping gene aliases to their normalized ID, or a path to a text file
            protein_alias (Union[str, Dict], optional): Dictionary mapping protein aliases to their normalized ID, or a path to a text file
            transcript_alias (Union[str, Dict], optional): Dictionary mapping transcript aliases to their normalized ID, or a path to a text file
        """
        self.df = read_gtf(gtf, result_type="pandas")  # TODO: switch to 'polars'?
        self.cds_fasta = [read_fasta(i) for i in cds]
        self.dna_fasta = [read_fasta(i) for i in dna]
        self.protein_fasta = [read_fasta(i) for i in peptide]
        self.rna_fasta = [read_fasta(i) for i in rna]

        if isinstance(canonical_transcript, str):
            self._canonical_transcript = txt_to_list(canonical_transcript)
        else:
            self._canonical_transcript = canonical_transcript

        if isinstance(contig_alias, str):
            self._contig_alias = tsv_to_dict(contig_alias)
        else:
            self._contig_alias = contig_alias

        if isinstance(exon_alias, str):
            self._exon_alias = tsv_to_dict(exon_alias)
        else:
            self._exon_alias = exon_alias

        if isinstance(gene_alias, str):
            self._gene_alias = tsv_to_dict(gene_alias)
        else:
            self._gene_alias = gene_alias

        if isinstance(protein_alias, str):
            self._protein_alias = tsv_to_dict(protein_alias)
        else:
            self._protein_alias = protein_alias

        if isinstance(transcript_alias, str):
            self._transcript_alias = tsv_to_dict(transcript_alias)
        else:
            self._transcript_alias = transcript_alias

    # ---------------------------------------------------------------------------------------------
    # Functions for loading variants
    # ---------------------------------------------------------------------------------------------
    def parse(self, string: str) -> List:
        """Parse a variant string into a variant object.

        Args:
            string (str): String representing a variant, in HGVS format

        Returns:
            List: One or more normalized variants
        """
        parsed = parse(string)
        result = self.variant(
            position_type=parsed["breakpoint1"]["position_type"],
            feature=parsed["breakpoint1"]["feature"],
            start=parsed["breakpoint1"]["start"],
            start_offset=parsed["breakpoint1"]["start_offset"],
            end=parsed["breakpoint1"]["end"],
            end_offset=parsed["breakpoint1"]["end_offset"],
            strand=parsed["breakpoint1"]["strand"],
            refseq=parsed["breakpoint1"]["refseq"],
            altseq=parsed["breakpoint1"]["altseq"],
            variant_type=parsed["breakpoint1"]["variant_type"],
            position_type2=parsed["breakpoint2"]["position_type"],
            feature2=parsed["breakpoint2"]["feature"],
            start2=parsed["breakpoint2"]["start"],
            start_offset2=parsed["breakpoint2"]["start_offset"],
            end2=parsed["breakpoint2"]["end"],
            end_offset2=parsed["breakpoint2"]["end_offset"],
            strand2=parsed["breakpoint2"]["strand"],
            refseq2=parsed["breakpoint2"]["refseq"],
            altseq2=parsed["breakpoint2"]["altseq"],
        )

        return result

    # TODO: add type hints
    def variant(
        self,
        *,
        position_type: str,
        feature: str,
        start: int,
        start_offset: Optional[int] = None,
        end: Optional[int] = None,
        end_offset: Optional[int] = None,
        strand: Optional[str] = None,
        refseq: Optional[str] = None,
        altseq: Optional[str] = None,
        variant_type: Optional[str] = None,
        position_type2: Optional[str] = None,
        feature2: Optional[str] = None,
        start2: Optional[int] = None,
        start_offset2: Optional[int] = None,
        end2: Optional[int] = None,
        end_offset2: Optional[int] = None,
        strand2: Optional[str] = None,
        refseq2: Optional[str] = None,
        altseq2: Optional[str] = None,
    ) -> List:
        """Initialize one or more position or variants objects from the given arguments.

        Args:
            position_type (str): Position type for `start` and `end`. One of 'cdna', 'dna', 'exon', 'protein', or 'rna'.
            feature (str): Feature such as a transcript ID or gene name.
            start (int): Start position.
            start_offset (int, optional): Offset from `start`.
            end (int, optional): End position. Defaults to `start`.
            end_offset (int, optional): Offset from `end`.
            strand (str, optional): Strand the feature is on. One of '+' or '-'.
            refseq (str, optional): Reference allele. Required if the given position represents a variant.
            altseq (str, optional): Alternate allele. Required if the given position represents a variant.
            variant_type (str, optional): Variant type. Required if `refseq` or `altseq` is given.
            position_type2 (str, optional): For fusions. Position type for `start2` and `end2`. One of 'cdna', 'dna', 'exon', 'protein', or 'rna'.
            feature2 (str, optional): For fusions. Feature such as a transcript ID or gene name.
            start2 (int, optional): For fusions. Start position.
            start_offset2 (int, optional): For fusions. Offset from `start2`.
            end2 (int, optional): For fusions. End position. Defaults to `start2`.
            end_offset2 (int, optional): For fusions. Offset from `end2`.
            strand2 (str, optional): For fusions. Strand the feature is on. One of '+' or '-'.
            refseq2 (str, optional): For fusions. Reference allele. Required if the given position represents a variant.
            altseq2 (str, optional): For fusions. Alternate allele. Required if the given position represents a variant.

        Raises:
            ValueError: Arguments are missing.

        Returns:
            List: Position or Variants objects
        """
        result: List[Any] = []

        # Set defaults for missing inputs
        start_offset = cast(int, start_offset or 0)
        end = cast(int, end or start)
        end_offset = cast(int, end_offset or 0)
        strand = cast(str, strand or "")
        refseq = cast(str, refseq or "")
        altseq = cast(str, altseq or "")
        variant_type = cast(str, variant_type or "")

        if (refseq or altseq) and not variant_type:
            raise ValueError("refseq and/or altseq given without a variant_type")

        strand_ = [strand] if strand else ["+", "-"]

        # Load a fusion
        # TODO: infer is a fusion if position_type2, etc is given?
        if variant_type == FUSION:
            fusion: Optional[Type[_Fusion]] = None

            # TODO: default to any exon for exon fusions?

            if position_type2:
                position_type2 = cast(str, position_type2)
            else:
                raise ValueError(f"missing argument required for {FUSION}: 'position_type2'")

            if feature2:
                feature2 = cast(str, feature2)
            else:
                raise ValueError(f"missing argument required for {FUSION}: 'feature2'")

            if start2:
                start2 = cast(int, start2)
            else:
                raise ValueError(f"missing argument required for {FUSION}: 'start2'")

            # TODO: support mixed type fusions?
            if position_type == CDNA and position_type2 == CDNA:
                fusion = CdnaFusion
            elif position_type == DNA and position_type2 == DNA:
                fusion = DnaFusion
            elif position_type == EXON and position_type2 == EXON:
                fusion = ExonFusion
            elif position_type == PROTEIN and position_type2 == PROTEIN:
                fusion = ProteinFusion
            elif position_type == RNA and position_type2 == RNA:
                fusion = RnaFusion

            if fusion:
                breakpoint1 = self.variant(
                    feature=feature,
                    start=start,
                    position_type=position_type,
                    start_offset=start_offset,
                    end=end,
                    end_offset=end_offset,
                    strand=strand,
                    refseq=refseq,
                    altseq=altseq,
                )
                breakpoint2 = self.variant(
                    feature=feature2,
                    start=start2,
                    position_type=position_type2,
                    start_offset=start_offset2,
                    end=end2,
                    end_offset=end_offset2,
                    strand=strand2,
                    refseq=refseq2,
                    altseq=altseq2,
                )

                return list(fusion(b1, b2) for b1, b2 in product(breakpoint1, breakpoint2))

        # Load a small variant
        elif variant_type or refseq or altseq:
            if position_type == CDNA:
                transcript_ids = self.transcript_ids(feature)
                result = self._cdna_to_cdna_variant(
                    transcript_ids, start, start_offset, end, end_offset, strand_, refseq, altseq
                )
            elif position_type == DNA:
                contig_ids = self.contig_ids(feature)
                result = self._dna_to_dna_variant(
                    contig_ids, start, start_offset, end, end_offset, strand_, refseq, altseq
                )
            elif position_type == PROTEIN:
                transcript_ids = self.transcript_ids(feature)
                result = self._protein_to_protein_variant(
                    transcript_ids, start, start_offset, end, end_offset, strand_, refseq, altseq
                )
            elif position_type == RNA:
                transcript_ids = self.transcript_ids(feature)
                result = self._rna_to_rna_variant(
                    transcript_ids, start, start_offset, end, end_offset, strand_, refseq, altseq
                )

        # Load a position
        else:
            if position_type == CDNA:
                transcript_ids = self.transcript_ids(feature)
                result = self._cdna_to_cdna(
                    transcript_ids, start, start_offset, end, end_offset, strand_
                )
            elif position_type == DNA:
                contig_ids = self.contig_ids(feature)
                result = self._dna_to_dna(contig_ids, start, start_offset, end, end_offset, strand_)
            elif position_type == EXON:
                transcript_ids = self.transcript_ids(feature)
                result = self._exon_to_exon(
                    transcript_ids, start, start_offset, end, end_offset, strand_
                )
            elif position_type == PROTEIN:
                transcript_ids = self.transcript_ids(feature)
                result = self._protein_to_protein(
                    transcript_ids, start, start_offset, end, end_offset, strand_
                )
            elif position_type == RNA:
                transcript_ids = self.transcript_ids(feature)
                result = self._rna_to_rna(
                    transcript_ids, start, start_offset, end, end_offset, strand_
                )

        return result

    # ---------------------------------------------------------------------------------------------
    # Functions for getting variant sequences
    # ---------------------------------------------------------------------------------------------
    # TODO: add type hints
    def sequence(self, position) -> str:
        """Return the sequence for the given position, inclusive.

        Args:
            position (Position): Position to retrieve sequence for

        Returns:
            str: Sequence

        Raises:
            ValueError: No method exists for getting a sequence for the given position type
        """
        # TODO: Get sequence for offset variants?
        if position.start_offset or position.end_offset:
            raise ValueError(f"Unable to get sequence for offset position {position}")

        if position.is_cdna:
            position = cast(CdnaPosition, position)
            return self._cds_sequence(position.transcript_id, position.start, position.end)
        elif position.is_dna:
            position = cast(DnaPosition, position)
            return self._dna_sequence(
                position.contig_id, position.start, position.end, position.strand
            )
        elif position.is_exon:
            raise NotImplementedError(f"Unable to get sequence for {position}")
        elif position.is_protein:
            position = cast(ProteinPosition, position)
            return self._protein_sequence(position.protein_id, position.start, position.end)
        elif position.is_rna:
            position = cast(RnaPosition, position)
            return self._rna_sequence(position.transcript_id, position.start, position.end)
        else:
            raise ValueError(f"Unable to get sequence for {position}")

    def _cds_sequence(self, transcript_id: str, start: int, end: int) -> str:
        """Return the sequence for the given position, inclusive.

        Args:
            transcript_id (str): Transcript ID
            start (int): Start position
            end (int): End position

        Returns:
            str: CDS sequence
        """
        return self._sequence(self.cds_fasta, transcript_id, start=start, end=end)

    def _dna_sequence(self, contig_id: str, start: int, end: int, strand: str) -> str:
        """Return the sequence for the given position, inclusive.

        Args:
            contig_id (str): Contig ID
            start (int): Start position
            end (int): End position
            strand (str): Strand ('+' or '-')

        Returns:
            str: DNA sequence
        """
        return self._sequence(self.dna_fasta, contig_id, start=start, end=end, strand=strand)

    def _protein_sequence(self, protein_id: str, start: int, end: int) -> str:
        """Return the sequence for the given position, inclusive.

        Args:
            protein_id (str): Protein ID
            start (int): Start position
            end (int): End position

        Returns:
            str: Protein sequence
        """
        return self._sequence(self.protein_fasta, protein_id, start=start, end=end)

    def _rna_sequence(self, transcript_id: str, start: int, end: int) -> str:
        """Return the sequence for the given position, inclusive.

        Args:
            transcript_id (str): Transcript ID
            start (int): Start position
            end (int): End position

        Returns:
            str: RNA sequence
        """
        return self._sequence(self.rna_fasta, transcript_id, start=start, end=end)

    def _sequence(
        self,
        fasta: List[Fasta],
        ref: str,
        start: Optional[int] = None,
        end: Optional[int] = None,
        strand: Optional[str] = None,
    ):
        for f in fasta:
            try:
                seq = f[ref]
                break
            except KeyError:
                continue
        else:
            raise KeyError(f"Sequence '{ref}' not found")

        # if no positions are given, return the whole sequence
        seqlen = len(seq)
        start = start if start is not None else 1
        end = end if end is not None else seqlen

        # validate that the given positions fall within the sequence
        if not (0 <= start <= seqlen):
            raise ValueError(f"Start must be from 1 to {seqlen} ({start})")
        if not (0 < end <= seqlen):
            raise ValueError(f"End must be from 1 to {seqlen} ({end})")

        # sanity check that the end position is after the start
        if end < start:
            raise ValueError(f"End must be >= start ({end} < {start})")

        subseq = seq[start - 1 : end]
        if strand == "-":
            subseq = reverse_complement(subseq)

        return subseq

    # ---------------------------------------------------------------------------------------------
    # Functions for mapping to other position types
    # ---------------------------------------------------------------------------------------------
    # TODO: add type hints
    def to_cdna(self, position) -> List:
        """Map a position to zero or more cDNA positions.

        Args:
            position (Position): Position or variant object.

        Returns:
            List[CdnaPosition]: If a variant was given, returns a variant, otherwise returns a position.
        """
        if isinstance(position, str):
            return self._string_to_cdna(position)
        else:
            return self._position_to_cdna(position)

    def to_dna(self, position) -> List:
        """Map a position to zero or more DNA positions.

        Args:
            position (Position): Position or variant object.

        Returns:
            List[DnaPosition]: If a variant was given, returns a variant, otherwise returns a position.
        """
        if isinstance(position, str):
            return self._string_to_dna(position)
        else:
            return self._position_to_dna(position)

    def to_exon(self, position) -> List:
        """Map a position to zero or more exon positions.

        Args:
            position (Position): Position or variant object.

        Returns:
            List[ExonPosition]: If a variant was given, returns a variant, otherwise returns a position.
        """
        if isinstance(position, str):
            return self._string_to_exon(position)
        else:
            return self._position_to_exon(position)

    def to_protein(self, position) -> List:
        """Map a position to zero or more protein positions.

        Args:
            position (Position): Position or variant object.

        Returns:
            List[ProteinPosition]: If a variant was given, returns a variant, otherwise returns a position.
        """
        if isinstance(position, str):
            return self._string_to_protein(position)
        else:
            return self._position_to_protein(position)

    def to_rna(self, position) -> List:
        """Map a position to zero or more RNA positions.

        Args:
            position (Position): Position or variant object.

        Returns:
            List[RnaPosition]: If a variant was given, returns a variant, otherwise returns a position.
        """
        if isinstance(position, str):
            return self._string_to_rna(position)
        else:
            return self._position_to_rna(position)

    def _string_to_cdna(self, string: str):
        return self._string_to(string, self._position_to_cdna)

    def _string_to_dna(self, string: str):
        return self._string_to(string, self._position_to_dna)

    def _string_to_exon(self, string: str):
        return self._string_to(string, self._position_to_exon)

    def _string_to_protein(self, string: str):
        return self._string_to(string, self._position_to_protein)

    def _string_to_rna(self, string: str):
        return self._string_to(string, self._position_to_rna)

    def _string_to(self, string: str, func: Callable):
        result = []

        for position in self.parse(string):
            result.extend(func(position))

        return result

    def _position_to_cdna(self, position):
        return self._position_to(
            position,
            fusionf=self._position_to_cdna,
            fusiont=CdnaFusion,
            cdnavf=self._cdna_to_cdna_variant,
            dnavf=self._dna_to_cdna_variant,
            proteinvf=self._protein_to_cdna_variant,
            rnavf=self._rna_to_cdna_variant,
            cdnaf=self._cdna_to_cdna,
            dnaf=self._dna_to_cdna,
            exonf=self._exon_to_cdna,
            proteinf=self._protein_to_cdna,
            rnaf=self._rna_to_cdna,
        )

    def _position_to_dna(self, position):
        return self._position_to(
            position,
            fusionf=self._position_to_dna,
            fusiont=DnaFusion,
            cdnavf=self._cdna_to_dna_variant,
            dnavf=self._dna_to_dna_variant,
            proteinvf=self._protein_to_dna_variant,
            rnavf=self._rna_to_dna_variant,
            cdnaf=self._cdna_to_dna,
            dnaf=self._dna_to_dna,
            exonf=self._exon_to_dna,
            proteinf=self._protein_to_dna,
            rnaf=self._rna_to_dna,
        )

    def _position_to_exon(self, position):
        return self._position_to(
            position,
            fusionf=self._position_to_exon,
            fusiont=ExonFusion,
            cdnavf=self._cdna_to_exon_variant,
            dnavf=self._dna_to_exon_variant,
            proteinvf=self._protein_to_exon_variant,
            rnavf=self._rna_to_exon_variant,
            cdnaf=self._cdna_to_exon,
            dnaf=self._dna_to_exon,
            exonf=self._exon_to_exon,
            proteinf=self._protein_to_exon,
            rnaf=self._rna_to_exon,
        )

    def _position_to_protein(self, position):
        return self._position_to(
            position,
            fusionf=self._position_to_protein,
            fusiont=ProteinFusion,
            cdnavf=self._cdna_to_protein_variant,
            dnavf=self._dna_to_protein_variant,
            proteinvf=self._protein_to_protein_variant,
            rnavf=self._rna_to_protein_variant,
            cdnaf=self._cdna_to_protein,
            dnaf=self._dna_to_protein,
            exonf=self._exon_to_protein,
            proteinf=self._protein_to_protein,
            rnaf=self._rna_to_protein,
        )

    def _position_to_rna(self, position):
        return self._position_to(
            position,
            fusionf=self._position_to_rna,
            fusiont=RnaFusion,
            cdnavf=self._cdna_to_rna_variant,
            dnavf=self._dna_to_rna_variant,
            proteinvf=self._protein_to_rna_variant,
            rnavf=self._rna_to_rna_variant,
            cdnaf=self._cdna_to_rna,
            dnaf=self._dna_to_rna,
            exonf=self._exon_to_rna,
            proteinf=self._protein_to_rna,
            rnaf=self._rna_to_rna,
        )

    def _position_to(
        self,
        position,
        fusionf: Callable,
        fusiont: Type,
        cdnavf: Callable,
        dnavf: Callable,
        proteinvf: Callable,
        rnavf: Callable,
        cdnaf: Callable,
        dnaf: Callable,
        exonf: Callable,
        proteinf: Callable,
        rnaf: Callable,
    ) -> List:
        if position.is_fusion:
            fusion = cast(_Fusion, position)
            breakpoint1 = fusionf(fusion.breakpoint1)
            breakpoint2 = fusionf(fusion.breakpoint2)
            return [fusiont(b1, b2) for b1, b2 in product(breakpoint1, breakpoint2)]
        elif position.is_small_variant:
            if position.is_cdna:
                position = cast(_CdnaSmallVariant, position)
                return cdnavf(
                    [position.transcript_id],
                    position.start,
                    position.start_offset,
                    position.end,
                    position.end_offset,
                    [position.strand],
                    position.refseq,
                    position.altseq,
                )
            elif position.is_dna:
                position = cast(_DnaSmallVariant, position)
                return dnavf(
                    [position.contig_id],
                    position.start,
                    position.start_offset,
                    position.end,
                    position.end_offset,
                    [position.strand],
                    position.refseq,
                    position.altseq,
                )
            elif position.is_protein:
                position = cast(_ProteinSmallVariant, position)
                return proteinvf(
                    [position.transcript_id],
                    position.start,
                    position.start_offset,
                    position.end,
                    position.end_offset,
                    [position.strand],
                    position.refseq,
                    position.altseq,
                )
            elif position.is_rna:
                position = cast(_RnaSmallVariant, position)
                return rnavf(
                    [position.transcript_id],
                    position.start,
                    position.start_offset,
                    position.end,
                    position.end_offset,
                    [position.strand],
                    position.refseq,
                    position.altseq,
                )
        else:
            if position.is_cdna:
                position = cast(CdnaPosition, position)
                return cdnaf(
                    [position.transcript_id],
                    position.start,
                    position.start_offset,
                    position.end,
                    position.end_offset,
                    [position.strand],
                )
            elif position.is_dna:
                position = cast(DnaPosition, position)
                return dnaf(
                    [position.contig_id],
                    position.start,
                    position.start_offset,
                    position.end,
                    position.end_offset,
                    [position.strand],
                )
            elif position.is_exon:
                position = cast(ExonPosition, position)
                return exonf(
                    [position.transcript_id],
                    position.start,
                    position.start_offset,
                    position.end,
                    position.end_offset,
                    [position.strand],
                )
            elif position.is_protein:
                position = cast(ProteinPosition, position)
                return proteinf(
                    [position.transcript_id],
                    position.start,
                    position.start_offset,
                    position.end,
                    position.end_offset,
                    [position.strand],
                )
            elif position.is_rna:
                position = cast(RnaPosition, position)
                return rnaf(
                    [position.transcript_id],
                    position.start,
                    position.start_offset,
                    position.end,
                    position.end_offset,
                    [position.strand],
                )

        raise AssertionError(f"Unknown position type for {position}")

    def _cdna_to_cdna(
        self,
        transcript_id: List[str],
        start: int,
        start_offset: int,
        end: int,
        end_offset: int,
        strand: List[str],
        include_stop: bool = True,
    ) -> List[CdnaPosition]:
        feature = [CDS, STOP_CODON] if include_stop else [CDS]

        def convert(n: int, offset: int):
            result = []

            # If the offset position can be normalized to a non-offset position, do so. Otherwise
            # just return an offset position.
            if offset:
                for dna in self._cdna_to_dna(
                    transcript_id,
                    start,
                    start_offset,
                    end,
                    end_offset,
                    strand,
                    include_stop=include_stop,
                ):
                    for cdna in self._dna_to_cdna(
                        [dna.contig_id],
                        dna.start,
                        dna.start_offset,
                        dna.end,
                        dna.end_offset,
                        [dna.strand],
                    ):
                        if cdna.transcript_id in transcript_id:
                            result.append(cdna)

                if result:
                    return result

            mask = (
                (self.df[TRANSCRIPT_ID].isin(transcript_id))
                & (self.df["cdna_start"] <= n)
                & (self.df["cdna_end"] >= n)
                & (self.df["strand"].isin(strand))
                & (self.df["feature"].isin(feature))
            )
            for _, cds in self.df[mask].iterrows():
                result.append(
                    CdnaPosition(
                        contig_id=cds.contig_id,
                        start=start,
                        start_offset=offset,
                        end=end,
                        end_offset=offset,
                        strand=cds.strand,
                        gene_id=cds.gene_id,
                        gene_name=cds.gene_name,
                        transcript_id=cds.transcript_id,
                        transcript_name=cds.transcript_name,
                        protein_id=cds.protein_id,
                    )
                )

            return result

        result_start = convert(start, start_offset)
        if start == end:
            return sorted(result_start)
        else:
            result_end = convert(end, end_offset)
            return join_positions(result_start, result_end, TRANSCRIPT_ID)

    def _cdna_to_cdna_variant(
        self,
        transcript_id: List[str],
        start: int,
        start_offset: int,
        end: int,
        end_offset: int,
        strand: List[str],
        refseq: str,
        altseq: str,
        include_stop: bool = True,
    ) -> List[_CdnaSmallVariant]:
        result = []

        for cdna in self._cdna_to_cdna(
            transcript_id, start, start_offset, end, end_offset, strand, include_stop=include_stop
        ):
            result.extend(self._variant_class_from_cdna(cdna, refseq, altseq))

        return result

    def _cdna_to_dna(
        self,
        transcript_id: List[str],
        start: int,
        start_offset: int,
        end: int,
        end_offset: int,
        strand: List[str],
        include_stop: bool = True,
    ) -> List[DnaPosition]:
        feature = [CDS, STOP_CODON] if include_stop else [CDS]

        def convert(n: int, offset: int):
            result = []

            mask = (
                (self.df[TRANSCRIPT_ID].isin(transcript_id))
                & (self.df["cdna_start"] <= n)
                & (self.df["cdna_end"] >= n)
                & (self.df["strand"].isin(strand))
                & (self.df["feature"].isin(feature))
            )
            for _, cds in self.df[mask].iterrows():
                if cds.strand == "-":
                    new_start = new_end = cds.end - (n - cds.cdna_start) - offset
                else:
                    new_start = new_end = cds.start + (n - cds.cdna_start) + offset

                # TODO: Check that new new_start is actually on the contig
                result.append(
                    DnaPosition(
                        contig_id=cds.contig_id,
                        start=new_start,
                        start_offset=0,
                        end=new_end,
                        end_offset=0,
                        strand=cds.strand,
                    )
                )

            return result

        result_start = convert(start, start_offset)
        if start == end:
            return sorted(result_start)
        else:
            result_end = convert(end, end_offset)
            return join_positions(result_start, result_end, CONTIG_ID)

    def _cdna_to_dna_variant(
        self,
        transcript_id: List[str],
        start: int,
        start_offset: int,
        end: int,
        end_offset: int,
        strand: List[str],
        refseq: str,
        altseq: str,
        include_stop: bool = True,
    ) -> List[_DnaSmallVariant]:
        result = []

        for dna in self._cdna_to_dna(
            transcript_id, start, start_offset, end, end_offset, strand, include_stop=include_stop
        ):
            result.extend(self._dna_small_variant_from_dna(dna, refseq, altseq))

        return result

    def _cdna_to_exon(
        self,
        transcript_id: List[str],
        start: int,
        start_offset: int,
        end: int,
        end_offset: int,
        strand: List[str],
        include_stop: bool = True,
    ) -> List[ExonPosition]:
        feature = [CDS, STOP_CODON] if include_stop else [CDS]

        def convert(n: int, offset: int):
            result = []

            # For an offset position, we need to calculate then equivalent DNA position then map
            # that to one or more exons. This is slower, so if there's no offset we can just map
            # directly to an exon.
            if offset:
                for dna in self._cdna_to_dna(
                    transcript_id,
                    start,
                    start_offset,
                    end,
                    end_offset,
                    strand,
                    include_stop=include_stop,
                ):
                    for exon in self._dna_to_exon(
                        [dna.contig_id],
                        dna.start,
                        dna.start_offset,
                        dna.end,
                        dna.end_offset,
                        [dna.strand],
                    ):
                        if exon.transcript_id in transcript_id:
                            result.append(exon)

                return result

            mask_cds = (
                (self.df[TRANSCRIPT_ID].isin(transcript_id))
                & (self.df["cdna_start"] <= n)
                & (self.df["cdna_end"] >= n)
                & (self.df["strand"].isin(strand))
                & (self.df["feature"].isin(feature))
            )
            for _, cds in self.df[mask_cds].iterrows():
                mask_exon = (
                    (self.df[TRANSCRIPT_ID].isin(transcript_id))
                    & (self.df["exon_number"] == cds.exon_number)
                    & (self.df["feature"] == EXON)
                )
                for _, exon_row in self.df[mask_exon].iterrows():
                    result.append(
                        ExonPosition(
                            contig_id=exon_row.contig_id,
                            start=int(exon_row.exon_number),
                            start_offset=0,
                            end=int(exon_row.exon_number),
                            end_offset=0,
                            strand=exon_row.strand,
                            gene_id=exon_row.gene_id,
                            gene_name=exon_row.gene_name,
                            transcript_id=exon_row.transcript_id,
                            transcript_name=exon_row.transcript_name,
                            exon_id=exon_row.exon_id,
                        )
                    )

            return result

        result_start = convert(start, start_offset)
        if start == end:
            return sorted(result_start)
        else:
            result_end = convert(end, end_offset)
            return join_positions(result_start, result_end, TRANSCRIPT_ID)

    def _cdna_to_exon_variant(
        self,
        transcript_id: List[str],
        start: int,
        start_offset: int,
        end: int,
        end_offset: int,
        strand: List[str],
        refseq: str,
        altseq: str,
        include_stop: bool = True,
    ) -> List[_ExonSmallVariant]:
        result = []

        for exon in self._cdna_to_exon(
            transcript_id, start, start_offset, end, end_offset, strand, include_stop=include_stop
        ):
            result.extend(self._exon_small_variant_from_exon(exon, refseq, altseq))

        return result

    def _cdna_to_protein(
        self,
        transcript_id: List[str],
        start: int,
        start_offset: int,
        end: int,
        end_offset: int,
        strand: List[str],
    ) -> List[ProteinPosition]:
        result = []

        for cdna in self._cdna_to_cdna(
            transcript_id, start, start_offset, end, end_offset, strand, include_stop=False
        ):
            # If the postion wasn't mapped to a non-offset position by `_cdna_to_cdna`, it means
            # that the position does not map to a protein.
            if cdna.start_offset or cdna.end_offset:
                continue

            # Convert the cDNA position to a protein position
            protein_start = calc_cdna_to_protein(cdna.start)
            protein_end = calc_cdna_to_protein(cdna.end)
            result.append(
                ProteinPosition.copy_from(
                    cdna, start=protein_start, start_offset=0, end=protein_end, end_offset=0
                )
            )

        return result

    def _cdna_to_protein_variant(
        self,
        transcript_id: List[str],
        start: int,
        start_offset: int,
        end: int,
        end_offset: int,
        strand: List[str],
        refseq: str,
        altseq: str,
    ) -> List[_ProteinSmallVariant]:
        # TODO: A lot of this function is duplicated from _cdna_to_protein()
        result = []

        for cdna in self._cdna_to_cdna(
            transcript_id, start, start_offset, end, end_offset, strand, include_stop=False
        ):
            # If the postion wasn't mapped to a non-offset position by `_cdna_to_cdna`, it means
            # that the position does not map to a protein.
            if cdna.start_offset or cdna.end_offset:
                continue

            # Convert the cDNA position to a protein position
            protein_start = calc_cdna_to_protein(cdna.start)
            protein_end = calc_cdna_to_protein(cdna.end)
            protein = ProteinPosition.copy_from(
                cdna, start=protein_start, start_offset=0, end=protein_end, end_offset=0
            )
            result.extend(self._protein_small_variant(cdna, protein, refseq, altseq))

        return result

    def _cdna_to_rna(
        self,
        transcript_id: List[str],
        start: int,
        start_offset: int,
        end: int,
        end_offset: int,
        strand: List[str],
        include_stop: bool = True,
    ) -> List[RnaPosition]:
        feature = [CDS, STOP_CODON] if include_stop else [CDS]

        def convert(n: int, offset: int):
            result = []

            # If the offset position can be normalized to a non-offset position, do so. Otherwise
            # just return an offset
            if offset:
                for dna in self._cdna_to_dna(
                    transcript_id,
                    start,
                    start_offset,
                    end,
                    end_offset,
                    strand,
                    include_stop=include_stop,
                ):
                    for rna in self._dna_to_rna(
                        [dna.contig_id],
                        dna.start,
                        dna.start_offset,
                        dna.end,
                        dna.end_offset,
                        [dna.strand],
                    ):
                        if rna.transcript_id in transcript_id:
                            result.append(rna)

                if result:
                    return result

            mask = (
                (self.df[TRANSCRIPT_ID].isin(transcript_id))
                & (self.df["cdna_start"] <= n)
                & (self.df["cdna_end"] >= n)
                & (self.df["strand"].isin(strand))
                & (self.df["feature"].isin(feature))
            )
            for _, cds in self.df[mask].iterrows():
                new_start = new_end = cds.transcript_start + (n - cds.cdna_start)
                result.append(
                    RnaPosition(
                        contig_id=cds.contig_id,
                        start=new_start,
                        start_offset=offset,
                        end=new_end,
                        end_offset=offset,
                        strand=cds.strand,
                        gene_id=cds.gene_id,
                        gene_name=cds.gene_name,
                        transcript_id=cds.transcript_id,
                        transcript_name=cds.transcript_name,
                    )
                )

            return result

        result_start = convert(start, start_offset)
        if start == end:
            return sorted(result_start)
        else:
            result_end = convert(end, end_offset)
            return join_positions(result_start, result_end, TRANSCRIPT_ID)

    def _cdna_to_rna_variant(
        self,
        transcript_id: List[str],
        start: int,
        start_offset: int,
        end: int,
        end_offset: int,
        strand: List[str],
        refseq: str,
        altseq: str,
        include_stop: bool = True,
    ) -> List[_RnaSmallVariant]:
        result = []

        for rna in self._cdna_to_rna(
            transcript_id, start, start_offset, end, end_offset, strand, include_stop=include_stop
        ):
            result.extend(self._rna_small_variant_from_rna(rna, refseq, altseq))

        return result

    def _dna_to_cdna(
        self,
        contig_id: List[str],
        start: int,
        start_offset: int,
        end: int,
        end_offset: int,
        strand: List[str],
        include_stop: bool = True,
    ) -> List[CdnaPosition]:
        feature = [CDS, STOP_CODON] if include_stop else [CDS]

        def convert(n: int, offset: int):
            result = []

            for strand_ in strand:
                if strand_ == "-":
                    n_ = n - offset
                else:
                    n_ = n + offset

                mask = (
                    (self.df[CONTIG_ID].isin(contig_id))
                    & (self.df["start"] <= n_)
                    & (self.df["end"] >= n_)
                    & (self.df["strand"] == strand_)
                    & (self.df["feature"].isin(feature))
                )
                for _, cds in self.df[mask].iterrows():
                    if cds.strand == "-":
                        new_start = new_end = cds.end - n_ + cds.cdna_start
                    else:
                        new_start = new_end = n_ - cds.start + cds.cdna_start

                    result.append(
                        CdnaPosition(
                            contig_id=cds.contig_id,
                            start=new_start,
                            start_offset=0,
                            end=new_end,
                            end_offset=0,
                            strand=cds.strand,
                            gene_id=cds.gene_id,
                            gene_name=cds.gene_name,
                            transcript_id=cds.transcript_id,
                            transcript_name=cds.transcript_name,
                            protein_id=cds.protein_id,
                        )
                    )

            return result

        result_start = convert(start, start_offset)
        if start == end:
            return sorted(result_start)
        else:
            result_end = convert(end, end_offset)
            return join_positions(result_start, result_end, TRANSCRIPT_ID)

    def _dna_to_cdna_variant(
        self,
        contig_id: List[str],
        start: int,
        start_offset: int,
        end: int,
        end_offset: int,
        strand: List[str],
        refseq: str,
        altseq: str,
        include_stop: bool = True,
    ) -> List[_CdnaSmallVariant]:
        result = []

        for cdna in self._dna_to_cdna(
            contig_id, start, start_offset, end, end_offset, strand, include_stop=include_stop
        ):
            result.extend(self._variant_class_from_cdna(cdna, refseq, altseq))

        return result

    def _dna_to_dna(
        self,
        contig_id: List[str],
        start: int,
        start_offset: int,
        end: int,
        end_offset: int,
        strand: List[str],
    ) -> List[DnaPosition]:
        result = []

        # TODO: Check that new new_start is actually on the contig
        for contig_id_, strand_ in product(contig_id, strand):
            if strand_ == "-":
                new_start = start - start_offset
                new_end = end - end_offset
            else:
                new_start = start + start_offset
                new_end = end + end_offset

            # Sort the start and end positions after adjusting by offsets
            new_start, new_end = sorted([new_start, new_end])

            result.append(
                DnaPosition(
                    contig_id=contig_id_,
                    start=new_start,
                    start_offset=0,
                    end=new_end,
                    end_offset=0,
                    strand=strand_,
                )
            )

        return result

    def _dna_to_dna_variant(
        self,
        contig_id: List[str],
        start: int,
        start_offset: int,
        end: int,
        end_offset: int,
        strand: List[str],
        refseq: str,
        altseq: str,
    ) -> List[_DnaSmallVariant]:
        result = []

        for dna in self._dna_to_dna(contig_id, start, start_offset, end, end_offset, strand):
            result.extend(self._dna_small_variant_from_dna(dna, refseq, altseq))

        return result

    def _dna_to_exon(
        self,
        contig_id: List[str],
        start: int,
        start_offset: int,
        end: int,
        end_offset: int,
        strand: List[str],
    ) -> List[ExonPosition]:
        def convert(n: int, offset: int):
            result = []

            for strand_ in strand:
                if strand_ == "-":
                    n_ = n - offset
                else:
                    n_ = n + offset

                mask = (
                    (self.df[CONTIG_ID].isin(contig_id))
                    & (self.df["start"] <= n_)
                    & (self.df["end"] >= n_)
                    & (self.df["strand"] == strand_)
                    & (self.df["feature"] == EXON)
                )
                for _, exon in self.df[mask].iterrows():
                    result.append(
                        ExonPosition(
                            contig_id=exon.contig_id,
                            start=int(exon.exon_number),
                            start_offset=0,
                            end=int(exon.exon_number),
                            end_offset=0,
                            strand=exon.strand,
                            gene_id=exon.gene_id,
                            gene_name=exon.gene_name,
                            transcript_id=exon.transcript_id,
                            transcript_name=exon.transcript_name,
                            exon_id=exon.exon_id,
                        )
                    )

            return result

        result_start = convert(start, start_offset)
        if start == end:
            return sorted(result_start)
        else:
            result_end = convert(end, end_offset)
            return join_positions(result_start, result_end, TRANSCRIPT_ID)

    def _dna_to_exon_variant(
        self,
        contig_id: List[str],
        start: int,
        start_offset: int,
        end: int,
        end_offset: int,
        strand: List[str],
        refseq: str,
        altseq: str,
    ) -> List[_ExonSmallVariant]:
        result = []

        for exon in self._dna_to_exon(contig_id, start, start_offset, end, end_offset, strand):
            result.extend(self._exon_small_variant_from_exon(exon, refseq, altseq))

        return result

    def _dna_to_protein(
        self,
        contig_id: List[str],
        start: int,
        start_offset: int,
        end: int,
        end_offset: int,
        strand: List[str],
    ) -> List[ProteinPosition]:
        result = []

        for cdna in self._dna_to_cdna(
            contig_id, start, start_offset, end, end_offset, strand, include_stop=False
        ):
            # Offset cDNA position are assumed to not map to a protein
            if cdna.start_offset or cdna.end_offset:
                continue

            pstart = calc_cdna_to_protein(cdna.start)
            pend = calc_cdna_to_protein(cdna.end)
            result.append(
                ProteinPosition.copy_from(
                    cdna, start=pstart, start_offset=0, end=pend, end_offset=0
                )
            )

        return result

    def _dna_to_protein_variant(
        self,
        contig_id: List[str],
        start: int,
        start_offset: int,
        end: int,
        end_offset: int,
        strand: List[str],
        refseq: str,
        altseq: str,
    ) -> List[_ProteinSmallVariant]:
        # TODO: A lot of this function is duplicated from _dna_to_protein()
        result = []

        for cdna in self._dna_to_cdna(
            contig_id, start, start_offset, end, end_offset, strand, include_stop=False
        ):
            # Offset cDNA position are assumed to not map to a protein
            if cdna.start_offset or cdna.end_offset:
                continue

            pstart = calc_cdna_to_protein(cdna.start)
            pend = calc_cdna_to_protein(cdna.end)
            protein = ProteinPosition.copy_from(
                cdna, start=pstart, start_offset=0, end=pend, end_offset=0
            )
            result.extend(self._protein_small_variant(cdna, protein, refseq, altseq))

        return result

    def _dna_to_rna(
        self,
        contig_id: List[str],
        start: int,
        start_offset: int,
        end: int,
        end_offset: int,
        strand: List[str],
    ) -> List[RnaPosition]:
        def convert(n: int, offset: int):
            result = []

            for strand_ in strand:
                if strand_ == "-":
                    n_ = n - offset
                else:
                    n_ = n + offset

                mask = (
                    (self.df[CONTIG_ID].isin(contig_id))
                    & (self.df["start"] <= n_)
                    & (self.df["end"] >= n_)
                    & (self.df["strand"] == strand_)
                    & (self.df["feature"] == EXON)
                )
                for _, exon in self.df[mask].iterrows():
                    if exon.strand == "-":
                        new_start = new_end = exon.end - n_ + exon.transcript_start
                    else:
                        new_start = new_end = n_ - exon.start + exon.transcript_start

                    result.append(
                        RnaPosition(
                            contig_id=exon.contig_id,
                            start=new_start,
                            start_offset=0,
                            end=new_end,
                            end_offset=0,
                            strand=exon.strand,
                            gene_id=exon.gene_id,
                            gene_name=exon.gene_name,
                            transcript_id=exon.transcript_id,
                            transcript_name=exon.transcript_name,
                        )
                    )

            return result

        result_start = convert(start, start_offset)
        if start == end:
            return sorted(result_start)
        else:
            result_end = convert(end, end_offset)
            return join_positions(result_start, result_end, TRANSCRIPT_ID)

    def _dna_to_rna_variant(
        self,
        contig_id: List[str],
        start: int,
        start_offset: int,
        end: int,
        end_offset: int,
        strand: List[str],
        refseq: str,
        altseq: str,
    ) -> List[_RnaSmallVariant]:
        result = []

        for rna in self._dna_to_rna(contig_id, start, start_offset, end, end_offset, strand):
            result.extend(self._rna_small_variant_from_rna(rna, refseq, altseq))

        return result

    def _exon_to_cdna(
        self,
        transcript_id: List[str],
        start: int,
        start_offset: int,
        end: int,
        end_offset: int,
        strand: List[str],
        include_stop: bool = True,
    ) -> List[CdnaPosition]:
        feature = [CDS, STOP_CODON] if include_stop else [CDS]

        def convert(n: int, offset: int):
            result = []

            # TODO: Is there a reasonable case where an exon position would have an offset?
            assert not offset, offset

            mask = (
                (self.df[TRANSCRIPT_ID].isin(transcript_id))
                & (self.df["exon_number"] == float(n))
                & (self.df["strand"].isin(strand))
                & (self.df["feature"].isin(feature))
            )
            for _, cds in self.df[mask].iterrows():
                result.append(
                    CdnaPosition(
                        contig_id=cds.contig_id,
                        start=cds.cdna_start,
                        start_offset=0,
                        end=cds.cdna_end,
                        end_offset=0,
                        strand=cds.strand,
                        gene_id=cds.gene_id,
                        gene_name=cds.gene_name,
                        transcript_id=cds.transcript_id,
                        transcript_name=cds.transcript_name,
                        protein_id=cds.protein_id,
                    )
                )

            return result

        result_start = convert(start, start_offset)
        if start == end:
            return sorted(result_start)
        else:
            result_end = convert(end, end_offset)
            return join_positions(result_start, result_end, TRANSCRIPT_ID)

    def _exon_to_dna(
        self,
        transcript_id: List[str],
        start: int,
        start_offset: int,
        end: int,
        end_offset: int,
        strand: List[str],
    ) -> List[DnaPosition]:
        def convert(n: int, offset: int):
            result = []

            # TODO: Is there a reasonable case where an exon position would have an offset?
            assert not offset, offset

            mask = (
                (self.df[TRANSCRIPT_ID].isin(transcript_id))
                & (self.df["exon_number"] == float(n))
                & (self.df["strand"].isin(strand))
                & (self.df["feature"] == EXON)
            )
            for _, exon in self.df[mask].iterrows():
                result.append(
                    DnaPosition(
                        contig_id=exon.contig_id,
                        start=exon.start,
                        start_offset=0,
                        end=exon.end,
                        end_offset=0,
                        strand=exon.strand,
                    )
                )

            return result

        result_start = convert(start, start_offset)
        if start == end:
            return sorted(result_start)
        else:
            result_end = convert(end, end_offset)
            return join_positions(result_start, result_end, CONTIG_ID)

    def _exon_to_exon(
        self,
        transcript_id: List[str],
        start: int,
        start_offset: int,
        end: int,
        end_offset: int,
        strand: List[str],
    ) -> List[ExonPosition]:
        def convert(n: int, offset):
            result = []

            # TODO: Is there a reasonable case where an exon position would have an offset?
            assert not offset, offset

            mask = (
                (self.df[TRANSCRIPT_ID].isin(transcript_id))
                & (self.df["exon_number"] == float(n))
                & (self.df["strand"].isin(strand))
                & (self.df["feature"] == EXON)
            )
            for _, exon in self.df[mask].iterrows():
                result.append(
                    ExonPosition(
                        contig_id=exon.contig_id,
                        start=int(exon.exon_number),
                        start_offset=0,
                        end=int(exon.exon_number),
                        end_offset=0,
                        strand=exon.strand,
                        gene_id=exon.gene_id,
                        gene_name=exon.gene_name,
                        transcript_id=exon.transcript_id,
                        transcript_name=exon.transcript_name,
                        exon_id=exon.exon_id,
                    )
                )

            return result

        result_start = convert(start, start_offset)
        if start == end:
            return sorted(result_start)
        else:
            result_end = convert(end, end_offset)
            return join_positions(result_start, result_end, TRANSCRIPT_ID)

    def _exon_to_protein(
        self,
        transcript_id: List[str],
        start: int,
        start_offset: int,
        end: int,
        end_offset: int,
        strand: List[str],
    ) -> List[ProteinPosition]:
        result = []

        for cdna in self._exon_to_cdna(
            transcript_id, start, start_offset, end, end_offset, strand, include_stop=False
        ):
            result.extend(
                self._cdna_to_protein(
                    [cdna.transcript_id],
                    cdna.start,
                    cdna.start_offset,
                    cdna.end,
                    cdna.end_offset,
                    [cdna.strand],
                )
            )

        return result

    def _exon_to_rna(
        self,
        transcript_id: List[str],
        start: int,
        start_offset: int,
        end: int,
        end_offset: int,
        strand: List[str],
    ) -> List[RnaPosition]:
        def convert(n: int, offset: int):
            result = []

            # TODO: Is there a reasonable case where an exon position would have an offset?
            assert not offset, offset

            mask = (
                (self.df[TRANSCRIPT_ID].isin(transcript_id))
                & (self.df["exon_number"] == float(n))
                & (self.df["strand"].isin(strand))
                & (self.df["feature"] == EXON)
            )
            for _, exon in self.df[mask].iterrows():
                result.append(
                    RnaPosition(
                        contig_id=exon.contig_id,
                        start=exon.transcript_start,
                        start_offset=0,
                        end=exon.transcript_end,
                        end_offset=0,
                        strand=exon.strand,
                        gene_id=exon.gene_id,
                        gene_name=exon.gene_name,
                        transcript_id=exon.transcript_id,
                        transcript_name=exon.transcript_name,
                    )
                )

            return result

        result_start = convert(start, start_offset)
        if start == end:
            return sorted(result_start)
        else:
            result_end = convert(end, end_offset)
            return join_positions(result_start, result_end, TRANSCRIPT_ID)

    def _protein_to_cdna(
        self,
        transcript_id: List[str],
        start: int,
        start_offset: int,
        end: int,
        end_offset: int,
        strand: List[str],
    ) -> List[CdnaPosition]:
        def convert(n: int):
            return ((n - 1) * 3) + 1

        # TODO: Is there a reasonable case where an protein position would have an offset?
        assert not start_offset, start_offset
        assert not end_offset, end_offset

        cdna_start = convert(start)
        cdna_end = convert(end) + 2

        return self._cdna_to_cdna(
            transcript_id, cdna_start, start_offset, cdna_end, end_offset, strand
        )

    def _protein_to_cdna_variant(
        self,
        transcript_id: List[str],
        start: int,
        start_offset: int,
        end: int,
        end_offset: int,
        strand: List[str],
        refseq: str,
        altseq: str,
    ) -> List[_CdnaSmallVariant]:
        result = []

        for cdna in self._protein_to_cdna(
            transcript_id, start, start_offset, end, end_offset, strand
        ):
            result.extend(self._variant_class_from_protein(cdna, refseq, altseq))

        return result

    def _protein_to_dna(
        self,
        transcript_id: List[str],
        start: int,
        start_offset: int,
        end: int,
        end_offset: int,
        strand: List[str],
    ) -> List[DnaPosition]:
        result = []

        for cdna in self._protein_to_cdna(
            transcript_id, start, start_offset, end, end_offset, strand
        ):
            result.extend(
                self._cdna_to_dna(
                    [cdna.transcript_id],
                    cdna.start,
                    cdna.start_offset,
                    cdna.end,
                    cdna.end_offset,
                    [cdna.strand],
                )
            )

        return sorted(set(result))

    def _protein_to_dna_variant(
        self,
        transcript_id: List[str],
        start: int,
        start_offset: int,
        end: int,
        end_offset: int,
        strand: List[str],
        refseq: str,
        altseq: str,
    ) -> List[_DnaSmallVariant]:
        result = []

        for cdna in self._protein_to_cdna_variant(
            transcript_id, start, start_offset, end, end_offset, strand, refseq, altseq
        ):
            result.extend(
                self._cdna_to_dna_variant(
                    [cdna.transcript_id],
                    cdna.start,
                    cdna.start_offset,
                    cdna.end,
                    cdna.end_offset,
                    [cdna.strand],
                    cdna.refseq,
                    cdna.altseq,
                    include_stop=True,
                )
            )

        return sorted(set(result))

    def _protein_to_exon(
        self,
        transcript_id: List[str],
        start: int,
        start_offset: int,
        end: int,
        end_offset: int,
        strand: List[str],
    ) -> List[ExonPosition]:
        result = []

        for cdna in self._protein_to_cdna(
            transcript_id, start, start_offset, end, end_offset, strand
        ):
            result.extend(
                self._cdna_to_exon(
                    [cdna.transcript_id],
                    cdna.start,
                    cdna.start_offset,
                    cdna.end,
                    cdna.end_offset,
                    [cdna.strand],
                )
            )

        return sorted(set(result))

    def _protein_to_exon_variant(
        self,
        transcript_id: List[str],
        start: int,
        start_offset: int,
        end: int,
        end_offset: int,
        strand: List[str],
        refseq: str,
        altseq: str,
    ) -> List[_ExonSmallVariant]:
        result = []

        for cdna in self._protein_to_cdna_variant(
            transcript_id, start, start_offset, end, end_offset, strand, refseq, altseq
        ):
            result.extend(
                self._cdna_to_exon_variant(
                    [cdna.transcript_id],
                    cdna.start,
                    cdna.start_offset,
                    cdna.end,
                    cdna.end_offset,
                    [cdna.strand],
                    cdna.refseq,
                    cdna.altseq,
                    include_stop=True,
                )
            )

        return sorted(set(result))

    def _protein_to_protein(
        self,
        transcript_id: List[str],
        start: int,
        start_offset: int,
        end: int,
        end_offset: int,
        strand: List[str],
    ) -> List[ProteinPosition]:
        result = []

        for cdna in self._protein_to_cdna(
            transcript_id, start, start_offset, end, end_offset, strand
        ):
            result.extend(
                self._cdna_to_protein(
                    [cdna.transcript_id],
                    cdna.start,
                    cdna.start_offset,
                    cdna.end,
                    cdna.end_offset,
                    [cdna.strand],
                )
            )

        return sorted(set(result))

    def _protein_to_protein_variant(
        self,
        transcript_id: List[str],
        start: int,
        start_offset: int,
        end: int,
        end_offset: int,
        strand: List[str],
        refseq: str,
        altseq: str,
    ) -> List[_ProteinSmallVariant]:
        result = []

        for cdna in self._protein_to_cdna_variant(
            transcript_id, start, start_offset, end, end_offset, strand, refseq, altseq
        ):
            result.extend(
                self._cdna_to_protein_variant(
                    [cdna.transcript_id],
                    cdna.start,
                    cdna.start_offset,
                    cdna.end,
                    cdna.end_offset,
                    [cdna.strand],
                    cdna.refseq,
                    cdna.altseq,
                )
            )

        return sorted(set(result))

    def _protein_to_rna(
        self,
        transcript_id: List[str],
        start: int,
        start_offset: int,
        end: int,
        end_offset: int,
        strand: List[str],
    ) -> List[RnaPosition]:
        result = []

        for cdna in self._protein_to_cdna(
            transcript_id, start, start_offset, end, end_offset, strand
        ):
            result.extend(
                self._cdna_to_rna(
                    [cdna.transcript_id],
                    cdna.start,
                    cdna.start_offset,
                    cdna.end,
                    cdna.end_offset,
                    [cdna.strand],
                )
            )

        return sorted(set(result))

    def _protein_to_rna_variant(
        self,
        transcript_id: List[str],
        start: int,
        start_offset: int,
        end: int,
        end_offset: int,
        strand: List[str],
        refseq: str,
        altseq: str,
    ) -> List[_RnaSmallVariant]:
        result = []

        for cdna in self._protein_to_cdna_variant(
            transcript_id, start, start_offset, end, end_offset, strand, refseq, altseq
        ):
            result.extend(
                self._cdna_to_rna_variant(
                    [cdna.transcript_id],
                    cdna.start,
                    cdna.start_offset,
                    cdna.end,
                    cdna.end_offset,
                    [cdna.strand],
                    cdna.refseq,
                    cdna.altseq,
                )
            )

        return sorted(set(result))

    def _rna_to_cdna(
        self,
        transcript_id: List[str],
        start: int,
        start_offset: int,
        end: int,
        end_offset: int,
        strand: List[str],
        include_stop: bool = True,
    ) -> List[CdnaPosition]:
        feature = [CDS, STOP_CODON] if include_stop else [CDS]

        def convert(n: int, offset: int):
            result = []

            # If the offset position can be normalized to a non-offset position, do so. Otherwise
            # just return an offset
            if offset:
                for dna in self._rna_to_dna(
                    transcript_id, start, start_offset, end, end_offset, strand
                ):
                    for cdna in self._dna_to_cdna(
                        [dna.contig_id],
                        dna.start,
                        dna.start_offset,
                        dna.end,
                        dna.end_offset,
                        [dna.strand],
                    ):
                        if cdna.transcript_id in transcript_id:
                            result.append(cdna)

                if result:
                    return result

            mask = (
                (self.df[TRANSCRIPT_ID].isin(transcript_id))
                & (self.df["transcript_start"] <= n)
                & (self.df["transcript_end"] >= n)
                & (self.df["strand"].isin(strand))
                & (self.df["feature"].isin(feature))
            )
            for _, cds in self.df[mask].iterrows():
                new_start = new_end = cds.cdna_start + (n - cds.transcript_start)
                result.append(
                    CdnaPosition(
                        contig_id=cds.contig_id,
                        start=new_start,
                        start_offset=offset,
                        end=new_end,
                        end_offset=offset,
                        strand=cds.strand,
                        gene_id=cds.gene_id,
                        gene_name=cds.gene_name,
                        transcript_id=cds.transcript_id,
                        transcript_name=cds.transcript_name,
                        protein_id=cds.protein_id,
                    )
                )

            return result

        result_start = convert(start, start_offset)
        if start == end:
            return sorted(result_start)
        else:
            result_end = convert(end, end_offset)
            return join_positions(result_start, result_end, TRANSCRIPT_ID)

    def _rna_to_cdna_variant(
        self,
        transcript_id: List[str],
        start: int,
        start_offset: int,
        end: int,
        end_offset: int,
        strand: List[str],
        refseq: str,
        altseq: str,
        include_stop: bool = True,
    ) -> List[_CdnaSmallVariant]:
        result = []

        for cdna in self._rna_to_cdna(
            transcript_id, start, start_offset, end, end_offset, strand, include_stop=include_stop
        ):
            result.extend(self._variant_class_from_cdna(cdna, refseq, altseq))

        return result

    def _rna_to_dna(
        self,
        transcript_id: List[str],
        start: int,
        start_offset: int,
        end: int,
        end_offset: int,
        strand: List[str],
    ) -> List[DnaPosition]:
        def convert(n: int, offset: int):
            result = []

            mask = (
                (self.df[TRANSCRIPT_ID].isin(transcript_id))
                & (self.df["transcript_start"] <= n)
                & (self.df["transcript_end"] >= n)
                & (self.df["strand"].isin(strand))
                & (self.df["feature"] == EXON)
            )
            exon_df = self.df[mask]
            for _, exon in exon_df.iterrows():
                if exon.strand == "-":
                    new_start = new_end = exon.end - (n - exon.transcript_start) - offset
                else:
                    new_start = new_end = exon.start + (n - exon.transcript_start) + offset

                # TODO: Check that new new_start is actually on the contig
                result.append(
                    DnaPosition(
                        contig_id=exon.contig_id,
                        start=new_start,
                        start_offset=0,
                        end=new_end,
                        end_offset=0,
                        strand=exon.strand,
                    )
                )

            return result

        result_start = convert(start, start_offset)
        if start == end:
            return sorted(result_start)
        else:
            result_end = convert(end, end_offset)
            return join_positions(result_start, result_end, CONTIG_ID)

    def _rna_to_dna_variant(
        self,
        transcript_id: List[str],
        start: int,
        start_offset: int,
        end: int,
        end_offset: int,
        strand: List[str],
        refseq: str,
        altseq: str,
    ) -> List[_DnaSmallVariant]:
        result = []

        for dna in self._rna_to_dna(transcript_id, start, start_offset, end, end_offset, strand):
            result.extend(self._dna_small_variant_from_dna(dna, refseq, altseq))

        return result

    def _rna_to_exon(
        self,
        transcript_id: List[str],
        start: int,
        start_offset: int,
        end: int,
        end_offset: int,
        strand: List[str],
    ) -> List[ExonPosition]:
        def convert(n: int, offset: int):
            result = []

            # For an offset position, we need to calculate then equivalent DNA position then map
            # that to one or more exons. This is slower, so if there's no offset we can just map
            # directly to an exon.
            if offset:
                for dna in self._rna_to_dna(
                    transcript_id, start, start_offset, end, end_offset, strand
                ):
                    for exon in self._dna_to_exon(
                        [dna.contig_id],
                        dna.start,
                        dna.start_offset,
                        dna.end,
                        dna.end_offset,
                        [dna.strand],
                    ):
                        if exon.transcript_id in transcript_id:
                            result.append(exon)

                return result

            mask = (
                (self.df[TRANSCRIPT_ID].isin(transcript_id))
                & (self.df["transcript_start"] <= n)
                & (self.df["transcript_end"] >= n)
                & (self.df["strand"].isin(strand))
                & (self.df["feature"] == EXON)
            )
            for _, exon_row in self.df[mask].iterrows():
                result.append(
                    ExonPosition(
                        contig_id=exon_row.contig_id,
                        start=int(exon_row.exon_number),
                        start_offset=0,
                        end=int(exon_row.exon_number),
                        end_offset=0,
                        strand=exon_row.strand,
                        gene_id=exon_row.gene_id,
                        gene_name=exon_row.gene_name,
                        transcript_id=exon_row.transcript_id,
                        transcript_name=exon_row.transcript_name,
                        exon_id=exon_row.exon_id,
                    )
                )

            return result

        result_start = convert(start, start_offset)
        if start == end:
            return sorted(result_start)
        else:
            result_end = convert(end, end_offset)
            return join_positions(result_start, result_end, TRANSCRIPT_ID)

    def _rna_to_exon_variant(
        self,
        transcript_id: List[str],
        start: int,
        start_offset: int,
        end: int,
        end_offset: int,
        strand: List[str],
        refseq: str,
        altseq: str,
    ) -> List[_ExonSmallVariant]:
        result = []

        for exon in self._rna_to_exon(transcript_id, start, start_offset, end, end_offset, strand):
            result.extend(self._exon_small_variant_from_exon(exon, refseq, altseq))

        return result

    def _rna_to_protein(
        self,
        transcript_id: List[str],
        start: int,
        start_offset: int,
        end: int,
        end_offset: int,
        strand: List[str],
    ) -> List[ProteinPosition]:
        result = []

        for cdna in self._rna_to_cdna(
            transcript_id, start, start_offset, end, end_offset, strand, include_stop=False
        ):
            result.extend(
                self._cdna_to_protein(
                    [cdna.transcript_id],
                    cdna.start,
                    cdna.start_offset,
                    cdna.end,
                    cdna.end_offset,
                    [cdna.strand],
                )
            )

        return sorted(set(result))

    def _rna_to_protein_variant(
        self,
        transcript_id: List[str],
        start: int,
        start_offset: int,
        end: int,
        end_offset: int,
        strand: List[str],
        refseq: str,
        altseq: str,
    ) -> List[_ProteinSmallVariant]:
        result = []

        for cdna in self._rna_to_cdna_variant(
            transcript_id,
            start,
            start_offset,
            end,
            end_offset,
            strand,
            refseq,
            altseq,
            include_stop=False,
        ):
            result.extend(
                self._cdna_to_protein_variant(
                    [cdna.transcript_id],
                    cdna.start,
                    cdna.start_offset,
                    cdna.end,
                    cdna.end_offset,
                    [cdna.strand],
                    cdna.refseq,
                    cdna.altseq,
                )
            )

        return sorted(set(result))

    def _rna_to_rna(
        self,
        transcript_id: List[str],
        start: int,
        start_offset: int,
        end: int,
        end_offset: int,
        strand: List[str],
    ) -> List[RnaPosition]:
        def convert(n: int, offset: int):
            result = []

            # If the offset position can be normalized to a non-offset position, do so. Otherwise
            # just return an offset
            if offset:
                for dna in self._rna_to_dna(
                    transcript_id, start, start_offset, end, end_offset, strand
                ):
                    for rna in self._dna_to_rna(
                        [dna.contig_id],
                        dna.start,
                        dna.start_offset,
                        dna.end,
                        dna.end_offset,
                        [dna.strand],
                    ):
                        if rna.transcript_id in transcript_id:
                            result.append(rna)

                return result

            mask = (
                (self.df[TRANSCRIPT_ID].isin(transcript_id))
                & (self.df["transcript_start"] <= n)
                & (self.df["transcript_end"] >= n)
                & (self.df["strand"].isin(strand))
                & (self.df["feature"] == EXON)
            )
            for _, exon in self.df[mask].iterrows():
                result.append(
                    RnaPosition(
                        contig_id=exon.contig_id,
                        start=start,
                        start_offset=offset,
                        end=end,
                        end_offset=offset,
                        strand=exon.strand,
                        gene_id=exon.gene_id,
                        gene_name=exon.gene_name,
                        transcript_id=exon.transcript_id,
                        transcript_name=exon.transcript_name,
                    )
                )

            return result

        result_start = convert(start, start_offset)
        if start == end:
            return sorted(result_start)
        else:
            result_end = convert(end, end_offset)
            return join_positions(result_start, result_end, TRANSCRIPT_ID)

    def _rna_to_rna_variant(
        self,
        transcript_id: List[str],
        start: int,
        start_offset: int,
        end: int,
        end_offset: int,
        strand: List[str],
        refseq: str,
        altseq: str,
    ) -> List[_RnaSmallVariant]:
        result = []

        for rna in self._rna_to_rna(transcript_id, start, start_offset, end, end_offset, strand):
            result.extend(self._rna_small_variant_from_rna(rna, refseq, altseq))

        return result

    def _variant_class_from_cdna(
        self, cdna: CdnaPosition, refseq: str, altseq: str
    ) -> List[_CdnaSmallVariant]:
        """Convert a cDNA position plus ref/alt nucleotides into a cDNA variant.

        Args:
            cdna (CdnaPosition): _Variant position
            refseq (str): feature allele
            altseq (str): alternate allele

        Raises:
            NotImplementedError: An unsupported combination of feature/alternate alleles was given
            ValueError: The given feature allele does not match the annotated feature allele

        Returns:
            List[_CdnaSmallVariant]: One or more cDNA variants
        """
        variant_list = []

        ref_annotated = self.sequence(cdna)
        for ref, alt in product(expand_nt(refseq), expand_nt(altseq)):
            # Assert that the given ref matches the annotated one
            if ref != ref_annotated:
                raise ValueError(
                    f"Given ref allele '{ref}' does not match annotated ref allele '{ref_annotated}'"
                )

            # Trim bases that are unchanged between the ref and alt alleles
            new_ref, new_alt, new_start, new_end = collapse_seq_change(ref, alt)
            start = cdna.start + new_start
            end = cdna.end - new_end

            # Determine the type of variant
            if is_deletion(new_ref, new_alt):
                variant = CdnaDeletion.copy_from(
                    cdna, start=start, end=end, refseq=new_ref, altseq=new_alt
                )
            elif is_delins(new_ref, new_alt):
                variant = CdnaDelins.copy_from(
                    cdna, start=start, end=end, refseq=new_ref, altseq=new_alt
                )
            elif is_duplication(new_ref, new_alt):
                variant = CdnaDuplication.copy_from(
                    cdna, start=start, end=end, refseq=new_ref, altseq=new_alt
                )
            elif is_insertion(new_ref, new_alt):
                variant = CdnaInsertion.copy_from(
                    cdna, start=start, end=end, refseq=new_ref, altseq=new_alt
                )
            elif is_substitution(new_ref, new_alt):
                variant = CdnaSubstitution.copy_from(
                    cdna, start=start, end=end, refseq=new_ref, altseq=new_alt
                )
            else:
                raise NotImplementedError(
                    f"Unrecognized sequence change '{new_ref}/{new_alt}' at {cdna}"
                )

            variant_list.append(variant)

        return variant_list

    def _variant_class_from_protein(
        self, cdna: CdnaPosition, refaa: str, altaa: str
    ) -> List[_CdnaSmallVariant]:
        """Convert a cDNA position plus ref/alt amino acids into a cDNA variant.

        Args:
            cdna (CdnaPosition): cDNA position
            refseq (str): Reference allele
            altaa (str): Alternate allele

        Raises:
            NotImplementedError: An unsupported combination of feature/alternate alleles was given

        Returns:
            List[_CdnaSmallVariant]: One or more cDNA variants
        """
        variant_list = []

        ref_annotated = self.sequence(cdna)
        for ref, alt in product(reverse_translate(refaa), reverse_translate(altaa)):
            # Assert that the given ref matches the annotated one
            if ref != ref_annotated:
                continue

            # For insertions, check that the sequence flanking the inserted sequence matches the ref
            if is_insertion(refaa, altaa) and not is_insertion(ref, alt):
                continue

            # Trim bases that are unchanged between the ref and alt alleles
            new_ref, new_alt, new_start, new_end = collapse_seq_change(ref, alt)
            start = cdna.start + new_start
            end = cdna.end - new_end

            # Determine the type of variant
            if is_deletion(new_ref, new_alt):
                variant = CdnaDeletion.copy_from(
                    cdna, start=start, end=end, refseq=new_ref, altseq=new_alt
                )
            elif is_delins(new_ref, new_alt):
                variant = CdnaDelins.copy_from(
                    cdna, start=start, end=end, refseq=new_ref, altseq=new_alt
                )
            elif is_duplication(new_ref, new_alt):
                variant = CdnaDuplication.copy_from(
                    cdna, start=start, end=end, refseq=new_ref, altseq=new_alt
                )
            elif is_insertion(new_ref, new_alt):
                variant = CdnaInsertion.copy_from(
                    cdna, start=start, end=end, refseq=new_ref, altseq=new_alt
                )
            elif is_substitution(new_ref, new_alt):
                variant = CdnaSubstitution.copy_from(
                    cdna, start=start, end=end, refseq=new_ref, altseq=new_alt
                )
            else:
                raise NotImplementedError(
                    f"Unrecognized sequence change '{new_ref}/{new_alt}' at {cdna}"
                )

            variant_list.append(variant)

        return variant_list

    def _dna_small_variant_from_dna(
        self, dna: DnaPosition, refseq: str, altseq: str
    ) -> List[_DnaSmallVariant]:
        """Convert a DNA position plus ref/alt nucleotides into a DNA variant.

        Args:
            dna (DnaPosition): DNA position
            refseq (str): Reference allele
            altseq (str): Alternate allele

        Raises:
            NotImplementedError: An unsupported combination of feature/alternate alleles was given

        Returns:
            List[_DnaSmallVariant]: One or more DNA variants
        """
        variant_list = []

        # TODO: DNA sequence get is slow
        ref_annotated = self.sequence(dna)
        for ref, alt in product(expand_nt(refseq), expand_nt(altseq)):
            # Assert that the given ref matches the annotated one
            if ref != ref_annotated:
                continue

            # Trim bases that are unchanged between the ref and alt alleles
            new_ref, new_alt, new_start, new_end = collapse_seq_change(ref, alt)
            start = dna.start + new_start
            end = dna.end - new_end

            # Determine the type of variant
            if is_deletion(new_ref, new_alt):
                variant = DnaDeletion.copy_from(
                    dna, start=start, end=end, refseq=new_ref, altseq=new_alt
                )
            elif is_delins(new_ref, new_alt):
                variant = DnaDelins.copy_from(
                    dna, start=start, end=end, refseq=new_ref, altseq=new_alt
                )
            elif is_duplication(new_ref, new_alt):
                variant = DnaDuplication.copy_from(
                    dna, start=start, end=end, refseq=new_ref, altseq=new_alt
                )
            elif is_insertion(new_ref, new_alt):
                variant = DnaInsertion.copy_from(
                    dna, start=start, end=end, refseq=new_ref, altseq=new_alt
                )
            elif is_substitution(new_ref, new_alt):
                variant = DnaSubstitution.copy_from(
                    dna, start=start, end=end, refseq=new_ref, altseq=new_alt
                )
            else:
                raise NotImplementedError(
                    f"Unrecognized sequence change '{new_ref}/{new_alt}' at {dna}"
                )

            variant_list.append(variant)

        return variant_list

    def _exon_small_variant_from_exon(
        self, exon: ExonPosition, refseq: str, altseq: str
    ) -> List[_ExonSmallVariant]:
        """Convert an exon position plus ref/alt nucleotides into an exon variant.

        Args:
            exon (ExonPosition): exon position
            refseq (str): Reference allele
            altseq (str): Alternate allele

        Raises:
            NotImplementedError: An unsupported combination of feature/alternate alleles was given

        Returns:
            List[_ExonSmallVariant]: One or more exon variants
        """
        raise NotImplementedError()  # TODO

    def _protein_small_variant(
        self, cdna: CdnaPosition, protein: ProteinPosition, refseq: str, altseq: str
    ) -> List[_ProteinSmallVariant]:
        """Convert a cDNA position plus ref/alt amino acids into a protein variant.

        Args:
            cdna (CdnaPosition): cDNA position
            refseq (str): Reference allele
            altseq (str): Alternate allele

        Raises:
            NotImplementedError: An unsupported combination of feature/alternate alleles was given

        Returns:
            List[_ProteinSmallVariant]: One or more protein variants
        """
        variant_list = []

        protein_refseq = self.sequence(protein)
        for cdna_ref, cdna_alt in product(expand_nt(refseq), expand_nt(altseq)):
            if is_frameshift(cdna_ref, cdna_alt):
                raise NotImplementedError()  # TODO

            for protein_alt in self.translate_cdna_variant(cdna, cdna_alt):
                # Trim bases that are unchanged between the ref and alt alleles
                new_ref, new_alt, new_start, new_end = collapse_seq_change(
                    protein_refseq, protein_alt
                )
                start = protein.start + new_start
                end = protein.end - new_end

                # Determine the type of variant
                if is_deletion(new_ref, new_alt):
                    variant = ProteinDeletion.copy_from(
                        protein, start=start, end=end, refseq=new_ref, altseq=new_alt
                    )
                elif is_delins(new_ref, new_alt):
                    variant = ProteinDelins.copy_from(
                        protein, start=start, end=end, refseq=new_ref, altseq=new_alt
                    )
                elif is_duplication(new_ref, new_alt):
                    variant = ProteinDuplication.copy_from(
                        protein, start=start, end=end, refseq=new_ref, altseq=new_alt
                    )
                elif is_insertion(new_ref, new_alt):
                    variant = ProteinInsertion.copy_from(
                        protein, start=start, end=end, refseq=new_ref, altseq=new_alt
                    )
                elif is_substitution(new_ref, new_alt):
                    variant = ProteinSubstitution.copy_from(
                        protein, start=start, end=end, refseq=new_ref, altseq=new_alt
                    )
                else:
                    raise NotImplementedError(
                        f"Unrecognized sequence change '{new_ref}/{new_alt}' at {protein}"
                    )

                variant_list.append(variant)

        return variant_list

    def _rna_small_variant_from_rna(
        self, rna: RnaPosition, refseq: str, altseq: str
    ) -> List[_RnaSmallVariant]:
        """Convert an RNA position plus ref/alt nucleotides into an RNA variant.

        Args:
            rna (RnaPosition): RNA position
            refseq (str): Reference allele
            altseq (str): Alternate allele

        Raises:
            NotImplementedError: An unsupported combination of feature/alternate alleles was given
            ValueError: The given feature allele does not match the annotated feature allele

        Returns:
            List[_RnaSmallVariant]: One or more RNA variants
        """
        variant_list = []

        ref_annotated = self.sequence(rna)
        for ref, alt in product(expand_nt(refseq), expand_nt(altseq)):
            # Assert that the given ref matches the annotated one
            if ref != ref_annotated:
                raise ValueError(
                    f"Given ref allele '{ref}' does not match annotated ref allele '{ref_annotated}' for {rna}"
                )

            # Trim bases that are unchanged between the ref and alt alleles
            new_ref, new_alt, new_start, new_end = collapse_seq_change(ref, alt)
            start = rna.start + new_start
            end = rna.end - new_end

            # Determine the type of variant
            if is_deletion(new_ref, new_alt):
                variant = RnaDeletion.copy_from(
                    rna, start=start, end=end, refseq=new_ref, altseq=new_alt
                )
            elif is_delins(new_ref, new_alt):
                variant = RnaDelins.copy_from(
                    rna, start=start, end=end, refseq=new_ref, altseq=new_alt
                )
            elif is_duplication(new_ref, new_alt):
                variant = RnaDuplication.copy_from(
                    rna, start=start, end=end, refseq=new_ref, altseq=new_alt
                )
            elif is_insertion(new_ref, new_alt):
                variant = RnaInsertion.copy_from(
                    rna, start=start, end=end, refseq=new_ref, altseq=new_alt
                )
            elif is_substitution(new_ref, new_alt):
                variant = RnaSubstitution.copy_from(
                    rna, start=start, end=end, refseq=new_ref, altseq=new_alt
                )
            else:
                raise NotImplementedError(
                    f"Unrecognized sequence change '{new_ref}/{new_alt}' at {rna}"
                )

            variant_list.append(variant)

        return variant_list

    # ---------------------------------------------------------------------------------------------
    # Functions for getting feature ID/names
    # ---------------------------------------------------------------------------------------------
    def contig_ids(self, feature: str = "") -> List[str]:
        """Return the contig IDs that map to the given feature. If no feature is given, return
        all contig IDs.

        Examples:
            >>> ensembl100.contig_ids("BRCA2")
            ['13']

        Args:
            feature (str, optional): Feature ID or name

        Returns:
            List[str]: Contig IDs
        """
        return self._query_feature(CONTIG_ID, feature)

    def exon_ids(self, feature: str = "") -> List[str]:
        """Return the exon IDs that map to the given feature. If no feature is given, return
        all exon IDs.

        Examples:
            >>> ensembl100.exon_ids("BRCA2")[:3]
            ['ENSE00000939167', 'ENSE00000939168', 'ENSE00000939169']

        Args:
            feature (str, optional): Feature ID or name

        Returns:
            List[str]: Exon IDs
        """
        return self._query_feature(EXON_ID, feature)

    def gene_ids(self, feature: str = "") -> List[str]:
        """Return the gene IDs that map to the given feature. If no feature is given, return
        all gene IDs.

        Examples:
            >>> ensembl100.gene_ids("BRCA2")
            ['ENSG00000139618']

        Args:
            feature (str, optional): Feature ID or name

        Returns:
            List[str]: Gene IDs
        """
        return self._query_feature(GENE_ID, feature)

    def gene_names(self, feature: str = "") -> List[str]:
        """Return the gene names that map to the given feature. If no feature is given, return
        all gene names.

        Examples:
            >>> ensembl100.gene_names("ENSG00000139618")
            ['BRCA2']

        Args:
            feature (str, optional): Feature ID or name

        Returns:
            List[str]: Gene names
        """
        return self._query_feature(GENE_NAME, feature)

    def protein_ids(self, feature: str = "") -> List[str]:
        """Return the protein IDs that map to the given feature. If no feature is given, return
        all protein IDs.

        Examples:
            >>> ensembl100.protein_ids("BRCA2")[:3]
            ['ENSP00000369497', 'ENSP00000433168', 'ENSP00000434898']

        Args:
            feature (str, optional): Feature ID or name

        Returns:
            List[str]: Protein IDs
        """
        return self._query_feature(PROTEIN_ID, feature)

    def transcript_ids(self, feature: str = "") -> List[str]:
        """Return the transcript IDs that map to the given feature. If no feature is given, return
        all transcript IDs.

        Examples:
            >>> ensembl100.transcript_ids("BRCA2")[:3]
            ['ENST00000380152', 'ENST00000470094', 'ENST00000528762']

        Args:
            feature (str, optional): Feature ID or name

        Returns:
            List[str]: Transcript IDs
        """
        return self._query_feature(TRANSCRIPT_ID, feature)

    def transcript_names(self, feature: str = "") -> List[str]:
        """Return the transcript names that map to the given feature. If no feature is given, return
        all transcript names.

        Examples:
            >>> ensembl100.transcript_names("BRCA2")[:3]
            ['BRCA2-201', 'BRCA2-202', 'BRCA2-203']

        Args:
            feature (str, optional): Feature ID or name

        Returns:
            List[str]: Transcript names
        """
        return self._query_feature(TRANSCRIPT_NAME, feature)

    def _query_feature(self, key: str, feature: str = "") -> List[str]:
        if feature:
            parts = []

            for feature, feature_type in self.normalize_id(feature):
                if feature_type == CONTIG_ID:
                    func = self._query_contig_id
                elif feature_type == EXON_ID:
                    func = self._query_exon_id
                elif feature_type == GENE_ID:
                    func = self._query_gene_id
                elif feature_type == GENE_NAME:
                    func = self._query_gene_name
                elif feature_type == PROTEIN_ID:
                    func = self._query_protein_id
                elif feature_type == TRANSCRIPT_ID:
                    func = self._query_transcript_id
                elif feature_type == TRANSCRIPT_NAME:
                    func = self._query_transcript_name
                else:
                    raise ValueError(f"Unable to get {key} for {feature} ({feature_type})")

                parts.append(func(feature)[key])

            if parts:
                result = pd.concat(parts)
            else:
                result = pd.Series(dtype="object")
        else:
            result = self.df[key]

        return self._uniquify_series(result)

    def _query_contig_id(self, feature: Union[List[str], str]) -> pd.DataFrame:
        return self._query(feature, CONTIG_ID)

    def _query_exon_id(self, feature: Union[List[str], str]) -> pd.DataFrame:
        return self._query(feature, EXON_ID)

    def _query_gene_id(self, feature: Union[List[str], str]) -> pd.DataFrame:
        return self._query(feature, GENE_ID)

    def _query_gene_name(self, feature: Union[List[str], str]) -> pd.DataFrame:
        return self._query(feature, GENE_NAME)

    def _query_protein_id(self, feature: Union[List[str], str]) -> pd.DataFrame:
        return self._query(feature, PROTEIN_ID)

    def _query_transcript_id(self, feature: Union[List[str], str]) -> pd.DataFrame:
        return self._query(feature, TRANSCRIPT_ID)

    def _query_transcript_name(self, feature: Union[List[str], str]) -> pd.DataFrame:
        return self._query(feature, TRANSCRIPT_NAME)

    def _query(self, feature: Union[List[str], str], col: str) -> pd.DataFrame:
        feature = [feature] if isinstance(feature, str) else feature
        sudbf = self.df.loc[self.df[col].isin(feature)]

        return sudbf

    def _uniquify_series(self, series: pd.Series) -> List:
        return sorted(series.dropna().unique().tolist())

    # ---------------------------------------------------------------------------------------------
    # Functions for getting feature aliases
    # ---------------------------------------------------------------------------------------------
    def contig_alias(self, contig_id: str) -> List[str]:
        """List all aliases of the given contig ID.

        Args:
            contig_id (str): contig ID

        Returns:
            List[str]: All aliases of the given contig ID
        """
        return self._alias(contig_id, self._contig_alias)

    def exon_alias(self, exon_id: str) -> List[str]:
        """List all aliases of the given exon ID.

        Args:
            exon_id (str): exon ID

        Returns:
            List[str]: All aliases of the given exon ID
        """
        return self._alias(exon_id, self._exon_alias)

    def gene_alias(self, gene_id: str) -> List[str]:
        """List all aliases of the given gene ID.

        Args:
            gene_id (str): gene ID

        Returns:
            List[str]: All aliases of the given gene ID
        """
        return self._alias(gene_id, self._gene_alias)

    def protein_alias(self, protein_id: str) -> List[str]:
        """List all aliases of the given protein ID.

        Args:
            protein_id (str): protein ID

        Returns:
            List[str]: All aliases of the given protein ID
        """
        return self._alias(protein_id, self._protein_alias)

    def transcript_alias(self, transcript_id: str) -> List[str]:
        """List all aliases of the given transcript ID.

        Args:
            transcript_id (str): transcript ID

        Returns:
            List[str]: All aliases of the given transcript ID
        """
        return self._alias(transcript_id, self._transcript_alias)

    def _alias(self, feature: str, alias_dict: Dict[str, List[str]]) -> List[str]:
        # Try different variations on the feature ID until one is found
        for key in [
            feature,
            strip_version(feature),
            feature.lower(),
            strip_version(feature).lower(),
            feature.upper(),
            strip_version(feature).upper(),
        ]:
            if alias := alias_dict.get(key, []):
                # Coerce the alias to a list
                if isinstance(alias, str):
                    return [alias]
                else:
                    return list(alias)
        else:
            return []

    # ---------------------------------------------------------------------------------------------
    # Functions for normalizing feature ID/names
    # ---------------------------------------------------------------------------------------------
    def normalize_id(self, feature: str) -> List[Tuple[str, str]]:
        """Normalize an ID or name to the annotated equivalent(s).

        Examples:
            >>> ensembl100.normalize_id("BRCA2-201")
            [('BRCA2-201', 'transcript_name')]

        Args:
            feature (str): Feature ID or name

        Returns:
            List[Tuple[str, str]]: List of (normalized ID/name, normalized type (e.g. 'gene',
                'transcript', etc.))
        """
        normalized = []

        feature_type, result = self._normalize_id(feature)
        if feature_type:
            normalized = self._uniquify_series(result[feature_type])

        return [(i, feature_type) for i in normalized]

    def _normalize_id(self, feature: str) -> Tuple[str, pd.DataFrame]:
        feature_type = ""
        result = pd.DataFrame()

        for key, func in [
            (CONTIG_ID, self._normalize_contig_id),
            (EXON_ID, self._normalize_exon_id),
            (GENE_ID, self._normalize_gene_id),
            (GENE_NAME, self._normalize_gene_name),
            (PROTEIN_ID, self._normalize_protein_id),
            (TRANSCRIPT_ID, self._normalize_transcript_id),
            (TRANSCRIPT_NAME, self._normalize_transcript_name),
        ]:
            result = func(feature)
            if not result.empty:
                feature_type = key
                break

        return feature_type, result

    def _normalize_contig_id(self, feature: str) -> pd.DataFrame:
        featurel = [feature] + self.contig_alias(feature)

        return self._query_contig_id(featurel)

    def _normalize_exon_id(self, feature: str) -> pd.DataFrame:
        featurel = [feature] + self.exon_alias(feature)

        return self._query_exon_id(featurel)

    def _normalize_gene_id(self, feature: str) -> pd.DataFrame:
        featurel = [feature] + self.gene_alias(feature)

        return self._query_gene_id(featurel)

    def _normalize_gene_name(self, feature: str) -> pd.DataFrame:
        featurel = [feature] + self.gene_alias(feature)

        return self._query_gene_name(featurel)

    def _normalize_protein_id(self, feature: str) -> pd.DataFrame:
        featurel = [feature] + self.protein_alias(feature)

        return self._query_protein_id(featurel)

    def _normalize_transcript_id(self, feature: str) -> pd.DataFrame:
        featurel = [feature] + self.transcript_alias(feature)

        return self._query_transcript_id(featurel)

    def _normalize_transcript_name(self, feature: str) -> pd.DataFrame:
        featurel = [feature] + self.transcript_alias(feature)

        return self._query_transcript_name(featurel)

    # ---------------------------------------------------------------------------------------------
    # Functions for checking feature type
    # ---------------------------------------------------------------------------------------------
    def is_contig(self, feature: str) -> bool:
        """Check if the given ID or name is a contig.

        Args:
            feature (str): Feature ID or name

        Returns:
            bool: True if is a contig else False
        """
        return any((i[1] == CONTIG_ID for i in self.normalize_id(feature)))

    def is_exon(self, feature: str) -> bool:
        """Check if the given ID or name is an exon.

        Args:
            feature (str): Feature ID or name

        Returns:
            bool: True if is an exon else False
        """
        return any((i[1] == EXON_ID for i in self.normalize_id(feature)))

    def is_gene(self, feature: str) -> bool:
        """Check if the given ID or name is a gene.

        Args:
            feature (str): Feature ID or name

        Returns:
            bool: True if is an gene else False
        """
        return any((i[1] in (GENE_ID, GENE_NAME) for i in self.normalize_id(feature)))

    def is_protein(self, feature: str) -> bool:
        """Check if the given ID or name is a protein.

        Args:
            feature (str): Feature ID or name

        Returns:
            bool: True if is an protein else False
        """
        return any((i[1] == PROTEIN_ID for i in self.normalize_id(feature)))

    def is_transcript(self, feature: str) -> bool:
        """Check if the given ID or name is a transcript.

        Args:
            feature (str): Feature ID or name

        Returns:
            bool: True if is an transcript else False
        """
        return any((i[1] in (TRANSCRIPT_ID, TRANSCRIPT_NAME) for i in self.normalize_id(feature)))

    # ---------------------------------------------------------------------------------------------
    # Functions for checking if a transcript is a canonical transcript
    # ---------------------------------------------------------------------------------------------
    def is_canonical_transcript(self, transcript_id: str) -> bool:
        """Check if the given transcript ID is the canonical transcript for its gene.

        Args:
            transcript_id (str): transcript ID

        Returns:
            bool: True if the transcript is the canonical transcript else False
        """
        return transcript_id in self._canonical_transcript

    # ---------------------------------------------------------------------------------------------
    # <feature>
    # ---------------------------------------------------------------------------------------------
    def cdna(self, feature: str, canonical: bool = False) -> List[CdnaPosition]:
        """Return the cDNA position(s) matching the given feature ID or name.

        Args:
            feature (str): Feature ID or name
            canonical (bool, optional): Only consider the canonical transcript when mapping. Defaults to False.

        Returns:
            List[CdnaPosition]: One or more cDNA positions.
        """
        result = []

        transcript_ids = self.transcript_ids(feature)
        mask = (self.df[TRANSCRIPT_ID].isin(transcript_ids)) & (self.df["feature"] == CDNA)
        for _, cdna in self.df[mask].iterrows():
            if canonical and not self.is_canonical_transcript(cdna.transcript_id):
                continue

            result.append(
                CdnaPosition(
                    contig_id=cdna.contig_id,
                    start=cdna.cdna_start,
                    start_offset=0,
                    end=cdna.cdna_end,
                    end_offset=0,
                    strand=cdna.strand,
                    gene_id=cdna.gene_id,
                    gene_name=cdna.gene_name,
                    transcript_id=cdna.transcript_id,
                    transcript_name=cdna.transcript_name,
                    protein_id=cdna.protein_id,
                )
            )

        return sorted(set(result))

    def dna(self, feature: str) -> List[DnaPosition]:
        """Return the DNA position(s) matching the given feature ID or name.

        Args:
            feature (str): Feature ID or name
            canonical (bool, optional): Only consider the canonical transcript when mapping. Defaults to False.

        Returns:
            List[DnaPosition]: One or more DNA positions.
        """
        result = []

        # get the strand of the original feature
        _, df = self._normalize_id(feature)
        strand_list = self._uniquify_series(df["strand"])

        for contig_id in self.contig_ids(feature):
            for fasta in self.dna_fasta:
                try:
                    contig_seq = fasta[contig_id]
                    break
                except KeyError:
                    continue
            else:
                raise KeyError(f"Sequence '{contig_id}' not found")

            start = 1
            end = len(contig_seq)
            for strand in strand_list:
                result.append(
                    DnaPosition(
                        contig_id=contig_id,
                        start=start,
                        start_offset=0,
                        end=end,
                        end_offset=0,
                        strand=strand,
                    )
                )

        return sorted(set(result))

    def exon(self, feature: str, canonical: bool = False) -> List[ExonPosition]:
        """Return the exon position(s) matching the given feature ID or name.

        Args:
            feature (str): Feature ID or name
            canonical (bool, optional): Only consider the canonical transcript when mapping. Defaults to False.

        Returns:
            List[ExonPosition]: One or more exon positions.
        """
        result = []

        exon_ids = self.exon_ids(feature)
        mask = (self.df[EXON_ID].isin(exon_ids)) & (self.df["feature"] == EXON)
        for _, exon in self.df[mask].iterrows():
            if canonical and not self.is_canonical_transcript(exon.transcript_id):
                continue

            result.append(
                ExonPosition(
                    contig_id=exon.contig_id,
                    start=int(exon.exon_number),
                    start_offset=0,
                    end=int(exon.exon_number),
                    end_offset=0,
                    strand=exon.strand,
                    gene_id=exon.gene_id,
                    gene_name=exon.gene_name,
                    transcript_id=exon.transcript_id,
                    transcript_name=exon.transcript_name,
                    exon_id=exon.exon_id,
                )
            )

        return sorted(set(result))

    def gene(self, feature: str) -> List[DnaPosition]:
        """Return the gene position(s) matching the given feature ID or name.

        Args:
            feature (str): Feature ID or name
            canonical (bool, optional): Only consider the canonical transcript when mapping. Defaults to False.

        Returns:
            List[DnaPosition]: One or more DNA positions.
        """
        result = []

        gene_ids = self.gene_ids(feature)
        mask = (self.df[GENE_ID].isin(gene_ids)) & (self.df["feature"] == "gene")
        for _, gene in self.df[mask].iterrows():
            result.append(
                DnaPosition(
                    contig_id=gene.contig_id,
                    start=gene.start,
                    start_offset=0,
                    end=gene.end,
                    end_offset=0,
                    strand=gene.strand,
                )
            )

        return sorted(set(result))

    def protein(self, feature: str, canonical: bool = False) -> List[ProteinPosition]:
        """Return the protein position(s) matching the given feature ID or name.

        Args:
            feature (str): Feature ID or name
            canonical (bool, optional): Only consider the canonical transcript when mapping. Defaults to False.

        Returns:
            List[ProteinPosition]: One or more protein positions.
        """
        result = []
        for cdna in self.cdna(feature, canonical=canonical):
            # Convert the cDNA position to a protein position
            protein_start = calc_cdna_to_protein(cdna.start)
            protein_end = calc_cdna_to_protein(cdna.end)
            result.append(
                ProteinPosition.copy_from(
                    cdna, start=protein_start, start_offset=0, end=protein_end, end_offset=0
                )
            )

        return sorted(set(result))

    def rna(self, feature: str, canonical: bool = False) -> List[RnaPosition]:
        """Return the RNA position(s) matching the given feature ID or name.

        Args:
            feature (str): Feature ID or name
            canonical (bool, optional): Only consider the canonical transcript when mapping. Defaults to False.

        Returns:
            List[RnaPosition]: One or more RNA positions.
        """
        result = []

        transcript_ids = self.transcript_ids(feature)
        mask = (self.df[TRANSCRIPT_ID].isin(transcript_ids)) & (self.df["feature"] == "transcript")
        for _, transcript in self.df[mask].iterrows():
            if canonical and not self.is_canonical_transcript(transcript.transcript_id):
                continue

            result.append(
                RnaPosition(
                    contig_id=transcript.contig_id,
                    start=transcript.transcript_start,
                    start_offset=0,
                    end=transcript.transcript_end,
                    end_offset=0,
                    strand=transcript.strand,
                    gene_id=transcript.gene_id,
                    gene_name=transcript.gene_name,
                    transcript_id=transcript.transcript_id,
                    transcript_name=transcript.transcript_name,
                )
            )

        return sorted(set(result))

    # ---------------------------------------------------------------------------------------------
    # Utility functions
    # ---------------------------------------------------------------------------------------------
    def translate_cdna_variant(self, cdna: CdnaPosition, cdna_altseq: str) -> List[str]:
        """Return the mutated protein sequence, given a cDNA position and alt allele.

        Args:
            cdna (CdnaPosition): cDNA position object
            cdna_altseq (str): cDNA alternate allele

        Returns:
            List[str]: protein alternate allele(s)
        """
        pep_altseq_set = set()

        # If no alt, assume we're talking about a deletion variant
        if not cdna_altseq:
            return [""]

        # Get the codon sequence
        codon_start_offset = (cdna.start - 1) % 3
        codon_start = cdna.start - codon_start_offset
        codon_end_offset = 2 - ((cdna.end - 1) % 3)
        codon_end = cdna.end + codon_end_offset
        # Create a new CdnaPosition encompassing the whole codon(s)
        codon = CdnaPosition.copy_from(cdna, start=codon_start, end=codon_end)
        codon_refseq = self.sequence(codon)
        # Assert that the codon sequence is divisible by 3
        assert len(codon_refseq) % 3 == 0

        # Mutate the codon sequence
        codon_refseq_left = codon_refseq[:codon_start_offset]
        codon_refseq_right = codon_refseq[-codon_end_offset:] if codon_end_offset else ""
        for i in expand_nt(cdna_altseq):
            codon_altseq = codon_refseq_left + i + codon_refseq_right
            # Assert that the altered codon sequence is divisible by 3
            assert len(codon_altseq) % 3 == 0, codon_altseq
            pep_altseq = "".join(AMINO_ACID_TABLE[codon] for codon in split_by_codon(codon_altseq))
            pep_altseq_set.add(pep_altseq)

        return sorted(pep_altseq_set)


def join_positions(start: List, end: List, merge_on: str) -> List:
    """Return the combination of two list of position or variant objects - one of start positions
    and one of end positions - into one list by the given key (e.g. 'transcript_id').

    All positions must be of the same class.

    Args:
        start (List[position or small variant]): Start positions
        end (List[position or small variant]): End positions
        merge_on (str): Attribute to merge on (e.g. 'transcript_id')

    Returns:
        List[position or small variant]: New position or variant objects
    """
    result = set()

    for start_pos, end_pos in product(start, end):
        # Make sure that the positions are sorted by position
        if start_pos.start > end_pos.start:
            start_pos, end_pos = end_pos, start_pos
        # Make a new position from the start of the 1st position and end of the 2nd
        start_key = getattr(start_pos, merge_on)
        end_key = getattr(end_pos, merge_on)
        if start_key == end_key:
            assert start_pos.__class__ == end_pos.__class__
            new = start_pos.copy_from(start_pos, end=end_pos.end, end_offset=end_pos.end_offset)
            result.add(new)

    return sorted(result)  # type: ignore
