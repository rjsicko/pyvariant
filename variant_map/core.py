"""Core logic for handling annotations and mapping positions/variants between different types."""
from __future__ import annotations

from itertools import product
from typing import Dict, List, Optional, Tuple, Union, cast

import pandas as pd
from gtfparse import read_gtf
from pyfaidx import Fasta

from .constants import (
    CONTIG_ID,
    EXON_ID,
    GENE_ID,
    GENE_NAME,
    PROTEIN_ID,
    TRANSCRIPT_ID,
    TRANSCRIPT_NAME,
)
from .files import read_fasta, tsv_to_dict, txt_to_list
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
    _Position,
    _ProteinSmallVariant,
    _RnaSmallVariant,
    _SmallVariant,
)
from .tables import AMINO_ACID_TABLE
from .types import Position
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
        canonical_transcript: str = "",
        contig_alias: str = "",
        exon_alias: str = "",
        gene_alias: str = "",
        protein_alias: str = "",
        transcript_alias: str = "",
    ):
        """
        Args:
            gtf (str): Path to a GTF files with feature annotations
            cds (List[str]): List of paths to a FASTA files of CDS sequences
            dna (List[str]): List of paths to a FASTA files of DNA sequences
            peptide (List[str]): List of paths to a FASTA files of peptide sequences
            rna (List[str]): List of paths to a FASTA files of RNA sequences
            canonical_transcript (str, optional): Path to a text file of canonical transcript IDs.
                Defaults to "".
            contig_alias (str, optional): Path to a TSV file mapping contig IDs to their alias(es).
                Defaults to "".
            exon_alias (str, optional): Path to a TSV file mapping contig IDs to their alias(es).
                Defaults to "".
            gene_alias (str, optional): Path to a TSV file mapping contig IDs to their alias(es).
                Defaults to "".
            protein_alias (str, optional): Path to a TSV file mapping contig IDs to their alias(es).
                Defaults to "".
            transcript_alias (str, optional): Path to a TSV file mapping contig IDs to their alias(es).
                Defaults to "".
        """
        self.df = read_gtf(gtf, result_type="pandas")  # TODO: switch to 'polars'?
        self.cds_fasta = [read_fasta(i) for i in cds]
        self.dna_fasta = [read_fasta(i) for i in dna]
        self.protein_fasta = [read_fasta(i) for i in peptide]
        self.rna_fasta = [read_fasta(i) for i in rna]
        self._canonical_transcript = txt_to_list(canonical_transcript)
        self._contig_alias = tsv_to_dict(contig_alias)
        self._exon_alias = tsv_to_dict(exon_alias)
        self._gene_alias = tsv_to_dict(gene_alias)
        self._protein_alias = tsv_to_dict(protein_alias)
        self._transcript_alias = tsv_to_dict(transcript_alias)

    # ---------------------------------------------------------------------------------------------
    # all_<feature_symbol>s
    # ---------------------------------------------------------------------------------------------
    def all_contig_ids(self) -> List[str]:
        """List all contig (chromosome) IDs.

        Examples:
            >>> ensembl100.all_contig_ids()[:3]
            ['1', '10', '11']

        Returns:
            List[str]: Contig IDs
        """
        return self._uniquify_series(self.df[CONTIG_ID])

    def all_exon_ids(self) -> List[str]:
        """List all exon IDs.

        Examples:
            >>> ensembl100.all_exon_ids()[:3]
            ['ENSE00000000001', 'ENSE00000000002', 'ENSE00000000003']

        Returns:
            List[str]: Exon IDs
        """
        return self._uniquify_series(self.df[EXON_ID])

    def all_gene_ids(self) -> List[str]:
        """List all gene IDs.

        Examples:
            >>> ensembl100.all_gene_ids()[:3]
            ['ENSG00000000003', 'ENSG00000000005', 'ENSG00000000419']

        Returns:
            List[str]: Gene IDs
        """
        return self._uniquify_series(self.df[GENE_ID])

    def all_gene_names(self) -> List[str]:
        """List all gene names.

        Examples:
            >>> ensembl100.all_gene_names()[:3]
            ['A1BG', 'A1BG-AS1', 'A1CF']

        Returns:
            List[str]: Gene names
        """
        return self._uniquify_series(self.df[GENE_NAME])

    def all_protein_ids(self) -> List[str]:
        """List all protein IDs.

        Examples:
            >>> ensembl100.all_protein_ids()[:3]
            ['ENSP00000000233', 'ENSP00000000412', 'ENSP00000000442']

        Returns:
            List[str]: Protein IDs
        """
        return self._uniquify_series(self.df[PROTEIN_ID])

    def all_transcript_ids(self) -> List[str]:
        """List all transcript IDs.

        Examples:
            >>> ensembl100.all_transcript_ids()[:3]
            ['ENST00000000233', 'ENST00000000412', 'ENST00000000442']

        Returns:
            List[str]: Transcript IDs
        """
        return self._uniquify_series(self.df[TRANSCRIPT_ID])

    def all_transcript_names(self) -> List[str]:
        """List all transcript names.

        Examples:
            >>> ensembl100.all_transcript_names()[:3]
            ['A1BG-201', 'A1BG-202', 'A1BG-203']

        Returns:
            List[str]: Transcript names
        """
        return self._uniquify_series(self.df[TRANSCRIPT_NAME])

    def _uniquify_series(self, series: pd.Series) -> List:
        return sorted(series.dropna().unique().tolist())

    # ---------------------------------------------------------------------------------------------
    # <feature_symbol>s
    # ---------------------------------------------------------------------------------------------
    def contig_ids(self, feature: str) -> List[str]:
        """Given an ID or name, return the corresponding contig IDs.

        Examples:
            >>> ensembl100.contig_ids("BRCA2")
            ['13']

        Args:
            feature (str): Feature ID or name

        Returns:
            List[str]: Contig IDs
        """
        return self._query_feature(CONTIG_ID, feature)

    def exon_ids(self, feature: str) -> List[str]:
        """Given an ID or name, return the corresponding exon IDs.

        Examples:
            >>> ensembl100.exon_ids("BRCA2")[:3]
            ['ENSE00000939167', 'ENSE00000939168', 'ENSE00000939169']

        Args:
            feature (str): Feature ID or name

        Returns:
            List[str]: Exon IDs
        """
        return self._query_feature(EXON_ID, feature)

    def gene_ids(self, feature: str) -> List[str]:
        """Given an ID or name, return the corresponding gene IDs.

        Examples:
            >>> ensembl100.gene_ids("BRCA2")
            ['ENSG00000139618']

        Args:
            feature (str): Feature ID or name

        Returns:
            List[str]: Gene IDs
        """
        return self._query_feature(GENE_ID, feature)

    def gene_names(self, feature: str) -> List[str]:
        """Given an ID or name, return the corresponding gene names.

        Examples:
            >>> ensembl100.gene_names("ENSG00000139618")
            ['BRCA2']

        Args:
            feature (str): Feature ID or name

        Returns:
            List[str]: Gene names
        """
        return self._query_feature(GENE_NAME, feature)

    def protein_ids(self, feature: str) -> List[str]:
        """Given an ID or name, return the corresponding protein IDs.

        Examples:
            >>> ensembl100.protein_ids("BRCA2")[:3]
            ['ENSP00000369497', 'ENSP00000433168', 'ENSP00000434898']

        Args:
            feature (str): Feature ID or name

        Returns:
            List[str]: Protein IDs
        """
        return self._query_feature(PROTEIN_ID, feature)

    def transcript_ids(self, feature: str) -> List[str]:
        """Given an ID or name, return the corresponding transcript IDs.

        Examples:
            >>> ensembl100.transcript_ids("BRCA2")[:3]
            ['ENST00000380152', 'ENST00000470094', 'ENST00000528762']

        Args:
            feature (str): Feature ID or name

        Returns:
            List[str]: Transcript IDs
        """
        return self._query_feature(TRANSCRIPT_ID, feature)

    def transcript_names(self, feature: str) -> List[str]:
        """Given an ID or name, return the corresponding transcript names.

        Examples:
            >>> ensembl100.transcript_names("BRCA2")[:3]
            ['BRCA2-201', 'BRCA2-202', 'BRCA2-203']

        Args:
            feature (str): Feature ID or name

        Returns:
            List[str]: Transcript names
        """
        return self._query_feature(TRANSCRIPT_NAME, feature)

    def _query_feature(self, key: str, feature: str) -> List[str]:
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

    # ---------------------------------------------------------------------------------------------
    # normalize_id
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
    # is_<feature>
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
    # <feature>_sequence
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
    # <feature>_alias
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
        for key in [feature, strip_version(feature)]:
            if alias := alias_dict.get(key, []):
                return alias
        else:
            return []

    # ---------------------------------------------------------------------------------------------
    # is_canonical_transcript
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
        mask = (self.df[TRANSCRIPT_ID].isin(transcript_ids)) & (self.df["feature"] == "cdna")
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
        mask = (self.df[EXON_ID].isin(exon_ids)) & (self.df["feature"] == "exon")
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
    # <feature>_to_<feature>
    # ---------------------------------------------------------------------------------------------
    # TODO: add type hints
    def to_cdna(self, position) -> List:
        """Map a position to zero or more cDNA positions.

        Args:
            position (Position): Position or variant object.

        Returns:
            List[CdnaPosition]: If a variant was given, returns a variant, otherwise returns a position.
        """
        if position.is_fusion:
            fusion = cast(_Fusion, position)
            breakpoint1 = self.to_cdna(fusion.breakpoint1)
            breakpoint2 = self.to_cdna(fusion.breakpoint2)
            return [CdnaFusion(b1, b2) for b1, b2 in product(breakpoint1, breakpoint2)]
        elif position.is_small_variant:
            position = cast(_SmallVariant, position)
            if position.is_cdna:
                return self._cdna_to_cdna_variant(position)
            elif position.is_dna:
                return self._dna_to_cdna_variant(position)
            elif position.is_protein:
                return self._protein_to_cdna_variant(position)
            elif position.is_rna:
                return self._rna_to_cdna_variant(position)
        else:
            position = cast(_Position, position)
            if position.is_cdna:
                return self._cdna_to_cdna(position)
            elif position.is_dna:
                return self._dna_to_cdna(position)
            elif position.is_exon:
                return self._exon_to_cdna(position)
            elif position.is_protein:
                return self._protein_to_cdna(position)
            elif position.is_rna:
                return self._rna_to_cdna(position)

        raise AssertionError(f"Unknown position type for {position}")

    def to_dna(self, position) -> List:
        """Map a position to zero or more cDNA positions.

        Args:
            position (Position): Position or variant object.

        Returns:
            List[DnaPosition]: If a variant was given, returns a variant, otherwise returns a position.
        """
        if position.is_fusion:
            fusion = cast(_Fusion, position)
            breakpoint1 = self.to_dna(fusion.breakpoint1)
            breakpoint2 = self.to_dna(fusion.breakpoint2)
            return [DnaFusion(b1, b2) for b1, b2 in product(breakpoint1, breakpoint2)]
        elif position.is_small_variant:
            position = cast(_SmallVariant, position)
            if position.is_cdna:
                return self._cdna_to_dna_variant(position)
            elif position.is_dna:
                return self._dna_to_dna_variant(position)
            elif position.is_protein:
                return self._protein_to_dna_variant(position)
            elif position.is_rna:
                return self._rna_to_dna_variant(position)
        else:
            position = cast(_Position, position)
            if position.is_cdna:
                return self._cdna_to_dna(position)
            elif position.is_dna:
                return self._dna_to_dna(position)
            elif position.is_exon:
                return self._exon_to_dna(position)
            elif position.is_protein:
                return self._protein_to_dna(position)
            elif position.is_rna:
                return self._rna_to_dna(position)

        raise AssertionError(f"Unknown position type for {position}")

    def to_exon(self, position) -> List:
        """Map a position to zero or more cDNA positions.

        Args:
            position (Position): Position or variant object.

        Returns:
            List[ExonPosition]: If a variant was given, returns a variant, otherwise returns a position.
        """
        if position.is_fusion:
            fusion = cast(_Fusion, position)
            breakpoint1 = self.to_exon(fusion.breakpoint1)
            breakpoint2 = self.to_exon(fusion.breakpoint2)
            return [ExonFusion(b1, b2) for b1, b2 in product(breakpoint1, breakpoint2)]
        elif position.is_small_variant:
            position = cast(_SmallVariant, position)
            if position.is_cdna:
                return self._cdna_to_exon_variant(position)
            elif position.is_dna:
                return self._dna_to_exon_variant(position)
            elif position.is_protein:
                return self._protein_to_exon_variant(position)
            elif position.is_rna:
                return self._rna_to_exon_variant(position)
        else:
            position = cast(_Position, position)
            if position.is_cdna:
                return self._cdna_to_exon(position)
            elif position.is_dna:
                return self._dna_to_exon(position)
            elif position.is_exon:
                return self._exon_to_exon(position)
            elif position.is_protein:
                return self._protein_to_exon(position)
            elif position.is_rna:
                return self._rna_to_exon(position)

        raise AssertionError(f"Unknown position type for {position}")

    def to_protein(self, position) -> List:
        """Map a position to zero or more cDNA positions.

        Args:
            position (Position): Position or variant object.

        Returns:
            List[ProteinPosition]: If a variant was given, returns a variant, otherwise returns a position.
        """
        if position.is_fusion:
            fusion = cast(_Fusion, position)
            breakpoint1 = self.to_protein(fusion.breakpoint1)
            breakpoint2 = self.to_protein(fusion.breakpoint2)
            return [ProteinFusion(b1, b2) for b1, b2 in product(breakpoint1, breakpoint2)]
        elif position.is_small_variant:
            position = cast(_SmallVariant, position)
            if position.is_cdna:
                return self._cdna_to_protein_variant(position)
            elif position.is_dna:
                return self._dna_to_protein_variant(position)
            elif position.is_protein:
                return self._protein_to_protein_variant(position)
            elif position.is_rna:
                return self._rna_to_protein_variant(position)
        else:
            position = cast(_Position, position)
            if position.is_cdna:
                return self._cdna_to_protein(position)
            elif position.is_dna:
                return self._dna_to_protein(position)
            elif position.is_exon:
                return self._exon_to_protein(position)
            elif position.is_protein:
                return self._protein_to_protein(position)
            elif position.is_rna:
                return self._rna_to_protein(position)

        raise AssertionError(f"Unknown position type for {position}")

    def to_rna(self, position) -> List:
        """Map a position to zero or more cDNA positions.

        Args:
            position (Position): Position or variant object.

        Returns:
            List[RnaPosition]: If a variant was given, returns a variant, otherwise returns a position.
        """
        if position.is_fusion:
            fusion = cast(_Fusion, position)
            breakpoint1 = self.to_rna(fusion.breakpoint1)
            breakpoint2 = self.to_rna(fusion.breakpoint2)
            return [RnaFusion(b1, b2) for b1, b2 in product(breakpoint1, breakpoint2)]
        elif position.is_small_variant:
            position = cast(_SmallVariant, position)
            if position.is_cdna:
                return self._cdna_to_rna_variant(position)
            elif position.is_dna:
                return self._dna_to_rna_variant(position)
            elif position.is_protein:
                return self._protein_to_rna_variant(position)
            elif position.is_rna:
                return self._rna_to_rna_variant(position)
        else:
            position = cast(_Position, position)
            if position.is_cdna:
                return self._cdna_to_rna(position)
            elif position.is_dna:
                return self._dna_to_rna(position)
            elif position.is_exon:
                return self._exon_to_rna(position)
            elif position.is_protein:
                return self._protein_to_rna(position)
            elif position.is_rna:
                return self._rna_to_rna(position)

        raise AssertionError(f"Unknown position type for {position}")

    def _cdna_to_cdna(
        self, position: CdnaPosition, include_stop: bool = True
    ) -> List[CdnaPosition]:
        feature = ["CDS", "stop_codon"] if include_stop else ["CDS"]

        def convert(n: int, offset: int):
            result = []

            # If the offset position can be normalized to a non-offset position, do so. Otherwise
            # just return an offset position.
            if offset:
                for dna in self._cdna_to_dna(position, include_stop=include_stop):
                    for cdna in self._dna_to_cdna(dna):
                        if cdna.transcript_id == position.transcript_id:
                            result.append(cdna)

                if result:
                    return result

            mask = (
                (self.df[TRANSCRIPT_ID] == position.transcript_id)
                & (self.df["cdna_start"] <= n)
                & (self.df["cdna_end"] >= n)
                & (self.df["strand"] == position.strand)
                & (self.df["feature"].isin(feature))
            )
            for _, cds in self.df[mask].iterrows():
                result.append(
                    CdnaPosition(
                        contig_id=cds.contig_id,
                        start=position.start,
                        start_offset=offset,
                        end=position.end,
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

        result_start = convert(position.start, position.start_offset)
        if position.start == position.end:
            return sorted(result_start)
        else:
            result_end = convert(position.end, position.end_offset)
            return join_positions(result_start, result_end, TRANSCRIPT_ID)

    def _cdna_to_cdna_variant(
        self, position: _CdnaSmallVariant, include_stop: bool = True
    ) -> List[_CdnaSmallVariant]:
        result = []

        for cdna in self._cdna_to_cdna(position, include_stop=include_stop):
            result.extend(
                self._cdna_small_variant_from_cdna(cdna, position.refseq, position.altseq)
            )

        return result

    def _cdna_to_dna(self, position: CdnaPosition, include_stop: bool = True) -> List[DnaPosition]:
        feature = ["CDS", "stop_codon"] if include_stop else ["CDS"]

        def convert(n: int, offset: int):
            result = []

            mask = (
                (self.df[TRANSCRIPT_ID] == position.transcript_id)
                & (self.df["cdna_start"] <= n)
                & (self.df["cdna_end"] >= n)
                & (self.df["strand"] == position.strand)
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

        result_start = convert(position.start, position.start_offset)
        if position.start == position.end:
            return sorted(result_start)
        else:
            result_end = convert(position.end, position.end_offset)
            return join_positions(result_start, result_end, CONTIG_ID)

    def _cdna_to_dna_variant(
        self, position: _CdnaSmallVariant, include_stop: bool = True
    ) -> List[_DnaSmallVariant]:
        result = []

        for dna in self._cdna_to_dna(position, include_stop=include_stop):
            result.extend(self._dna_small_variant_from_dna(dna, position.refseq, position.altseq))

        return result

    def _cdna_to_exon(
        self, position: CdnaPosition, include_stop: bool = True
    ) -> List[ExonPosition]:
        feature = ["CDS", "stop_codon"] if include_stop else ["CDS"]

        def convert(n: int, offset: int):
            result = []

            # For an offset position, we need to calculate then equivalent DNA position then map
            # that to one or more exons. This is slower, so if there's no offset we can just map
            # directly to an exon.
            if offset:
                for dna in self._cdna_to_dna(position, include_stop=include_stop):
                    for exon in self._dna_to_exon(dna):
                        if exon.transcript_id == position.transcript_id:
                            result.append(exon)

                return result

            mask_cds = (
                (self.df[TRANSCRIPT_ID] == position.transcript_id)
                & (self.df["cdna_start"] <= n)
                & (self.df["cdna_end"] >= n)
                & (self.df["strand"] == position.strand)
                & (self.df["feature"].isin(feature))
            )
            for _, cds in self.df[mask_cds].iterrows():
                mask_exon = (
                    (self.df[TRANSCRIPT_ID] == position.transcript_id)
                    & (self.df["exon_number"] == cds.exon_number)
                    & (self.df["feature"] == "exon")
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

        result_start = convert(position.start, position.start_offset)
        if position.start == position.end:
            return sorted(result_start)
        else:
            result_end = convert(position.end, position.end_offset)
            return join_positions(result_start, result_end, TRANSCRIPT_ID)

    def _cdna_to_exon_variant(
        self, position: _CdnaSmallVariant, include_stop: bool = True
    ) -> List[_ExonSmallVariant]:
        result = []

        for exon in self._cdna_to_exon(position, include_stop=include_stop):
            result.extend(
                self._exon_small_variant_from_exon(exon, position.refseq, position.altseq)
            )

        return result

    def _cdna_to_protein(self, position: CdnaPosition) -> List[ProteinPosition]:
        result = []

        for cdna in self._cdna_to_cdna(position, include_stop=False):
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

    def _cdna_to_protein_variant(self, position: _CdnaSmallVariant) -> List[_ProteinSmallVariant]:
        # TODO: A lot of this function is duplicated from _cdna_to_protein()
        result = []

        for cdna in self._cdna_to_cdna(position, include_stop=False):
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
            result.extend(
                self._protein_small_variant(cdna, protein, position.refseq, position.altseq)
            )

        return result

    def _cdna_to_rna(self, position: CdnaPosition, include_stop: bool = True) -> List[RnaPosition]:
        feature = ["CDS", "stop_codon"] if include_stop else ["CDS"]

        def convert(n: int, offset: int):
            result = []

            # If the offset position can be normalized to a non-offset position, do so. Otherwise
            # just return an offset position.
            if offset:
                for dna in self._cdna_to_dna(position, include_stop=include_stop):
                    for rna in self._dna_to_rna(dna):
                        if rna.transcript_id == position.transcript_id:
                            result.append(rna)

                if result:
                    return result

            mask = (
                (self.df[TRANSCRIPT_ID] == position.transcript_id)
                & (self.df["cdna_start"] <= n)
                & (self.df["cdna_end"] >= n)
                & (self.df["strand"] == position.strand)
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

        result_start = convert(position.start, position.start_offset)
        if position.start == position.end:
            return sorted(result_start)
        else:
            result_end = convert(position.end, position.end_offset)
            return join_positions(result_start, result_end, TRANSCRIPT_ID)

    def _cdna_to_rna_variant(
        self, position: _CdnaSmallVariant, include_stop: bool = True
    ) -> List[_RnaSmallVariant]:
        result = []

        for rna in self._cdna_to_rna(position, include_stop=include_stop):
            result.extend(self._rna_small_variant_from_rna(rna, position.refseq, position.altseq))

        return result

    def _dna_to_cdna(self, position: DnaPosition, include_stop: bool = True) -> List[CdnaPosition]:
        feature = ["CDS", "stop_codon"] if include_stop else ["CDS"]

        def convert(n: int, offset: int):
            result = []

            if position.strand == "-":
                n_ = n - offset
            else:
                n_ = n + offset

            mask = (
                (self.df[CONTIG_ID] == position.contig_id)
                & (self.df["start"] <= n_)
                & (self.df["end"] >= n_)
                & (self.df["strand"] == position.strand)
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

        result_start = convert(position.start, position.start_offset)
        if position.start == position.end:
            return sorted(result_start)
        else:
            result_end = convert(position.end, position.end_offset)
            return join_positions(result_start, result_end, TRANSCRIPT_ID)

    def _dna_to_cdna_variant(
        self, position: _DnaSmallVariant, include_stop: bool = True
    ) -> List[_CdnaSmallVariant]:
        result = []

        for cdna in self._dna_to_cdna(position, include_stop=include_stop):
            result.extend(
                self._cdna_small_variant_from_cdna(cdna, position.refseq, position.altseq)
            )

        return result

    def _dna_to_dna(self, position: DnaPosition) -> List[DnaPosition]:
        result = []

        if position.strand == "-":
            new_start = position.start - position.start_offset
            new_end = position.end - position.end_offset
        else:
            new_start = position.start + position.start_offset
            new_end = position.end + position.end_offset

        # Sort the start and end positions after adjusting by offsets
        new_start, new_end = sorted([new_start, new_end])

        # TODO: Check that new new_start is actually on the contig
        result.append(
            DnaPosition(
                contig_id=position.contig_id,
                start=new_start,
                start_offset=0,
                end=new_end,
                end_offset=0,
                strand=position.strand,
            )
        )

        return result

    def _dna_to_dna_variant(self, position: _DnaSmallVariant) -> List[_DnaSmallVariant]:
        result = []

        for dna in self._dna_to_dna(position):
            result.extend(self._dna_small_variant_from_dna(dna, position.refseq, position.altseq))

        return result

    def _dna_to_exon(self, position: DnaPosition) -> List[ExonPosition]:
        def convert(n: int, offset: int):
            result = []

            if position.strand == "-":
                n_ = n - offset
            else:
                n_ = n + offset

            mask = (
                (self.df[CONTIG_ID] == position.contig_id)
                & (self.df["start"] <= n_)
                & (self.df["end"] >= n_)
                & (self.df["strand"] == position.strand)
                & (self.df["feature"] == "exon")
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

        result_start = convert(position.start, position.start_offset)
        if position.start == position.end:
            return sorted(result_start)
        else:
            result_end = convert(position.end, position.end_offset)
            return join_positions(result_start, result_end, TRANSCRIPT_ID)

    def _dna_to_exon_variant(self, position: _DnaSmallVariant) -> List[_ExonSmallVariant]:
        result = []

        for exon in self._dna_to_exon(position):
            result.extend(
                self._exon_small_variant_from_exon(exon, position.refseq, position.altseq)
            )

        return result

    def _dna_to_protein(self, position: DnaPosition) -> List[ProteinPosition]:
        result = []

        for cdna in self._dna_to_cdna(position, include_stop=False):
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

    def _dna_to_protein_variant(self, position: _DnaSmallVariant) -> List[_ProteinSmallVariant]:
        # TODO: A lot of this function is duplicated from _dna_to_protein()
        result = []

        for cdna in self._dna_to_cdna(position, include_stop=False):
            # Offset cDNA position are assumed to not map to a protein
            if cdna.start_offset or cdna.end_offset:
                continue

            pstart = calc_cdna_to_protein(cdna.start)
            pend = calc_cdna_to_protein(cdna.end)
            protein = ProteinPosition.copy_from(
                cdna, start=pstart, start_offset=0, end=pend, end_offset=0
            )
            result.extend(
                self._protein_small_variant(cdna, protein, position.refseq, position.altseq)
            )

        return result

    def _dna_to_rna(self, position: DnaPosition) -> List[RnaPosition]:
        def convert(n: int, offset: int):
            result = []

            if position.strand == "-":
                n_ = n - offset
            else:
                n_ = n + offset

            mask = (
                (self.df[CONTIG_ID] == position.contig_id)
                & (self.df["start"] <= n_)
                & (self.df["end"] >= n_)
                & (self.df["strand"] == position.strand)
                & (self.df["feature"] == "exon")
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

        result_start = convert(position.start, position.start_offset)
        if position.start == position.end:
            return sorted(result_start)
        else:
            result_end = convert(position.end, position.end_offset)
            return join_positions(result_start, result_end, TRANSCRIPT_ID)

    def _dna_to_rna_variant(self, position: _DnaSmallVariant) -> List[_RnaSmallVariant]:
        result = []

        for rna in self._dna_to_rna(position):
            result.extend(self._rna_small_variant_from_rna(rna, position.refseq, position.altseq))

        return result

    def _exon_to_cdna(
        self, position: ExonPosition, include_stop: bool = True
    ) -> List[CdnaPosition]:
        feature = ["CDS", "stop_codon"] if include_stop else ["CDS"]

        def convert(n: int, offset: int):
            result = []

            # TODO: Is there a reasonable case where an exon position would have an offset?
            assert not offset, offset

            mask = (
                (self.df[TRANSCRIPT_ID] == position.transcript_id)
                & (self.df["exon_number"] == float(n))
                & (self.df["strand"] == position.strand)
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

        result_start = convert(position.start, position.start_offset)
        if position.start == position.end:
            return sorted(result_start)
        else:
            result_end = convert(position.end, position.end_offset)
            return join_positions(result_start, result_end, TRANSCRIPT_ID)

    def _exon_to_dna(self, position: ExonPosition) -> List[DnaPosition]:
        def convert(n: int, offset: int):
            result = []

            # TODO: Is there a reasonable case where an exon position would have an offset?
            assert not offset, offset

            mask = (
                (self.df[TRANSCRIPT_ID] == position.transcript_id)
                & (self.df["exon_number"] == float(n))
                & (self.df["strand"] == position.strand)
                & (self.df["feature"] == "exon")
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

        result_start = convert(position.start, position.start_offset)
        if position.start == position.end:
            return sorted(result_start)
        else:
            result_end = convert(position.end, position.end_offset)
            return join_positions(result_start, result_end, CONTIG_ID)

    def _exon_to_exon(self, position: ExonPosition) -> List[ExonPosition]:
        def convert(n: int, offset):
            result = []

            # TODO: Is there a reasonable case where an exon position would have an offset?
            assert not offset, offset

            mask = (
                (self.df[TRANSCRIPT_ID] == position.transcript_id)
                & (self.df["exon_number"] == float(n))
                & (self.df["strand"] == position.strand)
                & (self.df["feature"] == "exon")
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

        result_start = convert(position.start, position.start_offset)
        if position.start == position.end:
            return sorted(result_start)
        else:
            result_end = convert(position.end, position.end_offset)
            return join_positions(result_start, result_end, TRANSCRIPT_ID)

    def _exon_to_protein(self, position: ExonPosition) -> List[ProteinPosition]:
        result = []

        for cdna in self._exon_to_cdna(position, include_stop=False):
            result.extend(self._cdna_to_protein(cdna))

        return result

    def _exon_to_rna(self, position: ExonPosition) -> List[RnaPosition]:
        def convert(n: int, offset: int):
            result = []

            # TODO: Is there a reasonable case where an exon position would have an offset?
            assert not offset, offset

            mask = (
                (self.df[TRANSCRIPT_ID] == position.transcript_id)
                & (self.df["exon_number"] == float(n))
                & (self.df["strand"] == position.strand)
                & (self.df["feature"] == "exon")
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

        result_start = convert(position.start, position.start_offset)
        if position.start == position.end:
            return sorted(result_start)
        else:
            result_end = convert(position.end, position.end_offset)
            return join_positions(result_start, result_end, TRANSCRIPT_ID)

    def _protein_to_cdna(self, position: ProteinPosition) -> List[CdnaPosition]:
        def convert(n: int):
            return ((n - 1) * 3) + 1

        # TODO: Is there a reasonable case where an protein position would have an offset?
        assert not position.start_offset, position.start_offset
        assert not position.end_offset, position.end_offset

        cdna_start = convert(position.start)
        cdna_end = convert(position.end) + 2

        return [CdnaPosition.copy_from(position, start=cdna_start, end=cdna_end)]

    def _protein_to_cdna_variant(self, position: _ProteinSmallVariant) -> List[_CdnaSmallVariant]:
        result = []

        for cdna in self._protein_to_cdna(position):
            result.extend(
                self._cdna_small_variant_from_protein(cdna, position.refseq, position.altseq)
            )

        return result

    def _protein_to_dna(self, position: ProteinPosition) -> List[DnaPosition]:
        result = []

        for cdna in self._protein_to_cdna(position):
            result.extend(self._cdna_to_dna(cdna))

        return sorted(set(result))

    def _protein_to_dna_variant(self, position: _ProteinSmallVariant) -> List[_DnaSmallVariant]:
        result = []

        for cdna in self._protein_to_cdna_variant(position):
            result.extend(self._cdna_to_dna_variant(cdna, include_stop=True))

        return sorted(set(result))

    def _protein_to_exon(self, position: ProteinPosition) -> List[ExonPosition]:
        result = []

        for cdna in self._protein_to_cdna(position):
            result.extend(self._cdna_to_exon(cdna))

        return sorted(set(result))

    def _protein_to_exon_variant(self, position: _ProteinSmallVariant) -> List[_ExonSmallVariant]:
        result = []

        for cdna in self._protein_to_cdna_variant(position):
            result.extend(self._cdna_to_exon_variant(cdna, include_stop=True))

        return sorted(set(result))

    def _protein_to_protein(self, position: ProteinPosition) -> List[ProteinPosition]:
        result = []

        for cdna in self._protein_to_cdna(position):
            result.extend(self._cdna_to_protein(cdna))

        return sorted(set(result))

    def _protein_to_protein_variant(
        self, position: _ProteinSmallVariant
    ) -> List[_ProteinSmallVariant]:
        result = []

        for cdna in self._protein_to_cdna_variant(position):
            result.extend(self._cdna_to_protein_variant(cdna))

        return sorted(set(result))

    def _protein_to_rna(self, position: ProteinPosition) -> List[RnaPosition]:
        result = []

        for cdna in self._protein_to_cdna(position):
            result.extend(self._cdna_to_rna(cdna))

        return sorted(set(result))

    def _protein_to_rna_variant(self, position: _ProteinSmallVariant) -> List[_RnaSmallVariant]:
        result = []

        for cdna in self._protein_to_cdna_variant(position):
            result.extend(self._cdna_to_rna_variant(cdna))

        return sorted(set(result))

    def _rna_to_cdna(self, position: RnaPosition, include_stop: bool = True) -> List[CdnaPosition]:
        feature = ["CDS", "stop_codon"] if include_stop else ["CDS"]

        def convert(n: int, offset: int):
            result = []

            # If the offset position can be normalized to a non-offset position, do so. Otherwise
            # just return an offset position.
            if offset:
                for dna in self._rna_to_dna(position):
                    for cdna in self._dna_to_cdna(dna):
                        if cdna.transcript_id == position.transcript_id:
                            result.append(cdna)

                if result:
                    return result

            mask = (
                (self.df[TRANSCRIPT_ID] == position.transcript_id)
                & (self.df["transcript_start"] <= n)
                & (self.df["transcript_end"] >= n)
                & (self.df["strand"] == position.strand)
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

        result_start = convert(position.start, position.start_offset)
        if position.start == position.end:
            return sorted(result_start)
        else:
            result_end = convert(position.end, position.end_offset)
            return join_positions(result_start, result_end, TRANSCRIPT_ID)

    def _rna_to_cdna_variant(
        self, position: _RnaSmallVariant, include_stop: bool = True
    ) -> List[_CdnaSmallVariant]:
        result = []

        for cdna in self._rna_to_cdna(position, include_stop=include_stop):
            result.extend(
                self._cdna_small_variant_from_cdna(cdna, position.refseq, position.altseq)
            )

        return result

    def _rna_to_dna(self, position: RnaPosition) -> List[DnaPosition]:
        def convert(n: int, offset: int):
            result = []

            mask = (
                (self.df[TRANSCRIPT_ID] == position.transcript_id)
                & (self.df["transcript_start"] <= n)
                & (self.df["transcript_end"] >= n)
                & (self.df["strand"] == position.strand)
                & (self.df["feature"] == "exon")
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

        result_start = convert(position.start, position.start_offset)
        if position.start == position.end:
            return sorted(result_start)
        else:
            result_end = convert(position.end, position.end_offset)
            return join_positions(result_start, result_end, CONTIG_ID)

    def _rna_to_dna_variant(self, position: _RnaSmallVariant) -> List[_DnaSmallVariant]:
        result = []

        for dna in self._rna_to_dna(position):
            result.extend(self._dna_small_variant_from_dna(dna, position.refseq, position.altseq))

        return result

    def _rna_to_exon(self, position: RnaPosition) -> List[ExonPosition]:
        def convert(n: int, offset: int):
            result = []

            # For an offset position, we need to calculate then equivalent DNA position then map
            # that to one or more exons. This is slower, so if there's no offset we can just map
            # directly to an exon.
            if offset:
                for dna in self._rna_to_dna(position):
                    for exon in self._dna_to_exon(dna):
                        if exon.transcript_id == position.transcript_id:
                            result.append(exon)

                return result

            mask = (
                (self.df[TRANSCRIPT_ID] == position.transcript_id)
                & (self.df["transcript_start"] <= n)
                & (self.df["transcript_end"] >= n)
                & (self.df["strand"] == position.strand)
                & (self.df["feature"] == "exon")
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

        result_start = convert(position.start, position.start_offset)
        if position.start == position.end:
            return sorted(result_start)
        else:
            result_end = convert(position.end, position.end_offset)
            return join_positions(result_start, result_end, TRANSCRIPT_ID)

    def _rna_to_exon_variant(self, position: _RnaSmallVariant) -> List[_ExonSmallVariant]:
        result = []

        for exon in self._rna_to_exon(position):
            result.extend(
                self._exon_small_variant_from_exon(exon, position.refseq, position.altseq)
            )

        return result

    def _rna_to_protein(self, position: RnaPosition) -> List[ProteinPosition]:
        result = []

        for cdna in self._rna_to_cdna(position, include_stop=False):
            result.extend(self._cdna_to_protein(cdna))

        return sorted(set(result))

    def _rna_to_protein_variant(self, position: _RnaSmallVariant) -> List[_ProteinSmallVariant]:
        result = []

        for cdna in self._rna_to_cdna_variant(position, include_stop=False):
            result.extend(self._cdna_to_protein_variant(cdna))

        return sorted(set(result))

    def _rna_to_rna(self, position: RnaPosition) -> List[RnaPosition]:
        def convert(n: int, offset: int):
            result = []

            # If the offset position can be normalized to a non-offset position, do so. Otherwise
            # just return an offset position.
            if offset:
                for dna in self._rna_to_dna(position):
                    for rna in self._dna_to_rna(dna):
                        if rna.transcript_id == position.transcript_id:
                            result.append(rna)

                return result

            mask = (
                (self.df[TRANSCRIPT_ID] == position.transcript_id)
                & (self.df["transcript_start"] <= n)
                & (self.df["transcript_end"] >= n)
                & (self.df["strand"] == position.strand)
                & (self.df["feature"] == "exon")
            )
            for _, exon in self.df[mask].iterrows():
                result.append(
                    RnaPosition(
                        contig_id=exon.contig_id,
                        start=position.start,
                        start_offset=offset,
                        end=position.end,
                        end_offset=offset,
                        strand=exon.strand,
                        gene_id=exon.gene_id,
                        gene_name=exon.gene_name,
                        transcript_id=exon.transcript_id,
                        transcript_name=exon.transcript_name,
                    )
                )

            return result

        result_start = convert(position.start, position.start_offset)
        if position.start == position.end:
            return sorted(result_start)
        else:
            result_end = convert(position.end, position.end_offset)
            return join_positions(result_start, result_end, TRANSCRIPT_ID)

    def _rna_to_rna_variant(self, position: _RnaSmallVariant) -> List[_RnaSmallVariant]:
        result = []

        for rna in self._rna_to_rna(position):
            result.extend(self._rna_small_variant_from_rna(rna, position.refseq, position.altseq))

        return result

    def _cdna_small_variant_from_cdna(
        self, cdna: CdnaPosition, refseq: str, altseq: str
    ) -> List[_CdnaSmallVariant]:
        """Convert a cDNA position plus ref/alt nucleotides into a cDNA variant.

        Args:
            cdna (CdnaPosition): _Variant position
            refseq (str): reference allele
            altseq (str): alternate allele

        Raises:
            NotImplementedError: An unsupported combination of reference/alternate alleles was given
            ValueError: The given reference allele does not match the annotated reference allele

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

    def _cdna_small_variant_from_protein(
        self, cdna: CdnaPosition, refaa: str, altaa: str
    ) -> List[_CdnaSmallVariant]:
        """Convert a cDNA position plus ref/alt amino acids into a cDNA variant.

        Args:
            cdna (CdnaPosition): cDNA position
            refseq (str): Reference allele
            altaa (str): Alternate allele

        Raises:
            NotImplementedError: An unsupported combination of reference/alternate alleles was given

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
            NotImplementedError: An unsupported combination of reference/alternate alleles was given

        Returns:
            List[_DnaSmallVariant]: One or more DNA variants
        """
        variant_list = []

        # TODO: DNA sequence get is slow
        # ref_annotated = self.sequence(dna)
        for ref, alt in product(expand_nt(refseq), expand_nt(altseq)):
            # # Assert that the given ref matches the annotated one
            # if ref != ref_annotated:
            #     raise ValueError(
            #         f"Given ref allele '{ref}' does not match annotated ref allele '{ref_annotated}' for {dna}"
            #     )

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
            NotImplementedError: An unsupported combination of reference/alternate alleles was given

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
            NotImplementedError: An unsupported combination of reference/alternate alleles was given

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
            NotImplementedError: An unsupported combination of reference/alternate alleles was given
            ValueError: The given reference allele does not match the annotated reference allele

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


def join_positions(start: List[Position], end: List[Position], merge_on: str) -> List[Position]:
    """Return the combination of two list of position or variant objects - one of start positions
    and one of end positions - into one list by the given key (e.g. 'transcript_id').

    All positions must be of the same class.

    Args:
        start (List[PositionOrSmallVariantType]): Start positions
        end (List[PositionOrSmallVariantType]): End positions
        merge_on (str): Attribute to merge on (e.g. 'transcript_id')

    Returns:
        List[PositionOrSmallVariantType]: New position or variant objects
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
