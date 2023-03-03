from __future__ import annotations

from itertools import product
from typing import Callable, Dict, List, Optional, Tuple, Union

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
    CdnaPosition,
    DnaPosition,
    ExonPosition,
    MappablePositionOrSmallVariant,
    ProteinPosition,
    RnaPosition,
    _CdnaSmallVariant,
    _DnaSmallVariant,
    _ExonSmallVariant,
    _ProteinSmallVariant,
    _RnaSmallVariant,
)
from .tables import AMINO_ACID_TABLE
from .utils import (
    calc_cdna_to_protein,
    expand_nt,
    reverse_complement,
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
    def cds_sequence(
        self, transcript_id: str, start: Optional[int] = None, end: Optional[int] = None
    ) -> str:
        """Return the CDS sequence between the given position(s), inclusive.

        Examples:
            >>> ensembl100.cds_sequence("ENST00000380152", 10, end=20)
            'GGATCCAAAGA'

        Args:
            transcript_id (str): transcript ID
            start (Optional[int], optional): Start position of the sequence to return. Defaults to None.
            end (Optional[int], optional): End position of the sequence to return. Defaults to the same as `start`.

        Returns:
            str: The CDS sequence
        """
        return self._sequence(self.cds_fasta, transcript_id, start=start, end=end)

    def dna_sequence(
        self,
        contig_id: str,
        start: Optional[int] = None,
        end: Optional[int] = None,
        strand: str = "+",
    ) -> str:
        """Return the DNA sequence between the given position(s), inclusive.

        Examples:
            >>> ensembl100.dna_sequence("13", 32315086, end=32315096)
            'AAGCTTTTGTA'

        Args:
            contig_id (str): contig ID
            start (Optional[int], optional): Start position of the sequence to return. Defaults to None.
            end (Optional[int], optional): End position of the sequence to return. Defaults to the same as `start`.

        Returns:
            str: The DNA sequence
        """
        return self._sequence(self.dna_fasta, contig_id, start=start, end=end, strand=strand)

    def protein_sequence(
        self, protein_id: str, start: Optional[int] = None, end: Optional[int] = None
    ) -> str:
        """Return the protein sequence between the given position(s), inclusive.

        Examples:
            >>> ensembl100.protein_sequence("ENSP00000369497", 1, end=3)
            'MPI'

        Args:
            protein_id (str): protein ID
            start (Optional[int], optional): Start position of the sequence to return. Defaults to None.
            end (Optional[int], optional): End position of the sequence to return. Defaults to the same as `start`.

        Returns:
            str: The protein sequence
        """
        return self._sequence(self.protein_fasta, protein_id, start=start, end=end)

    def rna_sequence(
        self, transcript_id: str, start: Optional[int] = None, end: Optional[int] = None
    ) -> str:
        """Return the RNA sequence between the given position(s), inclusive.

        Examples:
            >>> ensembl100.rna_sequence("ENST00000380152", 1, end=10)
            'GGGCTTGTGG'

        Args:
            transcript_id (str): transcript ID
            start (Optional[int], optional): Start position of the sequence to return. Defaults to None.
            end (Optional[int], optional): End position of the sequence to return. Defaults to the same as `start`.

        Returns:
            str: The RNA sequence
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
            canonical (bool, optional): Only map to the canonical transcript. Defaults to False.

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
                    _data=self,
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
            canonical (bool, optional): Only map to the canonical transcript. Defaults to False.

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
                        _data=self,
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
            canonical (bool, optional): Only map to the canonical transcript. Defaults to False.

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
                    _data=self,
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
            canonical (bool, optional): Only map to the canonical transcript. Defaults to False.

        Returns:
            List[DnaPosition]: One or more DNA positions.
        """
        result = []

        gene_ids = self.gene_ids(feature)
        mask = (self.df[GENE_ID].isin(gene_ids)) & (self.df["feature"] == "gene")
        for _, gene in self.df[mask].iterrows():
            result.append(
                DnaPosition(
                    _data=self,
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
            canonical (bool, optional): Only map to the canonical transcript. Defaults to False.

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
            canonical (bool, optional): Only map to the canonical transcript. Defaults to False.

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
                    _data=self,
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
    def cdna_to_cdna(
        self,
        feature: str,
        start: int,
        end: Optional[int] = None,
        strand: Optional[str] = None,
        start_offset: Optional[int] = None,
        end_offset: Optional[int] = None,
        refseq: str = "",
        altseq: str = "",
        canonical: bool = False,
    ) -> Union[List[CdnaPosition], List[_CdnaSmallVariant]]:
        """Map a cDNA position to zero or more cDNA positions.

        Args:
            feature (str): Feature ID or name
            start (int): Start position
            end (Optional[int], optional): End position. Defaults to the same as `start`.
            strand (Optional[str], optional): Strand ('+' or '-'). Defaults to checking either.
            start_offset (Optional[int], optional): Nucleotide offset from `start`. Defaults to 0.
            end_offset (Optional[int], optional): Nucleotide offset from `end`. Defaults to 0.
            refseq (str, optional): Reference allele. Defaults to "".
            altseq (str, optional): Alternate allele. Defaults to "".
            canonical (bool, optional): Only map to the canonical transcript. Defaults to False.

        Returns:
            Union[List[CdnaPosition], List[_CdnaSmallVariant]]:
                If `refseq` or `altseq` was given, returns a variant. Otherwise returns a position.
        """
        if refseq or altseq:
            return self._map_variant(
                self.transcript_ids,
                self._cdna_to_cdna_variant,
                feature,
                start,
                refseq,
                altseq,
                start_offset,
                end,
                end_offset,
                strand,
                canonical,
            )
        else:
            return self._map(
                self.transcript_ids,
                self._cdna_to_cdna,
                feature,
                start,
                start_offset,
                end,
                end_offset,
                strand,
                canonical,
            )

    def cdna_to_dna(
        self,
        feature: str,
        start: int,
        end: Optional[int] = None,
        strand: Optional[str] = None,
        start_offset: Optional[int] = None,
        end_offset: Optional[int] = None,
        refseq: str = "",
        altseq: str = "",
        canonical: bool = False,
    ) -> Union[List[DnaPosition], List[_DnaSmallVariant]]:
        """Map a cDNA position to zero or more DNA positions.

        Args:
            feature (str): Feature ID or name
            start (int): Start position
            end (Optional[int], optional): End position. Defaults to the same as `start`.
            strand (Optional[str], optional): Strand ('+' or '-'). Defaults to checking either.
            start_offset (Optional[int], optional): Nucleotide offset from `start`. Defaults to 0.
            end_offset (Optional[int], optional): Nucleotide offset from `end`. Defaults to 0.
            refseq (str, optional): Reference allele. Defaults to "".
            altseq (str, optional): Alternate allele. Defaults to "".
            canonical (bool, optional): Only map to the canonical transcript. Defaults to False.

        Returns:
            Union[List[DnaPosition], List[_DnaSmallVariant]]:
                If `refseq` or `altseq` was given, returns a variant. Otherwise returns a position.
        """
        if refseq or altseq:
            return self._map_variant(
                self.transcript_ids,
                self._cdna_to_dna_variant,
                feature,
                start,
                refseq,
                altseq,
                start_offset,
                end,
                end_offset,
                strand,
                canonical,
            )
        else:
            return self._map(
                self.transcript_ids,
                self._cdna_to_dna,
                feature,
                start,
                start_offset,
                end,
                end_offset,
                strand,
                canonical,
            )

    def cdna_to_exon(
        self,
        feature: str,
        start: int,
        end: Optional[int] = None,
        strand: Optional[str] = None,
        start_offset: Optional[int] = None,
        end_offset: Optional[int] = None,
        refseq: str = "",
        altseq: str = "",
        canonical: bool = False,
    ) -> Union[List[ExonPosition], List[_ExonSmallVariant]]:
        """Map a cDNA position to zero or more exon positions.

        Args:
            feature (str): Feature ID or name
            start (int): Start position
            end (Optional[int], optional): End position. Defaults to the same as `start`.
            strand (Optional[str], optional): Strand ('+' or '-'). Defaults to checking either.
            start_offset (Optional[int], optional): Nucleotide offset from `start`. Defaults to 0.
            end_offset (Optional[int], optional): Nucleotide offset from `end`. Defaults to 0.
            refseq (str, optional): Reference allele. Defaults to "".
            altseq (str, optional): Alternate allele. Defaults to "".
            canonical (bool, optional): Only map to the canonical transcript. Defaults to False.

        Returns:
            Union[List[ExonPosition], List[_ExonSmallVariant]]:
                If `refseq` or `altseq` was given, returns a variant. Otherwise returns a position.
        """
        if refseq or altseq:
            return self._map_variant(
                self.transcript_ids,
                self._cdna_to_exon_variant,
                feature,
                start,
                refseq,
                altseq,
                start_offset,
                end,
                end_offset,
                strand,
                canonical,
            )
        else:
            return self._map(
                self.transcript_ids,
                self._cdna_to_exon,
                feature,
                start,
                start_offset,
                end,
                end_offset,
                strand,
                canonical,
            )

    def cdna_to_protein(
        self,
        feature: str,
        start: int,
        end: Optional[int] = None,
        strand: Optional[str] = None,
        start_offset: Optional[int] = None,
        end_offset: Optional[int] = None,
        refseq: str = "",
        altseq: str = "",
        canonical: bool = False,
    ) -> Union[List[ProteinPosition], List[_ProteinSmallVariant]]:
        """Map a cDNA position to zero or more protein positions.

        Args:
            feature (str): Feature ID or name
            start (int): Start position
            end (Optional[int], optional): End position. Defaults to the same as `start`.
            strand (Optional[str], optional): Strand ('+' or '-'). Defaults to checking either.
            start_offset (Optional[int], optional): Nucleotide offset from `start`. Defaults to 0.
            end_offset (Optional[int], optional): Nucleotide offset from `end`. Defaults to 0.
            refseq (str, optional): Reference allele. Defaults to "".
            altseq (str, optional): Alternate allele. Defaults to "".
            canonical (bool, optional): Only map to the canonical transcript. Defaults to False.

        Returns:
            Union[List[ProteinPosition], List[_ProteinSmallVariant]]:
                If `refseq` or `altseq` was given, returns a variant. Otherwise returns a position.
        """
        if refseq or altseq:
            return self._map_variant(
                self.transcript_ids,
                self._cdna_to_protein_variant,
                feature,
                start,
                refseq,
                altseq,
                start_offset,
                end,
                end_offset,
                strand,
                canonical,
            )
        else:
            return self._map(
                self.transcript_ids,
                self._cdna_to_protein,
                feature,
                start,
                start_offset,
                end,
                end_offset,
                strand,
                canonical,
            )

    def cdna_to_rna(
        self,
        feature: str,
        start: int,
        end: Optional[int] = None,
        strand: Optional[str] = None,
        start_offset: Optional[int] = None,
        end_offset: Optional[int] = None,
        refseq: str = "",
        altseq: str = "",
        canonical: bool = False,
    ) -> Union[List[RnaPosition], List[_RnaSmallVariant]]:
        """Map a cDNA position to zero or more RNA positions.

        Args:
            feature (str): Feature ID or name
            start (int): Start position
            end (Optional[int], optional): End position. Defaults to the same as `start`.
            strand (Optional[str], optional): Strand ('+' or '-'). Defaults to checking either.
            start_offset (Optional[int], optional): Nucleotide offset from `start`. Defaults to 0.
            end_offset (Optional[int], optional): Nucleotide offset from `end`. Defaults to 0.
            refseq (str, optional): Reference allele. Defaults to "".
            altseq (str, optional): Alternate allele. Defaults to "".
            canonical (bool, optional): Only map to the canonical transcript. Defaults to False.

        Returns:
            Union[List[RnaPosition], List[_RnaSmallVariant]]:
                If `refseq` or `altseq` was given, returns a variant. Otherwise returns a position.
        """
        if refseq or altseq:
            return self._map_variant(
                self.transcript_ids,
                self._cdna_to_rna_variant,
                feature,
                start,
                refseq,
                altseq,
                start_offset,
                end,
                end_offset,
                strand,
                canonical,
            )
        else:
            return self._map(
                self.transcript_ids,
                self._cdna_to_rna,
                feature,
                start,
                start_offset,
                end,
                end_offset,
                strand,
                canonical,
            )

    def dna_to_cdna(
        self,
        feature: str,
        start: int,
        end: Optional[int] = None,
        strand: Optional[str] = None,
        start_offset: Optional[int] = None,
        end_offset: Optional[int] = None,
        refseq: str = "",
        altseq: str = "",
        canonical: bool = False,
    ) -> Union[List[CdnaPosition], List[_CdnaSmallVariant]]:
        """Map a DNA position to zero or more cDNA positions.

        Args:
            feature (str): Feature ID or name
            start (int): Start position
            end (Optional[int], optional): End position. Defaults to the same as `start`.
            strand (Optional[str], optional): Strand ('+' or '-'). Defaults to checking either.
            start_offset (Optional[int], optional): Nucleotide offset from `start`. Defaults to 0.
            end_offset (Optional[int], optional): Nucleotide offset from `end`. Defaults to 0.
            refseq (str, optional): Reference allele. Defaults to "".
            altseq (str, optional): Alternate allele. Defaults to "".
            canonical (bool, optional): Only map to the canonical transcript. Defaults to False.

        Returns:
            Union[List[CdnaPosition], List[_CdnaSmallVariant]]:
                If `refseq` or `altseq` was given, returns a variant. Otherwise returns a position.
        """
        if refseq or altseq:
            return self._map_variant(
                self.contig_ids,
                self._dna_to_cdna_variant,
                feature,
                start,
                refseq,
                altseq,
                start_offset,
                end,
                end_offset,
                strand,
                canonical,
            )
        else:
            return self._map(
                self.contig_ids,
                self._dna_to_cdna,
                feature,
                start,
                start_offset,
                end,
                end_offset,
                strand,
                canonical,
            )

    def dna_to_dna(
        self,
        feature: str,
        start: int,
        end: Optional[int] = None,
        strand: Optional[str] = None,
        start_offset: Optional[int] = None,
        end_offset: Optional[int] = None,
        refseq: str = "",
        altseq: str = "",
        canonical: bool = False,
    ) -> Union[List[DnaPosition], List[_DnaSmallVariant]]:
        """Map a DNA position to zero or more DNA positions.

        Args:
            feature (str): Feature ID or name
            start (int): Start position
            end (Optional[int], optional): End position. Defaults to the same as `start`.
            strand (Optional[str], optional): Strand ('+' or '-'). Defaults to checking either.
            start_offset (Optional[int], optional): Nucleotide offset from `start`. Defaults to 0.
            end_offset (Optional[int], optional): Nucleotide offset from `end`. Defaults to 0.
            refseq (str, optional): Reference allele. Defaults to "".
            altseq (str, optional): Alternate allele. Defaults to "".
            canonical (bool, optional): Only map to the canonical transcript. Defaults to False.

        Returns:
            Union[List[DnaPosition], List[_DnaSmallVariant]]:
                If `refseq` or `altseq` was given, returns a variant. Otherwise returns a position.
        """
        if refseq or altseq:
            return self._map_variant(
                self.contig_ids,
                self._dna_to_dna_variant,
                feature,
                start,
                refseq,
                altseq,
                start_offset,
                end,
                end_offset,
                strand,
                canonical,
            )
        else:
            return self._map(
                self.contig_ids,
                self._dna_to_dna,
                feature,
                start,
                start_offset,
                end,
                end_offset,
                strand,
                canonical,
            )

    def dna_to_exon(
        self,
        feature: str,
        start: int,
        end: Optional[int] = None,
        strand: Optional[str] = None,
        start_offset: Optional[int] = None,
        end_offset: Optional[int] = None,
        refseq: str = "",
        altseq: str = "",
        canonical: bool = False,
    ) -> Union[List[ExonPosition], List[_ExonSmallVariant]]:
        """Map a DNA position to zero or more exon positions.

        Args:
            feature (str): Feature ID or name
            start (int): Start position
            end (Optional[int], optional): End position. Defaults to the same as `start`.
            strand (Optional[str], optional): Strand ('+' or '-'). Defaults to checking either.
            start_offset (Optional[int], optional): Nucleotide offset from `start`. Defaults to 0.
            end_offset (Optional[int], optional): Nucleotide offset from `end`. Defaults to 0.
            refseq (str, optional): Reference allele. Defaults to "".
            altseq (str, optional): Alternate allele. Defaults to "".
            canonical (bool, optional): Only map to the canonical transcript. Defaults to False.

        Returns:
            Union[List[ExonPosition], List[_ExonSmallVariant]]:
                If `refseq` or `altseq` was given, returns a variant. Otherwise returns a position.
        """
        if refseq or altseq:
            return self._map_variant(
                self.contig_ids,
                self._dna_to_exon_variant,
                feature,
                start,
                refseq,
                altseq,
                start_offset,
                end,
                end_offset,
                strand,
                canonical,
            )
        else:
            return self._map(
                self.contig_ids,
                self._dna_to_exon,
                feature,
                start,
                start_offset,
                end,
                end_offset,
                strand,
                canonical,
            )

    def dna_to_protein(
        self,
        feature: str,
        start: int,
        end: Optional[int] = None,
        strand: Optional[str] = None,
        start_offset: Optional[int] = None,
        end_offset: Optional[int] = None,
        refseq: str = "",
        altseq: str = "",
        canonical: bool = False,
    ) -> Union[List[ProteinPosition], List[_ProteinSmallVariant]]:
        """Map a DNA position to zero or more protein positions.

        Args:
            feature (str): Feature ID or name
            start (int): Start position
            end (Optional[int], optional): End position. Defaults to the same as `start`.
            strand (Optional[str], optional): Strand ('+' or '-'). Defaults to checking either.
            start_offset (Optional[int], optional): Nucleotide offset from `start`. Defaults to 0.
            end_offset (Optional[int], optional): Nucleotide offset from `end`. Defaults to 0.
            refseq (str, optional): Reference allele. Defaults to "".
            altseq (str, optional): Alternate allele. Defaults to "".
            canonical (bool, optional): Only map to the canonical transcript. Defaults to False.

        Returns:
            Union[List[ProteinPosition], List[_ProteinSmallVariant]]:
                If `refseq` or `altseq` was given, returns a variant. Otherwise returns a position.
        """
        if refseq or altseq:
            return self._map_variant(
                self.contig_ids,
                self._dna_to_protein_variant,
                feature,
                start,
                refseq,
                altseq,
                start_offset,
                end,
                end_offset,
                strand,
                canonical,
            )
        else:
            return self._map(
                self.contig_ids,
                self._dna_to_protein,
                feature,
                start,
                start_offset,
                end,
                end_offset,
                strand,
                canonical,
            )

    def dna_to_rna(
        self,
        feature: str,
        start: int,
        end: Optional[int] = None,
        strand: Optional[str] = None,
        start_offset: Optional[int] = None,
        end_offset: Optional[int] = None,
        refseq: str = "",
        altseq: str = "",
        canonical: bool = False,
    ) -> Union[List[RnaPosition], List[_RnaSmallVariant]]:
        """Map a DNA position to zero or more RNA positions.

        Args:
            feature (str): Feature ID or name
            start (int): Start position
            end (Optional[int], optional): End position. Defaults to the same as `start`.
            strand (Optional[str], optional): Strand ('+' or '-'). Defaults to checking either.
            start_offset (Optional[int], optional): Nucleotide offset from `start`. Defaults to 0.
            end_offset (Optional[int], optional): Nucleotide offset from `end`. Defaults to 0.
            refseq (str, optional): Reference allele. Defaults to "".
            altseq (str, optional): Alternate allele. Defaults to "".
            canonical (bool, optional): Only map to the canonical transcript. Defaults to False.

        Returns:
            Union[List[RnaPosition], List[_RnaSmallVariant]]:
                If `refseq` or `altseq` was given, returns a variant. Otherwise returns a position.
        """
        if refseq or altseq:
            return self._map_variant(
                self.contig_ids,
                self._dna_to_rna_variant,
                feature,
                start,
                refseq,
                altseq,
                start_offset,
                end,
                end_offset,
                strand,
                canonical,
            )
        else:
            return self._map(
                self.contig_ids,
                self._dna_to_rna,
                feature,
                start,
                start_offset,
                end,
                end_offset,
                strand,
                canonical,
            )

    def exon_to_cdna(
        self,
        feature: str,
        start: Optional[int] = None,
        end: Optional[int] = None,
        strand: Optional[str] = None,
        start_offset: Optional[int] = None,
        end_offset: Optional[int] = None,
        canonical: bool = False,
    ) -> List[CdnaPosition]:
        """Map an exon position to zero or more cDNA positions.

        Args:
            feature (str): Feature ID or name
            start (int): Start position
            end (Optional[int], optional): End position. Defaults to the same as `start`.
            strand (Optional[str], optional): Strand ('+' or '-'). Defaults to checking either.
            start_offset (Optional[int], optional): Nucleotide offset from `start`. Defaults to 0.
            end_offset (Optional[int], optional): Nucleotide offset from `end`. Defaults to 0.
            canonical (bool, optional): Only map to the canonical transcript. Defaults to False.

        Returns:
            List[CdnaPosition]: Zero or more cDNA positions
        """
        return self._map_exon(
            self._exon_to_cdna, feature, start, start_offset, end, end_offset, strand, canonical
        )

    def exon_to_dna(
        self,
        feature: str,
        start: Optional[int] = None,
        end: Optional[int] = None,
        strand: Optional[str] = None,
        start_offset: Optional[int] = None,
        end_offset: Optional[int] = None,
        canonical: bool = False,
    ) -> List[DnaPosition]:
        """Map an exon position to zero or more DNA positions.

        Args:
            feature (str): Feature ID or name
            start (int): Start position
            end (Optional[int], optional): End position. Defaults to the same as `start`.
            strand (Optional[str], optional): Strand ('+' or '-'). Defaults to checking either.
            start_offset (Optional[int], optional): Nucleotide offset from `start`. Defaults to 0.
            end_offset (Optional[int], optional): Nucleotide offset from `end`. Defaults to 0.
            canonical (bool, optional): Only map to the canonical transcript. Defaults to False.

        Returns:
            List[DnaPosition]: Zero or more DNA positions
        """
        return self._map_exon(
            self._exon_to_dna, feature, start, start_offset, end, end_offset, strand, canonical
        )

    def exon_to_exon(
        self,
        feature: str,
        start: Optional[int] = None,
        end: Optional[int] = None,
        strand: Optional[str] = None,
        start_offset: Optional[int] = None,
        end_offset: Optional[int] = None,
        canonical: bool = False,
    ) -> List[ExonPosition]:
        """Map an exon position to zero or more exon positions.

        Args:
            feature (str): Feature ID or name
            start (int): Start position
            end (Optional[int], optional): End position. Defaults to the same as `start`.
            strand (Optional[str], optional): Strand ('+' or '-'). Defaults to checking either.
            start_offset (Optional[int], optional): Nucleotide offset from `start`. Defaults to 0.
            end_offset (Optional[int], optional): Nucleotide offset from `end`. Defaults to 0.
            canonical (bool, optional): Only map to the canonical transcript. Defaults to False.

        Returns:
            List[ExonPosition]: Zero or more exon positions
        """
        return self._map_exon(
            self._exon_to_exon, feature, start, start_offset, end, end_offset, strand, canonical
        )

    def exon_to_protein(
        self,
        feature: str,
        start: Optional[int] = None,
        end: Optional[int] = None,
        strand: Optional[str] = None,
        start_offset: Optional[int] = None,
        end_offset: Optional[int] = None,
        canonical: bool = False,
    ) -> List[ProteinPosition]:
        """Map an exon position to zero or more protein positions.

        Args:
            feature (str): Feature ID or name
            start (int): Start position
            end (Optional[int], optional): End position. Defaults to the same as `start`.
            strand (Optional[str], optional): Strand ('+' or '-'). Defaults to checking either.
            start_offset (Optional[int], optional): Nucleotide offset from `start`. Defaults to 0.
            end_offset (Optional[int], optional): Nucleotide offset from `end`. Defaults to 0.
            canonical (bool, optional): Only map to the canonical transcript. Defaults to False.

        Returns:
            List[ProteinPosition]: Zero or more protein positions
        """
        return self._map_exon(
            self._exon_to_protein, feature, start, start_offset, end, end_offset, strand, canonical
        )

    def exon_to_rna(
        self,
        feature: str,
        start: Optional[int] = None,
        end: Optional[int] = None,
        strand: Optional[str] = None,
        start_offset: Optional[int] = None,
        end_offset: Optional[int] = None,
        canonical: bool = False,
    ) -> List[RnaPosition]:
        """Map an exon position to zero or more RNA positions.

        Args:
            feature (str): Feature ID or name
            start (int): Start position
            end (Optional[int], optional): End position. Defaults to the same as `start`.
            strand (Optional[str], optional): Strand ('+' or '-'). Defaults to checking either.
            start_offset (Optional[int], optional): Nucleotide offset from `start`. Defaults to 0.
            end_offset (Optional[int], optional): Nucleotide offset from `end`. Defaults to 0.
            canonical (bool, optional): Only map to the canonical transcript. Defaults to False.

        Returns:
            List[RnaPosition]: Zero or more RNA positions
        """
        return self._map_exon(
            self._exon_to_rna, feature, start, start_offset, end, end_offset, strand, canonical
        )

    def protein_to_cdna(
        self,
        feature: str,
        start: int,
        end: Optional[int] = None,
        strand: Optional[str] = None,
        start_offset: Optional[int] = None,
        end_offset: Optional[int] = None,
        refseq: str = "",
        altseq: str = "",
        canonical: bool = False,
    ) -> Union[List[CdnaPosition], List[_CdnaSmallVariant]]:
        """Map a protein position to zero or more cDNA positions.

        Args:
            feature (str): Feature ID or name
            start (int): Start position
            end (Optional[int], optional): End position. Defaults to the same as `start`.
            strand (Optional[str], optional): Strand ('+' or '-'). Defaults to checking either.
            start_offset (Optional[int], optional): Nucleotide offset from `start`. Defaults to 0.
            end_offset (Optional[int], optional): Nucleotide offset from `end`. Defaults to 0.
            refseq (str, optional): Reference allele. Defaults to "".
            altseq (str, optional): Alternate allele. Defaults to "".
            canonical (bool, optional): Only map to the canonical transcript. Defaults to False.

        Returns:
            Union[List[CdnaPosition], List[_CdnaSmallVariant]]:
                If `refseq` or `altseq` was given, returns a variant. Otherwise returns a position.
        """
        if refseq or altseq:
            return self._map_variant(
                self.transcript_ids,
                self._protein_to_cdna_variant,
                feature,
                start,
                refseq,
                altseq,
                start_offset,
                end,
                end_offset,
                strand,
                canonical,
            )
        else:
            return self._map(
                self.transcript_ids,
                self._protein_to_cdna,
                feature,
                start,
                start_offset,
                end,
                end_offset,
                strand,
                canonical,
            )

    def protein_to_dna(
        self,
        feature: str,
        start: int,
        end: Optional[int] = None,
        strand: Optional[str] = None,
        start_offset: Optional[int] = None,
        end_offset: Optional[int] = None,
        refseq: str = "",
        altseq: str = "",
        canonical: bool = False,
    ) -> Union[List[DnaPosition], List[_DnaSmallVariant]]:
        """Map a protein position to zero or more DNA positions.

        Args:
            feature (str): Feature ID or name
            start (int): Start position
            end (Optional[int], optional): End position. Defaults to the same as `start`.
            strand (Optional[str], optional): Strand ('+' or '-'). Defaults to checking either.
            start_offset (Optional[int], optional): Nucleotide offset from `start`. Defaults to 0.
            end_offset (Optional[int], optional): Nucleotide offset from `end`. Defaults to 0.
            refseq (str, optional): Reference allele. Defaults to "".
            altseq (str, optional): Alternate allele. Defaults to "".
            canonical (bool, optional): Only map to the canonical transcript. Defaults to False.

        Returns:
            Union[List[DnaPosition], List[_DnaSmallVariant]]:
                If `refseq` or `altseq` was given, returns a variant. Otherwise returns a position.
        """
        if refseq or altseq:
            return self._map_variant(
                self.transcript_ids,
                self._protein_to_dna_variant,
                feature,
                start,
                refseq,
                altseq,
                start_offset,
                end,
                end_offset,
                strand,
                canonical,
            )
        else:
            return self._map(
                self.transcript_ids,
                self._protein_to_dna,
                feature,
                start,
                start_offset,
                end,
                end_offset,
                strand,
                canonical,
            )

    def protein_to_exon(
        self,
        feature: str,
        start: int,
        end: Optional[int] = None,
        strand: Optional[str] = None,
        start_offset: Optional[int] = None,
        end_offset: Optional[int] = None,
        refseq: str = "",
        altseq: str = "",
        canonical: bool = False,
    ) -> Union[List[ExonPosition], List[_ExonSmallVariant]]:
        """Map a protein position to zero or more exon positions.

        Args:
            feature (str): Feature ID or name
            start (int): Start position
            end (Optional[int], optional): End position. Defaults to the same as `start`.
            strand (Optional[str], optional): Strand ('+' or '-'). Defaults to checking either.
            start_offset (Optional[int], optional): Nucleotide offset from `start`. Defaults to 0.
            end_offset (Optional[int], optional): Nucleotide offset from `end`. Defaults to 0.
            refseq (str, optional): Reference allele. Defaults to "".
            altseq (str, optional): Alternate allele. Defaults to "".
            canonical (bool, optional): Only map to the canonical transcript. Defaults to False.

        Returns:
            Union[List[ExonPosition], List[_ExonSmallVariant]]:
                If `refseq` or `altseq` was given, returns a variant. Otherwise returns a position.
        """
        if refseq or altseq:
            return self._map_variant(
                self.transcript_ids,
                self._protein_to_exon_variant,
                feature,
                start,
                refseq,
                altseq,
                start_offset,
                end,
                end_offset,
                strand,
                canonical,
            )
        else:
            return self._map(
                self.transcript_ids,
                self._protein_to_exon,
                feature,
                start,
                start_offset,
                end,
                end_offset,
                strand,
                canonical,
            )

    def protein_to_protein(
        self,
        feature: str,
        start: int,
        end: Optional[int] = None,
        strand: Optional[str] = None,
        start_offset: Optional[int] = None,
        end_offset: Optional[int] = None,
        refseq: str = "",
        altseq: str = "",
        canonical: bool = False,
    ) -> Union[List[ProteinPosition], List[_ProteinSmallVariant]]:
        """Map a protein position to zero or more protein positions.

        Args:
            feature (str): Feature ID or name
            start (int): Start position
            end (Optional[int], optional): End position. Defaults to the same as `start`.
            strand (Optional[str], optional): Strand ('+' or '-'). Defaults to checking either.
            start_offset (Optional[int], optional): Nucleotide offset from `start`. Defaults to 0.
            end_offset (Optional[int], optional): Nucleotide offset from `end`. Defaults to 0.
            refseq (str, optional): Reference allele. Defaults to "".
            altseq (str, optional): Alternate allele. Defaults to "".
            canonical (bool, optional): Only map to the canonical transcript. Defaults to False.

        Returns:
            Union[List[ProteinPosition], List[_ProteinSmallVariant]]:
                If `refseq` or `altseq` was given, returns a variant. Otherwise returns a position.
        """
        if refseq or altseq:
            return self._map_variant(
                self.transcript_ids,
                self._protein_to_protein_variant,
                feature,
                start,
                refseq,
                altseq,
                start_offset,
                end,
                end_offset,
                strand,
                canonical,
            )
        else:
            return self._map(
                self.transcript_ids,
                self._protein_to_protein,
                feature,
                start,
                start_offset,
                end,
                end_offset,
                strand,
                canonical,
            )

    def protein_to_rna(
        self,
        feature: str,
        start: int,
        end: Optional[int] = None,
        strand: Optional[str] = None,
        start_offset: Optional[int] = None,
        end_offset: Optional[int] = None,
        refseq: str = "",
        altseq: str = "",
        canonical: bool = False,
    ) -> Union[List[RnaPosition], List[_RnaSmallVariant]]:
        """Map a protein position to zero or more RNA positions.

        Args:
            feature (str): Feature ID or name
            start (int): Start position
            end (Optional[int], optional): End position. Defaults to the same as `start`.
            strand (Optional[str], optional): Strand ('+' or '-'). Defaults to checking either.
            start_offset (Optional[int], optional): Nucleotide offset from `start`. Defaults to 0.
            end_offset (Optional[int], optional): Nucleotide offset from `end`. Defaults to 0.
            refseq (str, optional): Reference allele. Defaults to "".
            altseq (str, optional): Alternate allele. Defaults to "".
            canonical (bool, optional): Only map to the canonical transcript. Defaults to False.

        Returns:
            Union[List[RnaPosition], List[_RnaSmallVariant]]:
                If `refseq` or `altseq` was given, returns a variant. Otherwise returns a position.
        """
        if refseq or altseq:
            return self._map_variant(
                self.transcript_ids,
                self._protein_to_rna_variant,
                feature,
                start,
                refseq,
                altseq,
                start_offset,
                end,
                end_offset,
                strand,
                canonical,
            )
        else:
            return self._map(
                self.transcript_ids,
                self._protein_to_rna,
                feature,
                start,
                start_offset,
                end,
                end_offset,
                strand,
                canonical,
            )

    def rna_to_cdna(
        self,
        feature: str,
        start: int,
        end: Optional[int] = None,
        strand: Optional[str] = None,
        start_offset: Optional[int] = None,
        end_offset: Optional[int] = None,
        refseq: str = "",
        altseq: str = "",
        canonical: bool = False,
    ) -> Union[List[CdnaPosition], List[_CdnaSmallVariant]]:
        """Map an RNA position to zero or more cDNA positions.

        Args:
            feature (str): Feature ID or name
            start (int): Start position
            end (Optional[int], optional): End position. Defaults to the same as `start`.
            strand (Optional[str], optional): Strand ('+' or '-'). Defaults to checking either.
            start_offset (Optional[int], optional): Nucleotide offset from `start`. Defaults to 0.
            end_offset (Optional[int], optional): Nucleotide offset from `end`. Defaults to 0.
            refseq (str, optional): Reference allele. Defaults to "".
            altseq (str, optional): Alternate allele. Defaults to "".
            canonical (bool, optional): Only map to the canonical transcript. Defaults to False.

        Returns:
            Union[List[CdnaPosition], List[_CdnaSmallVariant]]:
                If `refseq` or `altseq` was given, returns a variant. Otherwise returns a position.
        """
        if refseq or altseq:
            return self._map_variant(
                self.transcript_ids,
                self._rna_to_cdna_variant,
                feature,
                start,
                refseq,
                altseq,
                start_offset,
                end,
                end_offset,
                strand,
                canonical,
            )
        else:
            return self._map(
                self.transcript_ids,
                self._rna_to_cdna,
                feature,
                start,
                start_offset,
                end,
                end_offset,
                strand,
                canonical,
            )

    def rna_to_dna(
        self,
        feature: str,
        start: int,
        end: Optional[int] = None,
        strand: Optional[str] = None,
        start_offset: Optional[int] = None,
        end_offset: Optional[int] = None,
        refseq: str = "",
        altseq: str = "",
        canonical: bool = False,
    ) -> Union[List[DnaPosition], List[_DnaSmallVariant]]:
        """Map an RNA position to zero or more DNA positions.

        Args:
            feature (str): Feature ID or name
            start (int): Start position
            end (Optional[int], optional): End position. Defaults to the same as `start`.
            strand (Optional[str], optional): Strand ('+' or '-'). Defaults to checking either.
            start_offset (Optional[int], optional): Nucleotide offset from `start`. Defaults to 0.
            end_offset (Optional[int], optional): Nucleotide offset from `end`. Defaults to 0.
            refseq (str, optional): Reference allele. Defaults to "".
            altseq (str, optional): Alternate allele. Defaults to "".
            canonical (bool, optional): Only map to the canonical transcript. Defaults to False.

        Returns:
            Union[List[DnaPosition], List[_DnaSmallVariant]]:
                If `refseq` or `altseq` was given, returns a variant. Otherwise returns a position.
        """
        if refseq or altseq:
            return self._map_variant(
                self.transcript_ids,
                self._rna_to_dna_variant,
                feature,
                start,
                refseq,
                altseq,
                start_offset,
                end,
                end_offset,
                strand,
                canonical,
            )
        else:
            return self._map(
                self.transcript_ids,
                self._rna_to_dna,
                feature,
                start,
                start_offset,
                end,
                end_offset,
                strand,
                canonical,
            )

    def rna_to_exon(
        self,
        feature: str,
        start: int,
        end: Optional[int] = None,
        strand: Optional[str] = None,
        start_offset: Optional[int] = None,
        end_offset: Optional[int] = None,
        refseq: str = "",
        altseq: str = "",
        canonical: bool = False,
    ) -> Union[List[ExonPosition], List[_ExonSmallVariant]]:
        """Map an RNA position to zero or more exon positions.

        Args:
            feature (str): Feature ID or name
            start (int): Start position
            end (Optional[int], optional): End position. Defaults to the same as `start`.
            strand (Optional[str], optional): Strand ('+' or '-'). Defaults to checking either.
            start_offset (Optional[int], optional): Nucleotide offset from `start`. Defaults to 0.
            end_offset (Optional[int], optional): Nucleotide offset from `end`. Defaults to 0.
            refseq (str, optional): Reference allele. Defaults to "".
            altseq (str, optional): Alternate allele. Defaults to "".
            canonical (bool, optional): Only map to the canonical transcript. Defaults to False.

        Returns:
            Union[List[ExonPosition], List[_ExonSmallVariant]]:
                If `refseq` or `altseq` was given, returns a variant. Otherwise returns a position.
        """
        if refseq or altseq:
            return self._map_variant(
                self.transcript_ids,
                self._rna_to_exon_variant,
                feature,
                start,
                refseq,
                altseq,
                start_offset,
                end,
                end_offset,
                strand,
                canonical,
            )
        else:
            return self._map(
                self.transcript_ids,
                self._rna_to_exon,
                feature,
                start,
                start_offset,
                end,
                end_offset,
                strand,
                canonical,
            )

    def rna_to_protein(
        self,
        feature: str,
        start: int,
        end: Optional[int] = None,
        strand: Optional[str] = None,
        start_offset: Optional[int] = None,
        end_offset: Optional[int] = None,
        refseq: str = "",
        altseq: str = "",
        canonical: bool = False,
    ) -> Union[List[ProteinPosition], List[_ProteinSmallVariant]]:
        """Map an RNA position to zero or more protein positions.

        Args:
            feature (str): Feature ID or name
            start (int): Start position
            end (Optional[int], optional): End position. Defaults to the same as `start`.
            strand (Optional[str], optional): Strand ('+' or '-'). Defaults to checking either.
            start_offset (Optional[int], optional): Nucleotide offset from `start`. Defaults to 0.
            end_offset (Optional[int], optional): Nucleotide offset from `end`. Defaults to 0.
            refseq (str, optional): Reference allele. Defaults to "".
            altseq (str, optional): Alternate allele. Defaults to "".
            canonical (bool, optional): Only map to the canonical transcript. Defaults to False.

        Returns:
            Union[List[ProteinPosition], List[_ProteinSmallVariant]]:
                If `refseq` or `altseq` was given, returns a variant. Otherwise returns a position.
        """
        if refseq or altseq:
            return self._map_variant(
                self.transcript_ids,
                self._rna_to_protein_variant,
                feature,
                start,
                refseq,
                altseq,
                start_offset,
                end,
                end_offset,
                strand,
                canonical,
            )
        else:
            return self._map(
                self.transcript_ids,
                self._rna_to_protein,
                feature,
                start,
                start_offset,
                end,
                end_offset,
                strand,
                canonical,
            )

    def rna_to_rna(
        self,
        feature: str,
        start: int,
        end: Optional[int] = None,
        strand: Optional[str] = None,
        start_offset: Optional[int] = None,
        end_offset: Optional[int] = None,
        refseq: str = "",
        altseq: str = "",
        canonical: bool = False,
    ) -> Union[List[RnaPosition], List[_RnaSmallVariant]]:
        """Map an RNA position to zero or more RNA positions.

        Args:
            feature (str): Feature ID or name
            start (int): Start position
            end (Optional[int], optional): End position. Defaults to the same as `start`.
            strand (Optional[str], optional): Strand ('+' or '-'). Defaults to checking either.
            start_offset (Optional[int], optional): Nucleotide offset from `start`. Defaults to 0.
            end_offset (Optional[int], optional): Nucleotide offset from `end`. Defaults to 0.
            refseq (str, optional): Reference allele. Defaults to "".
            altseq (str, optional): Alternate allele. Defaults to "".
            canonical (bool, optional): Only map to the canonical transcript. Defaults to False.

        Returns:
            Union[List[RnaPosition], List[_RnaSmallVariant]]:
                If `refseq` or `altseq` was given, returns a variant. Otherwise returns a position.
        """
        if refseq or altseq:
            return self._map_variant(
                self.transcript_ids,
                self._rna_to_rna_variant,
                feature,
                start,
                refseq,
                altseq,
                start_offset,
                end,
                end_offset,
                strand,
                canonical,
            )
        else:
            return self._map(
                self.transcript_ids,
                self._rna_to_rna,
                feature,
                start,
                start_offset,
                end,
                end_offset,
                strand,
                canonical,
            )

    def _cdna_to_cdna(
        self,
        transcript_ids: List[str],
        start: int,
        start_offset: int,
        end: int,
        end_offset: int,
        strand: List[str],
        canonical: bool,
        include_stop: bool = True,
    ) -> List[CdnaPosition]:
        feature = ["CDS", "stop_codon"] if include_stop else ["CDS"]

        if canonical:
            transcript_ids = [i for i in transcript_ids if self.is_canonical_transcript(i)]

        def convert(position: int, offset: int):
            result = []

            # If the offset position can be normalized to a non-offset position, do so. Otherwise
            # just return an offset position.
            if offset:
                for dna in self._cdna_to_dna(
                    transcript_ids,
                    position,
                    offset,
                    position,
                    offset,
                    strand,
                    canonical,
                    include_stop=include_stop,
                ):
                    for cdna in dna.to_cdna(canonical=canonical):
                        if cdna.transcript_id in transcript_ids:
                            result.append(cdna)

                if result:
                    return result

            mask = (
                (self.df[TRANSCRIPT_ID].isin(transcript_ids))
                & (self.df["cdna_start"] <= position)
                & (self.df["cdna_end"] >= position)
                & (self.df["strand"].isin(strand))
                & (self.df["feature"].isin(feature))
            )
            for _, cds in self.df[mask].iterrows():
                result.append(
                    CdnaPosition(
                        _data=self,
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
        transcript_ids: List[str],
        start: int,
        start_offset: int,
        end: int,
        end_offset: int,
        strand: List[str],
        refseq: str,
        altseq: str,
        canonical: bool,
        include_stop: bool = True,
    ) -> List[_CdnaSmallVariant]:
        result: List[_CdnaSmallVariant] = []

        for cdna in self._cdna_to_cdna(
            transcript_ids,
            start,
            start_offset,
            end,
            end_offset,
            strand,
            canonical,
            include_stop=include_stop,
        ):
            result.extend(_CdnaSmallVariant.from_cdna(cdna, refseq, altseq))

        return result

    def _cdna_to_dna(
        self,
        transcript_ids: List[str],
        start: int,
        start_offset: int,
        end: int,
        end_offset: int,
        strand: List[str],
        canonical: bool,
        include_stop: bool = True,
    ) -> List[DnaPosition]:
        feature = ["CDS", "stop_codon"] if include_stop else ["CDS"]

        if canonical:
            transcript_ids = [i for i in transcript_ids if self.is_canonical_transcript(i)]

        def convert(position: int, offset: int):
            result = []

            mask = (
                (self.df[TRANSCRIPT_ID].isin(transcript_ids))
                & (self.df["cdna_start"] <= position)
                & (self.df["cdna_end"] >= position)
                & (self.df["strand"].isin(strand))
                & (self.df["feature"].isin(feature))
            )
            for _, cds in self.df[mask].iterrows():
                if cds.strand == "-":
                    new_start = new_end = cds.end - (position - cds.cdna_start) - offset
                else:
                    new_start = new_end = cds.start + (position - cds.cdna_start) + offset

                # TODO: Check that new new_start is actually on the contig
                result.append(
                    DnaPosition(
                        _data=self,
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
        transcript_ids: List[str],
        start: int,
        start_offset: int,
        end: int,
        end_offset: int,
        strand: List[str],
        refseq: str,
        altseq: str,
        canonical: bool,
        include_stop: bool = True,
    ) -> List[_DnaSmallVariant]:
        result: List[_DnaSmallVariant] = []

        for dna in self._cdna_to_dna(
            transcript_ids,
            start,
            start_offset,
            end,
            end_offset,
            strand,
            canonical,
            include_stop=include_stop,
        ):
            result.extend(_DnaSmallVariant.from_dna(dna, refseq, altseq))

        return result

    def _cdna_to_exon(
        self,
        transcript_ids: List[str],
        start: int,
        start_offset: int,
        end: int,
        end_offset: int,
        strand: List[str],
        canonical: bool,
        include_stop: bool = True,
    ) -> List[ExonPosition]:
        feature = ["CDS", "stop_codon"] if include_stop else ["CDS"]

        if canonical:
            transcript_ids = [i for i in transcript_ids if self.is_canonical_transcript(i)]

        def convert(position: int, offset: int):
            result = []

            # For an offset position, we need to calculate then equivalent DNA position then map
            # that to one or more exons. This is slower, so if there's no offset we can just map
            # directly to an exon.
            if offset:
                for dna in self._cdna_to_dna(
                    transcript_ids,
                    position,
                    offset,
                    position,
                    offset,
                    strand,
                    canonical,
                    include_stop=include_stop,
                ):
                    for exon in dna.to_exon(canonical=canonical):
                        if exon.transcript_id in transcript_ids:
                            result.append(exon)

                return result

            mask_cds = (
                (self.df[TRANSCRIPT_ID].isin(transcript_ids))
                & (self.df["cdna_start"] <= position)
                & (self.df["cdna_end"] >= position)
                & (self.df["strand"].isin(strand))
                & (self.df["feature"].isin(feature))
            )
            for _, cds in self.df[mask_cds].iterrows():
                mask_exon = (
                    (self.df[TRANSCRIPT_ID].isin(transcript_ids))
                    & (self.df["exon_number"] == cds.exon_number)
                    & (self.df["feature"] == "exon")
                )
                for _, exon_row in self.df[mask_exon].iterrows():
                    result.append(
                        ExonPosition(
                            _data=self,
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
        transcript_ids: List[str],
        start: int,
        start_offset: int,
        end: int,
        end_offset: int,
        strand: List[str],
        refseq: str,
        altseq: str,
        canonical: bool,
        include_stop: bool = True,
    ) -> List[_ExonSmallVariant]:
        result: List[_ExonSmallVariant] = []

        for exon in self._cdna_to_exon(
            transcript_ids,
            start,
            start_offset,
            end,
            end_offset,
            strand,
            canonical,
            include_stop=include_stop,
        ):
            result.extend(_ExonSmallVariant.from_exon(exon, refseq, altseq))

        return result

    def _cdna_to_protein(
        self,
        transcript_ids: List[str],
        start: int,
        start_offset: int,
        end: int,
        end_offset: int,
        strand: List[str],
        canonical: bool,
    ) -> List[ProteinPosition]:
        result = []
        for cdna in self._cdna_to_cdna(
            transcript_ids,
            start,
            start_offset,
            end,
            end_offset,
            strand,
            canonical,
            include_stop=False,
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
        transcript_ids: List[str],
        start: int,
        start_offset: int,
        end: int,
        end_offset: int,
        strand: List[str],
        refseq: str,
        altseq: str,
        canonical: bool,
    ) -> List[_ProteinSmallVariant]:
        # TODO: A lot of this function is duplicated from _cdna_to_protein()
        result: List[_ProteinSmallVariant] = []
        for cdna in self._cdna_to_cdna(
            transcript_ids,
            start,
            start_offset,
            end,
            end_offset,
            strand,
            canonical,
            include_stop=False,
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
            result.extend(_ProteinSmallVariant.from_cdna(cdna, protein, refseq, altseq))

        return result

    def _cdna_to_rna(
        self,
        transcript_ids: List[str],
        start: int,
        start_offset: int,
        end: int,
        end_offset: int,
        strand: List[str],
        canonical: bool,
        include_stop: bool = True,
    ) -> List[RnaPosition]:
        feature = ["CDS", "stop_codon"] if include_stop else ["CDS"]

        if canonical:
            transcript_ids = [i for i in transcript_ids if self.is_canonical_transcript(i)]

        def convert(position: int, offset: int):
            result = []

            # If the offset position can be normalized to a non-offset position, do so. Otherwise
            # just return an offset position.
            if offset:
                for dna in self._cdna_to_dna(
                    transcript_ids,
                    position,
                    offset,
                    position,
                    offset,
                    strand,
                    canonical,
                    include_stop=include_stop,
                ):
                    for rna in dna.to_rna(canonical=canonical):
                        if rna.transcript_id in transcript_ids:
                            result.append(rna)

                if result:
                    return result

            mask = (
                (self.df[TRANSCRIPT_ID].isin(transcript_ids))
                & (self.df["cdna_start"] <= position)
                & (self.df["cdna_end"] >= position)
                & (self.df["strand"].isin(strand))
                & (self.df["feature"].isin(feature))
            )
            for _, cds in self.df[mask].iterrows():
                new_start = new_end = cds.transcript_start + (position - cds.cdna_start)
                result.append(
                    RnaPosition(
                        _data=self,
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
        transcript_ids: List[str],
        start: int,
        start_offset: int,
        end: int,
        end_offset: int,
        strand: List[str],
        refseq: str,
        altseq: str,
        canonical: bool,
        include_stop: bool = True,
    ) -> List[_RnaSmallVariant]:
        result: List[_RnaSmallVariant] = []

        for rna in self._cdna_to_rna(
            transcript_ids,
            start,
            start_offset,
            end,
            end_offset,
            strand,
            canonical,
            include_stop=include_stop,
        ):
            result.extend(_RnaSmallVariant.from_rna(rna, refseq, altseq))

        return result

    def _dna_to_cdna(
        self,
        contig_ids: List[str],
        start: int,
        start_offset: int,
        end: int,
        end_offset: int,
        strand: List[str],
        canonical: bool,
        include_stop: bool = True,
    ) -> List[CdnaPosition]:
        feature = ["CDS", "stop_codon"] if include_stop else ["CDS"]

        def convert(position: int, offset: int):
            result = []

            for strand_ in strand:
                if strand_ == "-":
                    position_ = position - offset
                else:
                    position_ = position + offset

                mask = (
                    (self.df[CONTIG_ID].isin(contig_ids))
                    & (self.df["start"] <= position_)
                    & (self.df["end"] >= position_)
                    & (self.df["strand"] == strand_)
                    & (self.df["feature"].isin(feature))
                )
                for _, cds in self.df[mask].iterrows():
                    if canonical and not self.is_canonical_transcript(cds.transcript_id):
                        continue

                    if cds.strand == "-":
                        new_start = new_end = cds.end - position_ + cds.cdna_start
                    else:
                        new_start = new_end = position_ - cds.start + cds.cdna_start

                    result.append(
                        CdnaPosition(
                            _data=self,
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
        contig_ids: List[str],
        start: int,
        start_offset: int,
        end: int,
        end_offset: int,
        strand: List[str],
        refseq: str,
        altseq: str,
        canonical: bool,
        include_stop: bool = True,
    ) -> List[_CdnaSmallVariant]:
        result: List[_CdnaSmallVariant] = []

        for cdna in self._dna_to_cdna(
            contig_ids,
            start,
            start_offset,
            end,
            end_offset,
            strand,
            canonical,
            include_stop=include_stop,
        ):
            result.extend(_CdnaSmallVariant.from_cdna(cdna, refseq, altseq))

        return result

    def _dna_to_dna(
        self,
        contig_ids: List[str],
        start: int,
        start_offset: int,
        end: int,
        end_offset: int,
        strand: List[str],
        canonical: bool,
    ) -> List[DnaPosition]:
        result = []

        # Sort the start and end positions after adjusting by offsets
        start, end = sorted([start, end])

        for contig_id_, strand_ in product(contig_ids, strand):
            if strand_ == "-":
                new_start = start - start_offset
                new_end = end - end_offset
            else:
                new_start = start + start_offset
                new_end = end + end_offset

            # TODO: Check that new new_start is actually on the contig
            result.append(
                DnaPosition(
                    _data=self,
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
        contig_ids: List[str],
        start: int,
        start_offset: int,
        end: int,
        end_offset: int,
        strand: List[str],
        refseq: str,
        altseq: str,
        canonical: bool,
    ) -> List[_DnaSmallVariant]:
        result: List[_DnaSmallVariant] = []

        for dna in self._dna_to_dna(
            contig_ids, start, start_offset, end, end_offset, strand, canonical
        ):
            result.extend(_DnaSmallVariant.from_dna(dna, refseq, altseq))

        return result

    def _dna_to_exon(
        self,
        contig_ids: List[str],
        start: int,
        start_offset: int,
        end: int,
        end_offset: int,
        strand: List[str],
        canonical: bool,
    ) -> List[ExonPosition]:
        def convert(position: int, offset: int):
            result = []

            for strand_ in strand:
                if strand_ == "-":
                    position_ = position - offset
                else:
                    position_ = position + offset

                mask = (
                    (self.df[CONTIG_ID].isin(contig_ids))
                    & (self.df["start"] <= position_)
                    & (self.df["end"] >= position_)
                    & (self.df["strand"] == strand_)
                    & (self.df["feature"].isin(["exon"]))
                )
                for _, exon in self.df[mask].iterrows():
                    if canonical and not self.is_canonical_transcript(exon.transcript_id):
                        continue

                    result.append(
                        ExonPosition(
                            _data=self,
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
        contig_ids: List[str],
        start: int,
        start_offset: int,
        end: int,
        end_offset: int,
        strand: List[str],
        refseq: str,
        altseq: str,
        canonical: bool,
    ) -> List[_ExonSmallVariant]:
        result: List[_ExonSmallVariant] = []

        for exon in self._dna_to_exon(
            contig_ids, start, start_offset, end, end_offset, strand, canonical
        ):
            result.extend(_ExonSmallVariant.from_exon(exon, refseq, altseq))

        return result

    def _dna_to_protein(
        self,
        contig_ids: List[str],
        start: int,
        start_offset: int,
        end: int,
        end_offset: int,
        strand: List[str],
        canonical: bool,
    ) -> List[ProteinPosition]:
        protein = []
        for cdna in self._dna_to_cdna(
            contig_ids, start, start_offset, end, end_offset, strand, canonical, include_stop=False
        ):
            # Offset cDNA position are assumed to not map to a protein
            if cdna.start_offset or cdna.end_offset:
                continue

            pstart = calc_cdna_to_protein(cdna.start)
            pend = calc_cdna_to_protein(cdna.end)
            protein.append(
                ProteinPosition.copy_from(
                    cdna, start=pstart, start_offset=0, end=pend, end_offset=0, _data=self
                )
            )

        return protein

    def _dna_to_protein_variant(
        self,
        contig_ids: List[str],
        start: int,
        start_offset: int,
        end: int,
        end_offset: int,
        strand: List[str],
        refseq: str,
        altseq: str,
        canonical: bool,
    ) -> List[_ProteinSmallVariant]:
        # TODO: A lot of this function is duplicated from _dna_to_protein()
        result: List[_ProteinSmallVariant] = []
        for cdna in self._dna_to_cdna(
            contig_ids, start, start_offset, end, end_offset, strand, canonical, include_stop=False
        ):
            # Offset cDNA position are assumed to not map to a protein
            if cdna.start_offset or cdna.end_offset:
                continue

            pstart = calc_cdna_to_protein(cdna.start)
            pend = calc_cdna_to_protein(cdna.end)
            protein = ProteinPosition.copy_from(
                cdna, start=pstart, start_offset=0, end=pend, end_offset=0, _data=self
            )
            result.extend(_ProteinSmallVariant.from_cdna(cdna, protein, refseq, altseq))

        return result

    def _dna_to_rna(
        self,
        contig_ids: List[str],
        start: int,
        start_offset: int,
        end: int,
        end_offset: int,
        strand: List[str],
        canonical: bool,
    ) -> List[RnaPosition]:
        def convert(position: int, offset: int):
            result = []

            for strand_ in strand:
                if strand_ == "-":
                    position_ = position - offset
                else:
                    position_ = position + offset

                mask = (
                    (self.df[CONTIG_ID].isin(contig_ids))
                    & (self.df["start"] <= position_)
                    & (self.df["end"] >= position_)
                    & (self.df["strand"] == strand_)
                    & (self.df["feature"] == "exon")
                )
                for _, exon in self.df[mask].iterrows():
                    if canonical and not self.is_canonical_transcript(exon.transcript_id):
                        continue

                    if exon.strand == "-":
                        new_start = new_end = exon.end - position_ + exon.transcript_start
                    else:
                        new_start = new_end = position_ - exon.start + exon.transcript_start

                    result.append(
                        RnaPosition(
                            _data=self,
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
        contig_ids: List[str],
        start: int,
        start_offset: int,
        end: int,
        end_offset: int,
        strand: List[str],
        refseq: str,
        altseq: str,
        canonical: bool,
    ) -> List[_RnaSmallVariant]:
        result: List[_RnaSmallVariant] = []

        for rna in self._dna_to_rna(
            contig_ids, start, start_offset, end, end_offset, strand, canonical
        ):
            result.extend(_RnaSmallVariant.from_rna(rna, refseq, altseq))

        return result

    def _exon_to_cdna(
        self,
        transcript_ids: List[str],
        start: int,
        start_offset: int,
        end: int,
        end_offset: int,
        strand: List[str],
        canonical: bool,
        include_stop: bool = True,
    ) -> List[CdnaPosition]:
        feature = ["CDS", "stop_codon"] if include_stop else ["CDS"]

        if canonical:
            transcript_ids = [i for i in transcript_ids if self.is_canonical_transcript(i)]

        def convert(position: int, offset: int):
            result = []

            # TODO: Is there a reasonable case where an exon position would have an offset?
            assert not offset, offset

            mask = (
                (self.df[TRANSCRIPT_ID].isin(transcript_ids))
                & (self.df["exon_number"] == float(position))
                & (self.df["strand"].isin(strand))
                & (self.df["feature"].isin(feature))
            )
            for _, cds in self.df[mask].iterrows():
                result.append(
                    CdnaPosition(
                        _data=self,
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
        transcript_ids: List[str],
        start: int,
        start_offset: int,
        end: int,
        end_offset: int,
        strand: List[str],
        canonical: bool,
    ) -> List[DnaPosition]:
        if canonical:
            transcript_ids = [i for i in transcript_ids if self.is_canonical_transcript(i)]

        def convert(position: int, offset: int):
            result = []

            # TODO: Is there a reasonable case where an exon position would have an offset?
            assert not offset, offset

            mask = (
                (self.df[TRANSCRIPT_ID].isin(transcript_ids))
                & (self.df["exon_number"] == float(position))
                & (self.df["strand"].isin(strand))
                & (self.df["feature"] == "exon")
            )
            for _, exon in self.df[mask].iterrows():
                result.append(
                    DnaPosition(
                        _data=self,
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
        transcript_ids: List[str],
        start: int,
        start_offset: int,
        end: int,
        end_offset: int,
        strand: List[str],
        canonical: bool,
    ) -> List[ExonPosition]:
        if canonical:
            transcript_ids = [i for i in transcript_ids if self.is_canonical_transcript(i)]

        def convert(position: int, offset):
            result = []

            # TODO: Is there a reasonable case where an exon position would have an offset?
            assert not offset, offset

            mask = (
                (self.df[TRANSCRIPT_ID].isin(transcript_ids))
                & (self.df["exon_number"] == float(position))
                & (self.df["strand"].isin(strand))
                & (self.df["feature"] == "exon")
            )
            for _, exon in self.df[mask].iterrows():
                result.append(
                    ExonPosition(
                        _data=self,
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
        transcript_ids: List[str],
        start: int,
        start_offset: int,
        end: int,
        end_offset: int,
        strand: List[str],
        canonical: bool,
    ) -> List[ProteinPosition]:
        result = []
        for cdna in self._exon_to_cdna(
            transcript_ids,
            start,
            start_offset,
            end,
            end_offset,
            strand,
            canonical,
            include_stop=False,
        ):
            result.extend(cdna.to_protein(canonical=canonical))

        return result

    def _exon_to_rna(
        self,
        transcript_ids: List[str],
        start: int,
        start_offset: int,
        end: int,
        end_offset: int,
        strand: List[str],
        canonical: bool,
    ) -> List[RnaPosition]:
        if canonical:
            transcript_ids = [i for i in transcript_ids if self.is_canonical_transcript(i)]

        def convert(position: int, offset: int):
            result = []

            # TODO: Is there a reasonable case where an exon position would have an offset?
            assert not offset, offset

            mask = (
                (self.df[TRANSCRIPT_ID].isin(transcript_ids))
                & (self.df["exon_number"] == float(position))
                & (self.df["strand"].isin(strand))
                & (self.df["feature"] == "exon")
            )
            for _, exon in self.df[mask].iterrows():
                result.append(
                    RnaPosition(
                        _data=self,
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
        transcript_ids: List[str],
        start: int,
        start_offset: int,
        end: int,
        end_offset: int,
        strand: List[str],
        canonical: bool,
    ) -> List[CdnaPosition]:
        def convert(position: int):
            return ((position - 1) * 3) + 1

        # TODO: Is there a reasonable case where an protein position would have an offset?
        assert not start_offset, start_offset
        assert not end_offset, end_offset

        cdna_start = convert(start)
        cdna_end = convert(end) + 2

        return self._cdna_to_cdna(
            transcript_ids, cdna_start, 0, cdna_end, 0, strand, canonical, include_stop=False
        )

    def _protein_to_cdna_variant(
        self,
        transcript_ids: List[str],
        start: int,
        start_offset: int,
        end: int,
        end_offset: int,
        strand: List[str],
        refseq: str,
        altseq: str,
        canonical: bool,
    ) -> List[_CdnaSmallVariant]:
        result: List[_CdnaSmallVariant] = []

        for cdna in self._protein_to_cdna(
            transcript_ids, start, start_offset, end, end_offset, strand, canonical
        ):
            result.extend(_CdnaSmallVariant.from_protein(cdna, refseq, altseq))

        return result

    def _protein_to_dna(
        self,
        transcript_ids: List[str],
        start: int,
        start_offset: int,
        end: int,
        end_offset: int,
        strand: List[str],
        canonical: bool,
    ) -> List[DnaPosition]:
        result = []
        for cdna in self._protein_to_cdna(
            transcript_ids, start, start_offset, end, end_offset, strand, canonical
        ):
            result.extend(cdna.to_dna(canonical=canonical))

        return sorted(set(result))

    def _protein_to_dna_variant(
        self,
        transcript_ids: List[str],
        start: int,
        start_offset: int,
        end: int,
        end_offset: int,
        strand: List[str],
        refseq: str,
        altseq: str,
        canonical: bool,
    ) -> List[_DnaSmallVariant]:
        result = []
        for cdna in self._protein_to_cdna_variant(
            transcript_ids, start, start_offset, end, end_offset, strand, refseq, altseq, canonical
        ):
            result.extend(cdna.to_dna(canonical=canonical))

        return sorted(set(result))

    def _protein_to_exon(
        self,
        transcript_ids: List[str],
        start: int,
        start_offset: int,
        end: int,
        end_offset: int,
        strand: List[str],
        canonical: bool,
    ) -> List[ExonPosition]:
        result = []
        for cdna in self._protein_to_cdna(
            transcript_ids, start, start_offset, end, end_offset, strand, canonical
        ):
            result.extend(cdna.to_exon(canonical=canonical))

        return sorted(set(result))

    def _protein_to_exon_variant(
        self,
        transcript_ids: List[str],
        start: int,
        start_offset: int,
        end: int,
        end_offset: int,
        strand: List[str],
        refseq: str,
        altseq: str,
        canonical: bool,
    ) -> List[_ExonSmallVariant]:
        result = []
        for cdna in self._protein_to_cdna_variant(
            transcript_ids, start, start_offset, end, end_offset, strand, refseq, altseq, canonical
        ):
            result.extend(cdna.to_exon(canonical=canonical))

        return sorted(set(result))

    def _protein_to_protein(
        self,
        transcript_ids: List[str],
        start: int,
        start_offset: int,
        end: int,
        end_offset: int,
        strand: List[str],
        canonical: bool,
    ) -> List[ProteinPosition]:
        result = []
        for cdna in self._protein_to_cdna(
            transcript_ids, start, start_offset, end, end_offset, strand, canonical
        ):
            result.extend(cdna.to_protein(canonical=canonical))

        return sorted(set(result))

    def _protein_to_protein_variant(
        self,
        transcript_ids: List[str],
        start: int,
        start_offset: int,
        end: int,
        end_offset: int,
        strand: List[str],
        refseq: str,
        altseq: str,
        canonical: bool,
    ) -> List[_ProteinSmallVariant]:
        result = []
        for cdna in self._protein_to_cdna_variant(
            transcript_ids, start, start_offset, end, end_offset, strand, refseq, altseq, canonical
        ):
            result.extend(cdna.to_protein(canonical=canonical))

        return sorted(set(result))

    def _protein_to_rna(
        self,
        transcript_ids: List[str],
        start: int,
        start_offset: int,
        end: int,
        end_offset: int,
        strand: List[str],
        canonical: bool,
    ) -> List[RnaPosition]:
        result = []
        for cdna in self._protein_to_cdna(
            transcript_ids, start, start_offset, end, end_offset, strand, canonical
        ):
            result.extend(cdna.to_rna(canonical=canonical))

        return sorted(set(result))

    def _protein_to_rna_variant(
        self,
        transcript_ids: List[str],
        start: int,
        start_offset: int,
        end: int,
        end_offset: int,
        strand: List[str],
        refseq: str,
        altseq: str,
        canonical: bool,
    ) -> List[_RnaSmallVariant]:
        result = []
        for cdna in self._protein_to_cdna_variant(
            transcript_ids, start, start_offset, end, end_offset, strand, refseq, altseq, canonical
        ):
            result.extend(cdna.to_rna(canonical=canonical))

        return sorted(set(result))

    def _rna_to_cdna(
        self,
        transcript_ids: List[str],
        start: int,
        start_offset: int,
        end: int,
        end_offset: int,
        strand: List[str],
        canonical: bool,
        include_stop: bool = True,
    ) -> List[CdnaPosition]:
        feature = ["CDS", "stop_codon"] if include_stop else ["CDS"]

        if canonical:
            transcript_ids = [i for i in transcript_ids if self.is_canonical_transcript(i)]

        def convert(position: int, offset: int):
            result = []

            # If the offset position can be normalized to a non-offset position, do so. Otherwise
            # just return an offset position.
            if offset:
                for dna in self._rna_to_dna(
                    transcript_ids, position, offset, position, offset, strand, canonical
                ):
                    for cdna in dna.to_cdna(canonical=canonical):
                        if cdna.transcript_id in transcript_ids:
                            result.append(cdna)

                if result:
                    return result

            mask = (
                (self.df[TRANSCRIPT_ID].isin(transcript_ids))
                & (self.df["transcript_start"] <= position)
                & (self.df["transcript_end"] >= position)
                & (self.df["strand"].isin(strand))
                & (self.df["feature"].isin(feature))
            )
            for _, cds in self.df[mask].iterrows():
                new_start = new_end = cds.cdna_start + (position - cds.transcript_start)
                result.append(
                    CdnaPosition(
                        _data=self,
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
        transcript_ids: List[str],
        start: int,
        start_offset: int,
        end: int,
        end_offset: int,
        strand: List[str],
        refseq: str,
        altseq: str,
        canonical: bool,
        include_stop: bool = True,
    ) -> List[_CdnaSmallVariant]:
        result: List[_CdnaSmallVariant] = []

        for cdna in self._rna_to_cdna(
            transcript_ids,
            start,
            start_offset,
            end,
            end_offset,
            strand,
            canonical,
            include_stop=include_stop,
        ):
            result.extend(_CdnaSmallVariant.from_cdna(cdna, refseq, altseq))

        return result

    def _rna_to_dna(
        self,
        transcript_ids: List[str],
        start: int,
        start_offset: int,
        end: int,
        end_offset: int,
        strand: List[str],
        canonical: bool,
    ) -> List[DnaPosition]:
        if canonical:
            transcript_ids = [i for i in transcript_ids if self.is_canonical_transcript(i)]

        def convert(position: int, offset: int):
            result = []

            mask = (
                (self.df[TRANSCRIPT_ID].isin(transcript_ids))
                & (self.df["transcript_start"] <= position)
                & (self.df["transcript_end"] >= position)
                & (self.df["strand"].isin(strand))
                & (self.df["feature"].isin(["exon"]))
            )
            exon_df = self.df[mask]
            for _, exon in exon_df.iterrows():
                if exon.strand == "-":
                    new_start = new_end = exon.end - (position - exon.transcript_start) - offset
                else:
                    new_start = new_end = exon.start + (position - exon.transcript_start) + offset

                # TODO: Check that new new_start is actually on the contig
                result.append(
                    DnaPosition(
                        _data=self,
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
        transcript_ids: List[str],
        start: int,
        start_offset: int,
        end: int,
        end_offset: int,
        strand: List[str],
        refseq: str,
        altseq: str,
        canonical: bool,
    ) -> List[_DnaSmallVariant]:
        result: List[_DnaSmallVariant] = []

        for dna in self._rna_to_dna(
            transcript_ids, start, start_offset, end, end_offset, strand, canonical
        ):
            result.extend(_DnaSmallVariant.from_dna(dna, refseq, altseq))

        return result

    def _rna_to_exon(
        self,
        transcript_ids: List[str],
        start: int,
        start_offset: int,
        end: int,
        end_offset: int,
        strand: List[str],
        canonical: bool,
    ) -> List[ExonPosition]:
        if canonical:
            transcript_ids = [i for i in transcript_ids if self.is_canonical_transcript(i)]

        def convert(position: int, offset: int):
            result = []

            # For an offset position, we need to calculate then equivalent DNA position then map
            # that to one or more exons. This is slower, so if there's no offset we can just map
            # directly to an exon.
            if offset:
                for dna in self._rna_to_dna(
                    transcript_ids, position, offset, position, offset, strand, canonical
                ):
                    for exon in dna.to_exon(canonical=canonical):
                        if exon.transcript_id in transcript_ids:
                            result.append(exon)

                return result

            mask = (
                (self.df[TRANSCRIPT_ID].isin(transcript_ids))
                & (self.df["transcript_start"] <= position)
                & (self.df["transcript_end"] >= position)
                & (self.df["strand"].isin(strand))
                & (self.df["feature"].isin(["exon"]))
            )
            for _, exon_row in self.df[mask].iterrows():
                result.append(
                    ExonPosition(
                        _data=self,
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
        transcript_ids: List[str],
        start: int,
        start_offset: int,
        end: int,
        end_offset: int,
        strand: List[str],
        refseq: str,
        altseq: str,
        canonical: bool,
    ) -> List[_ExonSmallVariant]:
        result: List[_ExonSmallVariant] = []

        for exon in self._rna_to_exon(
            transcript_ids, start, start_offset, end, end_offset, strand, canonical
        ):
            result.extend(_ExonSmallVariant.from_exon(exon, refseq, altseq))

        return result

    def _rna_to_protein(
        self,
        transcript_ids: List[str],
        start: int,
        start_offset: int,
        end: int,
        end_offset: int,
        strand: List[str],
        canonical: bool,
    ) -> List[ProteinPosition]:
        result = []
        for cdna in self._rna_to_cdna(
            transcript_ids,
            start,
            start_offset,
            end,
            end_offset,
            strand,
            canonical,
            include_stop=False,
        ):
            result.extend(cdna.to_protein(canonical=canonical))

        return sorted(set(result))

    def _rna_to_protein_variant(
        self,
        transcript_ids: List[str],
        start: int,
        start_offset: int,
        end: int,
        end_offset: int,
        strand: List[str],
        refseq: str,
        altseq: str,
        canonical: bool,
    ) -> List[_ProteinSmallVariant]:
        result = []
        for cdna in self._rna_to_cdna_variant(
            transcript_ids,
            start,
            start_offset,
            end,
            end_offset,
            strand,
            refseq,
            altseq,
            canonical,
            include_stop=False,
        ):
            result.extend(cdna.to_protein(canonical=canonical))

        return sorted(set(result))

    def _rna_to_rna(
        self,
        transcript_ids: List[str],
        start: int,
        start_offset: int,
        end: int,
        end_offset: int,
        strand: List[str],
        canonical: bool,
    ) -> List[RnaPosition]:
        if canonical:
            transcript_ids = [i for i in transcript_ids if self.is_canonical_transcript(i)]

        def convert(position: int, offset: int):
            result = []

            # If the offset position can be normalized to a non-offset position, do so. Otherwise
            # just return an offset position.
            if offset:
                for dna in self._rna_to_dna(
                    transcript_ids, position, offset, position, offset, strand, canonical
                ):
                    for rna in dna.to_rna(canonical=canonical):
                        if rna.transcript_id in transcript_ids:
                            result.append(rna)

                if result:
                    return result

            mask = (
                (self.df[TRANSCRIPT_ID].isin(transcript_ids))
                & (self.df["transcript_start"] <= position)
                & (self.df["transcript_end"] >= position)
                & (self.df["strand"].isin(strand))
                & (self.df["feature"] == "exon")
            )
            for _, exon in self.df[mask].iterrows():
                result.append(
                    RnaPosition(
                        _data=self,
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
        transcript_ids: List[str],
        start: int,
        start_offset: int,
        end: int,
        end_offset: int,
        strand: List[str],
        refseq: str,
        altseq: str,
        canonical: bool,
    ) -> List[_RnaSmallVariant]:
        result: List[_RnaSmallVariant] = []

        for rna in self._rna_to_rna(
            transcript_ids, start, start_offset, end, end_offset, strand, canonical
        ):
            result.extend(_RnaSmallVariant.from_rna(rna, refseq, altseq))

        return result

    def _normalize_map_args(
        self,
        start: Optional[int],
        start_offset: Optional[int],
        end: Optional[int],
        end_offset: Optional[int],
        strand: Optional[str],
    ) -> Tuple[Optional[int], int, Optional[int], int, List[str]]:
        start_ = start if start is not None else end
        end_ = end if end is not None else start
        start_offset_ = start_offset or 0
        end_offset_ = end_offset or 0
        strand_ = [strand] if strand is not None else ["+", "-"]

        return start_, start_offset_, end_, end_offset_, strand_

    def _map(
        self,
        idfunc: Callable,
        mapfunc: Callable,
        feature: str,
        start: int,
        start_offset: Optional[int],
        end: Optional[int],
        end_offset: Optional[int],
        strand: Optional[str],
        canonical: bool,
    ) -> List:
        start_, start_offset_, end_, end_offset_, strand_ = self._normalize_map_args(
            start, start_offset, end, end_offset, strand
        )
        feature_ = idfunc(feature)
        result = mapfunc(feature_, start_, start_offset_, end_, end_offset_, strand_, canonical)

        return sorted(set(result))

    def _map_variant(
        self,
        idfunc: Callable,
        mapfunc: Callable,
        feature: str,
        start: int,
        refseq: str,
        altseq: str,
        start_offset: Optional[int],
        end: Optional[int],
        end_offset: Optional[int],
        strand: Optional[str],
        canonical: bool,
    ) -> List:
        start_, start_offset_, end_, end_offset_, strand_ = self._normalize_map_args(
            start, start_offset, end, end_offset, strand
        )
        feature_ = idfunc(feature)
        result = mapfunc(
            feature_, start_, start_offset_, end_, end_offset_, strand_, refseq, altseq, canonical
        )

        return sorted(set(result))

    def _map_exon(
        self,
        mapfunc: Callable,
        feature: str,
        start: Optional[int],
        start_offset: Optional[int],
        end: Optional[int],
        end_offset: Optional[int],
        strand: Optional[str],
        canonical: bool,
    ) -> List:
        start_, start_offset_, end_, end_offset_, strand_ = self._normalize_map_args(
            start, start_offset, end, end_offset, strand
        )
        feature_ = self.transcript_ids(feature)

        if start_ is not None:
            positions = list(product(feature_, [start_], [end_], strand_))
        else:
            positions = []
            for exon_id in self.exon_ids(feature):
                positions.extend(self.exon_index(exon_id, transcript_ids=feature_))

        result = []
        for transcript_id, start__, end__, strand__ in positions:
            result.extend(
                mapfunc(
                    [transcript_id],
                    start__,
                    start_offset_,
                    end__,
                    end_offset_,
                    [strand__],
                    canonical,
                )
            )

        return sorted(set(result))

    # ---------------------------------------------------------------------------------------------
    # Utility functions
    # ---------------------------------------------------------------------------------------------
    def cds_offset(self, transcript_id: str) -> int:
        """Get the integer offset of the CDS from the start of the spliced RNA.

        Args:
            transcript_id (str): transcript ID

        Returns:
            int: Nucleotide offset of the cDNA start from the start of the transcript
        """
        rna = self.cdna_to_rna(transcript_id, 1)
        assert len(rna) == 1, f"Ambiguous transcript start for '{transcript_id}'"
        offset = rna[0].start - 1
        assert offset >= 0

        return offset

    def exon_index(
        self, exon_id: str, transcript_ids: Union[List[str], str] = []
    ) -> List[Tuple[str, int, int, str]]:
        """Get the exon number(s) for an exon ID on each transcript. Optionally, restrict to a
        specific transcript(s).

        Args:
            exon_id (str): exon ID
            transcript_ids (Union[List[str], str], optional): Restrict the search to the given
                transcript IDs. Defaults to [].

        Returns:
            List[Tuple[str, int, int, str]]: List of (exon ID, exon index, exon index, exon strand)
        """
        result = []
        restrict = [transcript_ids] if isinstance(transcript_ids, str) else transcript_ids

        mask = (self.df[EXON_ID] == exon_id) & (self.df["feature"] == "exon")
        for _, exon in self.df[mask].iterrows():
            if restrict and exon.transcript_id not in restrict:
                continue
            else:
                result.append((exon.transcript_id, exon.exon_number, exon.exon_number, exon.strand))

        return sorted(set(result))

    def translate_cds_variant(
        self, transcript_id: str, cdna_start: int, cdna_end: int, cdna_alt: str
    ) -> List[str]:
        """Return the mutated protein sequence, given a cDNA position and alt allele.

        Args:
            transcript_id (str): transcript ID
            cdna_start (int): cDNA start position
            cdna_end (int): cDNA end position
            cdna_alt (str): cDNA alternate allele

        Returns:
            List[str]: protein alternate allele(s)
        """
        pep_altseq_set = set()

        # If no alt, assume we're talking about a deletion variant
        if not cdna_alt:
            return [""]

        # Get the codon sequence
        codon_start_offset = (cdna_start - 1) % 3
        codon_start = cdna_start - codon_start_offset
        codon_end_offset = 2 - ((cdna_end - 1) % 3)
        codon_end = cdna_end + codon_end_offset
        codon_refseq = self.cds_sequence(transcript_id, codon_start, codon_end)

        # Assert that the codon sequence is divisible by 3
        assert len(codon_refseq) % 3 == 0

        # Mutate the codon sequence
        codon_refseq_left = codon_refseq[:codon_start_offset]
        codon_refseq_right = codon_refseq[-codon_end_offset:] if codon_end_offset else ""
        for i in expand_nt(cdna_alt):
            codon_altseq = codon_refseq_left + i + codon_refseq_right
            # Assert that the altered codon sequence is divisible by 3
            assert len(codon_altseq) % 3 == 0, codon_altseq
            pep_altseq = "".join(AMINO_ACID_TABLE[codon] for codon in split_by_codon(codon_altseq))
            pep_altseq_set.add(pep_altseq)

        return sorted(pep_altseq_set)


def join_positions(
    start: List[MappablePositionOrSmallVariant],
    end: List[MappablePositionOrSmallVariant],
    merge_on: str,
) -> List[MappablePositionOrSmallVariant]:
    """Return the combination of two list of position or variant objects - one of start positions
    and one of end positions - into one list by the given key (e.g. 'transcript_id').

    All positions must be of the same class.

    Args:
        start (List[MappablePositionOrSmallVariant]): Start positions
        end (List[MappablePositionOrSmallVariant]): End positions
        merge_on (str): Attribute to merge on (e.g. 'transcript_id')

    Returns:
        List[MappablePositionOrSmallVariant]: New position or variant objects
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
