"""Definitions for the `EnsemblRelease` class."""
from __future__ import annotations

from typing import List, Optional

from .core import Core
from .ensembl_cache import EnsemblCache
from .files import tsv_to_dict, txt_to_list


class EnsemblRelease(Core):
    """Handles converting between position types and retrieving information on biological features
    based on an Ensembl release.
    """

    def __init__(
        self,
        species: str,
        release: int,
        cache_dir: str = "",
        canonical_transcript: str = "",
        contig_alias: str = "",
        exon_alias: str = "",
        gene_alias: str = "",
        protein_alias: str = "",
        transcript_alias: str = "",
    ):
        """
        Args:
            species (str): Species name
            release (int): Ensembl release number
            cache_dir (str, optional): Path to the cache directory. Defaults to platform-specific user cache
            canonical_transcript (str, optional): Path to a list of canonical transcripts
            contig_alias (str, optional): Path to a tab-separated file of contig aliases
            exon_alias (str, optional): Path to a tab-separated file of exon aliases
            gene_alias (str, optional): Path to a tab-separated file of gene aliases
            protein_alias (str, optional): Path to a tab-separated file of protein aliases
            transcript_alias (str, optional): Path to a tab-separated file of transcript aliases
        """
        self.ensembl_cache = EnsemblCache(species, release, cache_dir=cache_dir)
        self.cache_dir = self.ensembl_cache.release_cache_dir
        self.reference = self.ensembl_cache.reference
        self.release = self.ensembl_cache.release
        self.species = self.ensembl_cache.species
        self.df = self.ensembl_cache.load_df()
        self.cds_fasta = []
        self.dna_fasta = [self.ensembl_cache.load_dna_fasta()]
        self.protein_fasta = [self.ensembl_cache.load_pep_fasta()]
        self.rna_fasta = [
            self.ensembl_cache.load_cdna_fasta(),
            self.ensembl_cache.load_ncrna_fasta(),
        ]
        self._canonical_transcript = txt_to_list(canonical_transcript)
        self._contig_alias = tsv_to_dict(contig_alias)
        self._exon_alias = tsv_to_dict(exon_alias)
        self._gene_alias = tsv_to_dict(gene_alias)
        self._protein_alias = tsv_to_dict(protein_alias)
        self._transcript_alias = tsv_to_dict(transcript_alias)

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}(species={self.species}, release={self.release})"

    def install(
        self,
        clean: bool = True,
        recache: bool = False,
        redownload: bool = False,
        restrict_genes: List[str] = [],
    ):
        """Download missing data, process, and cache.

        Args:
            clean (bool, optional): Delete temporary files. Defaults to True.
            recache (bool, optional): Overwrite any existing cache. Defaults to False.
            redownload (bool, optional): Redownload files from Ensembl. Defaults to False.
            restrict_genes (List[str], optional): Restrict cache to the specified genes. Defaults to [].
        """
        self.ensembl_cache.install(
            clean=clean, recache=recache, redownload=redownload, restrict_genes=restrict_genes
        )

    def cds_sequence(
        self, transcript_id: str, start: Optional[int] = None, end: Optional[int] = None
    ) -> str:
        """Return the CDS sequence between the given position(s), inclusive.

        Args:
            transcript_id (str): transcript ID
            start (Optional[int], optional): Start position of the sequence to return. Defaults to None.
            end (Optional[int], optional): End position of the sequence to return. Defaults to the same as `start`.

        Returns:
            str: The CDS sequence
        """
        offset = self.cds_offset(transcript_id)
        cds_start = start + offset if start is not None else offset
        cds_end = end + offset if end is not None else offset

        return self._sequence(self.rna_fasta, transcript_id, start=cds_start, end=cds_end)
