"""Definitions for the `EnsemblRelease` class."""
from __future__ import annotations

from typing import List, cast

from .core import Core
from .ensembl_cache import EnsemblCache
from .files import tsv_to_dict, txt_to_list
from .positions import CdnaPosition, DnaPosition, ProteinPosition, RnaPosition


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

        # Ensembl does not provide a CDS FASTA so we need to get the sequence from the RNA FASTA
        if position.is_cdna:
            if self.cds_fasta:
                return self._cds_sequence(position.transcript_id, position.start, position.end)
            else:
                position = cast(CdnaPosition, position)
                if rna := self.to_rna(position):
                    position = rna[0]

        if position.is_dna:
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
