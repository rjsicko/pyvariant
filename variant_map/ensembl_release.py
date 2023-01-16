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
        self.ensembl_cache = EnsemblCache(species, release, cache_dir=cache_dir)
        self.cache_dir = self.ensembl_cache.release_cache_dir
        self.reference = self.ensembl_cache.reference
        self.release = self.ensembl_cache.release
        self.species = self.ensembl_cache.species
        self.df = self.ensembl_cache.load_df()
        self.cds = []
        self.dna = [self.ensembl_cache.load_dna_fasta()]
        self.peptide = [self.ensembl_cache.load_pep_fasta()]
        self.rna = [self.ensembl_cache.load_cdna_fasta(), self.ensembl_cache.load_ncrna_fasta()]
        self.canonical_transcript = txt_to_list(canonical_transcript)
        self.contig_alias = tsv_to_dict(contig_alias)
        self.exon_alias = tsv_to_dict(exon_alias)
        self.gene_alias = tsv_to_dict(gene_alias)
        self.protein_alias = tsv_to_dict(protein_alias)
        self.transcript_alias = tsv_to_dict(transcript_alias)

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}(species={self.species}, release={self.release})"

    def install(
        self,
        clean: bool = True,
        recache: bool = False,
        redownload: bool = False,
        restrict_genes: List[str] = [],
    ):
        """Download missing data, process, and cache."""
        self.ensembl_cache.install(
            clean=clean, recache=recache, redownload=redownload, restrict_genes=restrict_genes
        )

    def cds_sequence(
        self, transcript_id: str, start: Optional[int] = None, end: Optional[int] = None
    ) -> str:
        """Return the nucleotide sequence at the given CDS coordinates."""
        offset = self.cds_offset(transcript_id)
        cds_start = start + offset if start is not None else offset
        cds_end = end + offset if end is not None else offset

        return self._get_sequence(self.rna, transcript_id, start=cds_start, end=cds_end)
