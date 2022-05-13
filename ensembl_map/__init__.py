"""A python package for mapping between chromosome, gene, exon, protein, and transcript coordinates.

This module uses cached Ensembl annotations to perform coordinate mapping, along with 
GSC annotations to map gene symbols or RefSeq IDs to Ensembl IDs. 

Notes:
    * Coordinates start at 1 (not 0)
    * Coordinates are relative to the positive strand of the contig
    * Mapping functions normalize names to their Ensembl ID(s)
    * Since genes commonly have multiple transcripts, coordinates are based on the GSC-
        defined "best transcript" if one is available (see: `select_best_transcript`)
    * Ensembl releases <= 75 use GRCh37; releases > 75 use GRCh38 (see: 
        https://github.com/Ensembl/ensembl-rest/wiki/GRCh37-REST-server)

References:
    * pyensembl: https://readthedocs.org/projects/pyensembl/downloads/pdf/latest/
"""
# from . import ensembl
# from .best_transcripts import get_best_transcript, is_best_transcript  # noqa: F401
# from .gsc import set_ensembl_release
# from .ensembl import (  # noqa: F401
#     get_cache_dir,
#     reference_name,
#     release,
#     release_cache_dir,
#     release_cache_files,
#     species,
# )
# from .ids import (  # noqa: F401
#     get_contig_ids,
#     get_exon_ids,
#     get_gene_ids,
#     get_gene_names,
#     get_protein_ids,
#     get_transcript_ids,
#     get_transcript_names,
#     is_contig,
#     is_exon,
#     is_gene,
#     is_protein,
#     is_transcript,
#     normalize_feature,
# )
# from .map import (  # noqa: F401
#     cds_to_cds,
#     cds_to_contig,
#     cds_to_exon,
#     cds_to_gene,
#     cds_to_protein,
#     cds_to_transcript,
#     contig_to_cds,
#     contig_to_contig,
#     contig_to_exon,
#     contig_to_gene,
#     contig_to_protein,
#     contig_to_transcript,
#     exon_to_cds,
#     exon_to_contig,
#     exon_to_exon,
#     exon_to_gene,
#     exon_to_protein,
#     exon_to_transcript,
#     gene_to_cds,
#     gene_to_contig,
#     gene_to_exon,
#     gene_to_gene,
#     gene_to_protein,
#     gene_to_transcript,
#     get_map_func,
#     protein_to_cds,
#     protein_to_contig,
#     protein_to_exon,
#     protein_to_gene,
#     protein_to_protein,
#     protein_to_transcript,
#     transcript_to_cds,
#     transcript_to_contig,
#     transcript_to_exon,
#     transcript_to_gene,
#     transcript_to_protein,
#     transcript_to_transcript,
# )
# from .normalize import normalize_cds, normalize_protein, normalize_transcript  # noqa: F401
# from .sequences import (  # noqa: F401
#     cds_sequence,
#     contig_sequence,
#     gene_sequence,
#     protein_sequence,
#     transcript_sequence,
# )
# from .utils import is_ensembl_id  # noqa: F401

# # load the default annotations
# set_ensembl_release()
