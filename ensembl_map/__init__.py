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
