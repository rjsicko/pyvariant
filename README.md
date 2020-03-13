# Ensembl Map

Map between gene, transcript, exon, CDS, and protein coordinates. This package converts coordinates based on Ensembl annotations.

## Usage

This package is intended to be imported and used as part of a larger workflow e.g.:

```python
>>> import ensembl_map
>>> ensembl_map.gene_to_protein('KIT', 55594077)
[('ENSP00000288135', 621, 621, '+'), ('ENSP00000390987', 617, 617, '+')]
```
