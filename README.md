# Ensembl Map

- [Ensembl Map](#ensembl-map)
  - [Description](#description)
  - [Quick Start](#quick-start)
  - [How to Install](#how-to-install)
  - [How to Build Annotations](#how-to-build-annotations)
  - [How to Use](#how-to-use)
  - [Examples](#examples)

## Description

This is a Python package for mapping sequence variants between chromosome, gene, exon, protein, and transcript representations.

## Quick Start

TODO

## How to Install

```pytest
pip install ensembl_map
```

## How to Build Annotations

TODO

## How to Use

```python
>>> from ensembl_map import EnsemblRelease
>>> release = EnsemblRelease(species='homo_sapiens', release=100)
```

```python
>>> positions = release.protein_to_cdna("ENSP00000309572", 525)
>>> positions
[CdnaPosition(_data=EnsemblRelease(species=homo_sapiens, release=100), contig_id='5', start=1573, end=1575, strand='-', gene_id='ENSG00000164362', gene_name='TERT', transcript_id='ENST00000310581', transcript_name='TERT-201', protein_id='ENSP00000309572')]
```

```python
>>> positions[0].to_dna()
[DnaPosition(_data=EnsemblRelease(species=homo_sapiens, release=100), contig_id='5', start=1282623, end=1293313, strand='-')]
```

## Examples

TODO
