# EnsemblMap

- [EnsemblMap](#ensemblmap)
  - [Example](#example)
  - [Installation](#installation)
    - [Download Annotations](#download-annotations)
    - [Custom Annotations](#custom-annotations)

EnsemblMap is a Python package for mapping between chromosome, cDNA, gene, exon, protein, and transcript representations.

## Example

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

## Installation

You can install EnsemblMap using [pip](https://pip.pypa.io/en/latest/quickstart.html):

```sh
pip install ensembl_map
```

This will also install any required Python packages.

Before using EnsemblMap, you will need to either download the annotation files or build them yourself.

### Download Annotations

TODO

### Custom Annotations

TODO
