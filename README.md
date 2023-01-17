# Variant-Map

## What is it?

**variant-map** is a Python package for mapping biological sequence variants (mutations) to their equivalent chromosome, cDNA, gene, exon, protein, and RNA positions.

## How to get it

The easiest way to get variant-map is using [pip](https://pip.pypa.io/en/latest/quickstart.html):

```sh
pip install variant-map
```

The source code is hosted on GitHub at: <https://github.com/mattdoug604/variant_map>

## How to use it

Before you can use variant-map, you will need to download the necessary genomic data.

To download and install Ensembl data, run:

```shell
variant-map install --species <species-name> --release <Ensembl-release-number>
```

For example:

```shell
variant-map install --species 'homo sapiens' --release 100
```

At the time of writing, installing a human dataset takes roughly 30-45 minutes and 1.5G of storage space. However, the actual time and space required to install a dataset will depend entirely on the size of the dataset, your computer, internet speed, etc.

By default, the data is downloaded to a [platform-specific data directory](https://pypi.org/project/appdirs/) that is generally only accessible by the user (e.g. `/home/<you>/.local/share/variant-map/`). If you want the data to be accessible to other users, or your home directory does have enough storage space, you may to specify a different directory to download to with the `--cache` option:

```shell
variant-map install --species homo_sapiens --release 100 --cache /path/to/cache/
```

For more options, run:

```shell
variant-map install --help
```

Alternatively, you can run the installation from inside a Python process:

```python
>>> from variant-map import EnsemblRelease
>>> ensembl100 = EnsemblRelease(species='homo_sapiens', release=100, cache_dir="/path/to/cache/")
>>> ensembl100.install()
```

Once the data is installed, the `EnsemblRelease` object provides methods for getting information about different features, getting the DNA/RNA/protein sequence at a position, and converting between positions, etc. Here are some examples:

Convert from a protein position to equivalent cDNA positions:

```python
>>> ensembl100.protein_to_cdna("TERT", 525)
[CdnaMappablePosition(_data=EnsemblRelease(species=homo_sapiens, release=100), contig_id='5', start=1573, end=1575, strand='-', gene_id='ENSG00000164362', gene_name='TERT', transcript_id='ENST00000310581', transcript_name='TERT-201', protein_id='ENSP00000309572'), CdnaMappablePosition(_data=EnsemblRelease(species=homo_sapiens, release=100), contig_id='5', start=1573, end=1575, strand='-', gene_id='ENSG00000164362', gene_name='TERT', transcript_id='ENST00000334602', transcript_name='TERT-202', protein_id='ENSP00000334346'), CdnaMappablePosition(_data=EnsemblRelease(species=homo_sapiens, release=100), contig_id='5', start=1573, end=1575, strand='-', gene_id='ENSG00000164362', gene_name='TERT', transcript_id='ENST00000460137', transcript_name='TERT-203', protein_id='ENSP00000425003'), CdnaMappablePosition(_data=EnsemblRelease(species=homo_sapiens, release=100), contig_id='5', start=1573, end=1575, strand='-', gene_id='ENSG00000164362', gene_name='TERT', transcript_id='ENST00000656021', transcript_name='TERT-206', protein_id='ENSP00000499759')]
```

Get a list of transcript IDs for a gene:

```python
>>> ensembl100.transcript_ids("BRCA2")
['ENST00000380152', 'ENST00000470094', 'ENST00000528762', 'ENST00000530893', 'ENST00000533776', 'ENST00000544455', 'ENST00000614259', 'ENST00000665585', 'ENST00000666593', 'ENST00000670614', 'ENST00000671466']
```

Get the DNA sequence at a specific chromosome position:

```python
>>> ensembl100.dna_sequence("5", 1293313, 1293323) 
'CTGGGCTCCT
```

For a complete list of methods, run:

```python
>>> help(EnsemblRelease)
```

## Notes

* This package follows [HGVS nomenclature](https://varnomen.hgvs.org/) recommendations for representing variants.

## License

This package is distributed with the [MIT](LICENSE) license.
