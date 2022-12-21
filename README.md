# ensembl_map

## What is it?

**ensembl_map** is a Python package for converting between equivalent chromosome, cDNA, gene, exon, protein, and transcript positions.

## How to get it

The easiest way to get ensembl_map is using [pip](https://pip.pypa.io/en/latest/quickstart.html):

```sh
pip install ensembl_map
```

The source code is hosted on GitHub at: <https://github.com/mattdoug604/ensembl_map>

## How to use it

Before you can use ensembl_map, you will need to download the necessary genomic data. 

To download and install Ensembl data, run:

```shell
ensembl_map install --species <species-name> --release <Ensembl-release-number>
```

For example:

```shell
ensembl_map install --species homo_sapiens --release 100
```

By default, the data is downloaded to a [platform-specific data directory](https://pypi.org/project/appdirs/) that is generally only accessible by the user. If you want the data to be accessible to other users, you may to specify a custom data directory with the `--cache` option:

```shell
ensembl_map install --species homo_sapiens --release 100 --cache /path/to/cache/
```

For more options, run:

```shell
ensembl_map install --help
```

Alternatively, you can run the installation from inside a Python process:

```python
>>> from ensembl_map import EnsemblRelease
>>> ensembl100 = EnsemblRelease(species='homo_sapiens', release=100, cache_dir="/path/to/cache/")
>>> ensembl100.install()
```

Once the data is installed, the `EnsemblRelease` object provides methods for getting information about different features, getting the DNA/RNA/protein sequence at a position, and converting between positions, etc. Here are some examples:

Convert from a protein position to a cDNA position:

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

## License

This package is distributed with the [MIT](LICENSE) license.

