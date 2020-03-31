EnsemblMap
==========

EnsemblMap is a python package for mapping between gene, transcript, exon, CDS, and protein coordinates.

Coordinates are based on [Ensembl](http://www.ensembl.org/) annotations which are downloaded and cached on the user's system using [pyensembl](https://github.com/openvax/pyensembl).

# Installation

EnsemblMap requires `python>=3`.

The latest release of EnsemblMap can be installed using [pip](https://pip.pypa.io/en/latest/quickstart/).

```shell
pip install ensembl_map
```

The development version of EnsemblMap can be cloned from git:

```shell
git clone https://github.com/mattdoug604/ensembl_map.git
cd ensembl_map
git checkout develop
pip install .
```

## Cache Location

When the package is run for the first time, Ensembl annotation are downloaded to a platform specific location. You can specify the directory where these files should be downloaded through the package itself:

```python
>>> import ensembl_map
>>> ensembl_map.set_cache_dir("/custom/cache/dir")
```

or, by setting the environmental variable 'PYENSEMBL_CACHE_DIR':

```bash
export PYENSEMBL_CACHE_DIR=/custom/cache/dir
```

# Usage

This package is intended to be imported and used as part of a larger workflow:

```python
>>> import ensembl_map
```

## Setting the Ensembl Release

By default, coordinates are mapped against human Ensembl release 99. A different release can be set by running:

```python
>>> ensembl_map.set_ensembl_release(69, "human")
```

Human Ensembl releases 54 through 99 are supported, along with some releases from other species.

## Coverting Coordinates

Typical usage involves converting from one feature type to another. This is done by calling a function with a name like `ensembl_map.<type>_to<type>()`:
* ```cds_to_exon(feature, start, end=None, raise_error=True)```
* ```cds_to_gene(feature, start, end=None, raise_error=True)```
* ```cds_to_protein(feature, start, end=None, raise_error=True)```
* ```cds_to_transcript(feature, start, end=None, raise_error=True)```
* ```exon_to_cds(feature, raise_error=True)```
* ```exon_to_gene(feature, raise_error=True)```
* ```exon_to_protein(feature, raise_error=True)```
* ```exon_to_transcript(feature, raise_error=True)```
* ```gene_to_cds(feature, start, end=None, raise_error=True)```
* ```gene_to_exon(feature, start, end=None, raise_error=True)```
* ```gene_to_protein(feature, start, end=None, raise_error=True)```
* ```gene_to_transcript(feature, start, end=None, raise_error=True)```
* ```protein_to_cds(feature, start, end=None, raise_error=True)```
* ```protein_to_exon(feature, start, end=None, raise_error=True)```
* ```protein_to_gene(feature, start, end=None, raise_error=True)```
* ```protein_to_transcript(feature, start, end=None, raise_error=True)```
* ```transcript_to_cds(feature, start, end=None, raise_error=True)```
* ```transcript_to_exon(feature, start, end=None, raise_error=True)```
* ```transcript_to_gene(feature, start, end=None, raise_error=True)```
* ```transcript_to_protein(feature, start, end=None, raise_error=True)```

Where `feature` is the name or Ensembl ID of the feature to convert, `start` is the position on that feature, and `end` is an optional second position.

The option, `raise_error`, raises an error if no match for `feature` was found in the Ensembl annotations (otherwise an empty list is returned).

Example:

```python
>>> ensembl_map.gene_to_protein('KIT', 55594077)
[Protein(biotype='protein_coding', contig='4', end=621, protein_id='ENSP00000288135', sequence='A', start=621, strand='+'), Protein(biotype='protein_coding', contig='4', end=617, protein_id='ENSP00000390987', sequence='A', start=617, strand='+')]
```


## Getting Sequences

This package can also be used to get the sequence of a CDS, protein, or transcript at a position:

```python
>>> ensembl_map.protein_sequence("ENSP00000288135", 1)
'M'
>>> ensembl_map.protein_sequence("ENSP00000288135", 112, 191)
'RDPAKLFLVDRSLYGKEDNDTLVRCPLTDPEVTNYSLKGCQGKPLPKDLRFIPDPKAGIMIKSVKRAYHRLCLHCSVDQE'
```

## Extends pyensembl

EnsemblMap acts as an API for pyensembl. The functionality of pyensembl can be accessed directly though:

```python
>>> ensembl_map.Ensembl().data
```

See the [pyensembl docs](https://readthedocs.org/projects/pyensembl/downloads/pdf/latest/) for more info.
