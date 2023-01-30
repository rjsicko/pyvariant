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

To download and install data from Ensembl, run:

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

Once the data is installed, the `EnsemblRelease` object provides methods for getting information about different features, getting the DNA/RNA/protein sequence at a position, and converting between positions, etc.

### Map between feature types

Here are a couple examples, with diagrams, to illustrate how positions on the DNA, cDNA, RNA, and protein (and exon) can be equivalent and how `variant-map` can be used map between them.

Examples are using human Ensembl 100. Diagrams are not to scale.

```text
Legend:
=== <- exon
--- <- intron/UTR
```

#### Example 1: KIT (ENSG00000157404.16, + strand)

```text
                 54657957   54658081           54695512
                        |          |           |
DNA:             5' <---============-----------======> 3'
Chromosome 4                  |                  |
                       54658023                  54695514

                        1                      2
                        |                      |
Exons:               5' ============-----------======> 3'
ENSE00000000233 (1)
ENSE00001032350 (2)

                        1        125           126
                        |          |           |
RNA:                 5' ============-----------======> 3'
ENST00000288135.6             |                  |
                             59                  128

                                  67           68
                                   |           |
cDNA:                      5' ======-----------======> 3'
ENST00000288135.6             |                  |
                              1                  70

                                  23
                                   |
Protein:                   5' ======-----------======> 3'
ENSP00000288135.6             |                  |
                              1                  24
```

Map DNA to RNA:

```python
>>> ensembl100.dna_to_rna("4", 54658081)
[RnaMappablePosition(_data=EnsemblRelease(species=homo_sapiens, release=100), contig_id='4', start=125, start_offset=0, end=125, end_offset=0, strand='+', gene_id='ENSG00000157404', gene_name='KIT', transcript_id='ENST00000288135', transcript_name='KIT-201'), RnaMappablePosition(_data=EnsemblRelease(species=homo_sapiens, release=100), contig_id='4', start=164, start_offset=0, end=164, end_offset=0, strand='+', gene_id='ENSG00000157404', gene_name='KIT', transcript_id='ENST00000412167', transcript_name='KIT-202'), RnaMappablePosition(_data=EnsemblRelease(species=homo_sapiens, release=100), contig_id='4', start=143, start_offset=0, end=143, end_offset=0, strand='+', gene_id='ENSG00000157404', gene_name='KIT', transcript_id='ENST00000514582', transcript_name='KIT-204')]
```

Map protein to cDNA:

```python
>>> ensembl100.protein_to_cdna("ENSP00000288135", 23)
[CdnaMappablePosition(_data=EnsemblRelease(species=homo_sapiens, release=100), contig_id='4', start=67, start_offset=0, end=69, end_offset=0, strand='+', gene_id='ENSG00000157404', gene_name='KIT', transcript_id='ENST00000288135', transcript_name='KIT-201', protein_id='ENSP00000288135')]
```

#### Example 2: TERT (ENSG00000164362.21, - strand)

```text
                        1294666     1294771    1295068
                              |           |          |
DNA:                 5' <======-----------============---> 3'
Chromosome 5                                |   |
                                      1294773   1294848

                              2                      1
                              |                      |
Exons:               3' <======-----------============ 5'
ENSE00003896691 (1)
ENSE00001197112 (2)
                            299           298        1
                              |           |          |
RNA:                 3' <======-----------============ 5'
ENST00000310581.10                          |   |
                                          296   80

                            220           219
                              |           |
cDNA:                3' <======-----------======= 5'
ENST00000310581.10                          |   |
                                          217   1

                             74
                              |
Protein              3' <======-----------======= 5'
ENSP00000309572.5                           |   |
                                           73   1
```

Map cDNA to RNA:

```python
>>> ensembl100.cdna_to_rna("ENST00000310581", 219, 220)
[RnaMappablePosition(_data=EnsemblRelease(species=homo_sapiens, release=100), contig_id='5', start=298, start_offset=0, end=299, end_offset=0, strand='-', gene_id='ENSG00000164362', gene_name='TERT', transcript_id='ENST00000310581', transcript_name='TERT-201')]
```

Map RNA to exon:

```python
>>> ensembl100.rna_to_exon("ENST00000310581", 299)
[ExonMappablePosition(_data=EnsemblRelease(species=homo_sapiens, release=100), contig_id='5', start=2, start_offset=0, end=2, end_offset=0, strand='-', gene_id='ENSG00000164362', gene_name='TERT', transcript_id='ENST00000310581', transcript_name='TERT-201', exon_id='ENSE00001197112')]
```

### Retrieve feature information

`variant-map` also has functions for retrieving general information about various features. For example:

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

## Variant naming standards

This package follows [HGVS nomenclature](https://varnomen.hgvs.org/) recommendations for representing variants.


### Offset positions

When describing variants in an intron or UTR, it can be more informative to describe the position relative to a transcript, rather than the genome. These positions are described with a "+" or "-". For example, "TERT:c.-125" means "a position 125 nucleotides 5’ of the ATG translation initiation codon."

Here are some relavent excerpts copied from [HGVS](https://varnomen.hgvs.org/bg-material/numbering/):

> untranslated region (UTR)
> * nucleotides upstream (5’) of the ATG-translation initiation codon (start) are marked with a “-” (minus) and numbered c.-1, c.-2, c.-3, etc. (i.e. going further upstream
> * nucleotides downstream (3’) of the translation termination codon (stop) are marked with a “\*” (asterisk) and numbered c.\*1, c.\*2, c.\*3, etc. (i.e. going further downstream)
>
> introns
> * nucleotides at the 5’ end of an intron are numbered relative to the last nucleotide of the directly upstream exon, followed by a “+” (plus) and their position in to the intron, like c.87+1, c.87+2, c.87+3, …
> * nucleotides at the 3’ end of an intron are numbered relative to the first nucleotide of the directly downstream exon, followed by a “-” (minus) and their position out of the intron, like …, c.88-3, c.88-2, c.88-1.
>   * in the middle of the intron nucleotide numbering changes from “+” (plus) to “-” (minus): e.g. …, c.87+676, c.87+677, c.87+678, c.88-678, c.88-677, c.88-676, …
>   * in introns with an uneven number of nucleotides the central nucleotide is numbered relative to the upstream exon followed by a “+” (plus): e.g. …, c.87+676, c.87+677, c.87+678, c.87+679, c.88-678, c.88-677, c.88-676, …
> * introns in the 5’UTR are numbered as normal introns, starting with “c.-” nucleotide numbers (c.-85+1, c.-85+2, c.-85+3, …, c.-84-3, c.-84-2, c.-84-1)
> * introns in the 3’UTR are numbered as normal introns, starting with “c.\*” nucleotide numbers (c.\*37+1, c.\*37+2, c.\*37+3, …, c.\*38-3, c.\*38-2, c.\*38-1)
> * NOTE: a coding DNA reference sequence does not contain intron or 5’ and 3’ gene flanking sequences and can therefore not be used as a reference to describe variants in these regions see Reference Sequences. Correct descriptions refer to a genomic reference sequence like LRG_199t1:c.357+1G>A, NC_000023.10(NM_004006.2):c.357+1G>A or NG_012232.1(NM_004006.2):c.357+1G>A.


## License

This package is distributed with the [MIT](LICENSE) license.
