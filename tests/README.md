# Notes on testing

Here are a couple diagrams which can guide you in making test cases. They show how positions on different features relate to each other. Examples are using Ensembl 100. Diagrams are not to scale.

```text
Legend:
=== <- exon
--- <- intron/UTR
```

## Example 1: KIT (ENSG00000157404.16, + strand)

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

## Example 2: TERT (ENSG00000164362.21, - strand)

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

## Mapping offset positions

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
