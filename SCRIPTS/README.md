## SCRIPTS

Lots of scripts that can be used to reproduce the analyses and figures.

### create_coMET-like_figure.R

Creating two figures:

- one with the visualization of chromosome positions of EWAS results and
genomic annotations of nearby regions (with the use of 
[`karyoploteR` package](https://bernatgel.github.io/karyoploter_tutorial/))
- second with a rotated half of the correlation matrix (with the help of 
[`ggcorrplot` package](https://rpkgs.datanovia.com/ggcorrplot/) and 
[`grid` package](https://www.stat.auckland.ac.nz/~paul/grid/grid.html))

### gg_qqplot.R

Script for calculating data and plotting QQplots. Adapted from:
https://slowkow.com/notes/ggplot2-qqplot/

### grabGenes.R

Script for extracting human genes and transcripts within a certain region. Uses
the following packages:

- org.Hs.eg.db
- AnnotationDbi
- Biobase
- regioneR
- GenomicRanges
- TxDb.Hsapiens.UCSC.hg19.knownGene
- karyoploteR

### grabRegulRegions.R

Script for extracting regulatory region annotation within a certain region.
Uses the following packages:

- biomaRt
- regioneR
