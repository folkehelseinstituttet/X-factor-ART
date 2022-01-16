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

> NOTE: all these packages are part of BioConductor

### grabRegulRegions.R

Script for extracting regulatory region annotation within a certain region.
Uses the following packages:

- biomaRt
- regioneR

> NOTE: all these packages are part of BioConductor

### extract_data_from_yin.sh

Extracting data from PDF file in Supplementary Materials of the publication:

[Yin, Y. *et al.* Impact of cytosine methylation on DNA binding specificities of human transcription factors. Science, Vol 356 (6337), 2017](http://www.sciencemag.org/lookup/doi/10.1126/science.aaj2239)

This script extracts the information about TFs (transcription factors)
classification based on their binding to methylated and unmethylated DNA
sequences. The data is saved as tab-delimited text file
[`DATA/extracted_lines_clean.txt`](../DATA/extracted_lines_clean.txt),
which is checked and cleaned in
[`create_TF_methyl_binding_dataset.R`](create_TF_methyl_binding_dataset.R)
and analyzed in [`check_TFs_binding_signif_CpGs.R`](check_TFs_binding_signif_CpGs.R).

The output of running the `extract_data_from_yin.sh` script is in
[`extract_data_from_yin.out`](extract_data_from_yin.out).

### create_TF_methyl_binding_dataset.R

Reading in the extracted lines [`DATA/extracted_lines_clean.txt`](../DATA/extracted_lines_clean.txt)
and cleaning the dataset. The final dataset is in
[`DATA/extracted_lines_clean_all.dat`](../DATA/extracted_lines_clean_all.dat).

> NOTE: The script includes detailed information about the data!

Output of running the script is in [`create_TF_methyl_binding_dataset.html`](create_TF_methyl_binding_dataset.html).

### check_TFs_binding_signif_CpGs.R

Checking the classification of the TFs (transcription factors) that were found
to possibly bind to the significant CpGs, based on:

1. data from [JASPAR 2022 db](https://jaspar.genereg.net) visualized in
[ensembl browser, GRCh37](http://grch37.ensembl.org/),
2. manual search in [MeDReader db](http://medreader.org/browse-tf), and
3. data extracted from Yin, Y, et al., Science, 2017, as described above.

The information collected in points *1.* and *2.* was manually entered into
text file [`DATA/TFs_binding_to_CpGs_JASPAR_ensembl.dat`](../DATA/TFs_binding_to_CpGs_JASPAR_ensembl.dat).

The output of running the .R script is in
[`check_TFs_binding_signif_CpGs.html`](check_TFs_binding_signif_CpGs.html).
