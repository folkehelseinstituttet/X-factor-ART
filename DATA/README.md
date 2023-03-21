## DATA

All the data and results that can be made public.

### XchromosomeResultsSexStratified.rds

All the results of running different XWASes (EWAS on the CpG probes on the X-chromosome) on the MoBa data.

List of `tibbles` named as follows:

- sex of the children,
- internal number of the model used to calculate the association,
- any adjustment variables used;
- e.g., `Girls.parents.1d.bw` means that these are results for girls-only
analysis, where the model (`1d`) includes correction for parental methylation
(`parents`), as well as birth weight and gestational age adjustments (`bw`).

Each tibble has the following variables:

|  _variable name_       |   _format_     |    _description_       |
|:-----------------------|:--------------:|:-----------------------|
| `cpg_id`               | character      | name of the CpG probe  |
| `effect_size`          | numeric        | effect size            |
| `std_err`              | numeric        | standard error         |
| `ps`                   | numeric        | p-value                |
| `ps_adj_BH`            | numeric        | FDR-adjusted p-value   |
| `CHR`                  | character      | chromosome (X)         |
| `MAPINFO`              | numeric        | coordinates in GRCh37  |
| `Strand`               | character      | 'F' (forward) or 'R' (reverse) |

### DMRFFResults.rds

All the results of running DMRff (finding differentially methylated regions) on the MoBa data.

The R-package can be found here: https://github.com/perishky/dmrff

List of `tibbles` named as described above.

Each tibble has the following variables:

|  _variable name_       |   _format_     |    _description_       |
|:-----------------------|:--------------:|:-----------------------|
| `chr`                  | character      | chromosome (X)         |
| `start`                | numeric        | start coordinate       |
| `end`                  | numeric        | end coordinate         |
| `n`                    | numeric        | number of CpG sites within the DMR |
| `estimate`             | numeric        |                        |
| `se`                   | numeric        | std. error             |
| `z`                    | numeric        | test statistic         |
| `p.value`              | numeric        | p-value                |
| `p.adjust`             | numeric        | FDR-adjusted p-value   |

### CpG_info_Xchrom_manifest.rds

`tibble` with the information about all the CpG probes that are on X chromosome
on the EPIC array. Taken from the original Illumina Manifest file (available
online at https://support.illumina.com/array/array_kits/infinium-methylationepic-beadchip-kit/downloads.html.

Some important fields are listed below:

|  _variable name_       |   _format_     |    _description_       |
|:-----------------------|:--------------:|:-----------------------|
| `Name`                 | character      | CpG id (equal to `IllmnID`) |
| `CHR`                  | character      | chromosome (X)         |
| `MAPINFO`              | numeric        | coordinates in GRCh37  |

### Genes and regulatory regions

All the files beginning with `cur_genes_` or `cur_regul_regs_` contain
a `GRanges` object with genes or regulatory regions, respectively, that holds
specific annotations in the regions around the indicated CpG. These `GRanges`
objects can be used in subsequent call to `karyoploteR` plotting functions.

These can be read in as follows:

```{r}
genes <- readRDS(
  here("DATA", "cur_genes_girls_model1_cg26175661.rds.txt")
)

regul_regs <- readRDS(
  here("DATA", "cur_regul_regs_girls_model1_cg26175661.rds.txt")
)
```

Analogously, there are `DMR_genes.rds` and `DMR_regul_regs.rds` that hold
`GRanges` objects with genes or regulatory regions within significant DMRs.
`DMR_genes_ensembl.rds` holds the same information as `DMR_genes.rds` but in
a tibble format, taken from ensembl through query via `biomaRt`.

For convenience, these data are also saved in tab-delimited text files beginning
with `genes_in_DMRs_` or `regul_regs_in_DMRs_`.

### Matrices of DNA methylation correlation

All the files beginning with `correlation_plot_` or `correlation_matrix_`
contain a `ggplot` object or a `tibble`, respectively, that holds the
correlation coefficients between all pairs of CpGs around the CpG indicated
in the file name.

These can be read in as follows:

```{r}
corr_plot <- readRDS(
  here("DATA", "correlation_plot_girls_model1_cg26175661.rds.txt")
)

corr_matrix <- readRDS(
  here("DATA", "correlation_matrix_girls_model1_cg26175661.rds.txt")
)
```

### Data cleaning and analysis of data on transcription factor (TF) binding to methylated DNA

As the supplementary material to the original publication by
Yin, Y., et al. in Science (356(6337), 2017. https://doi.org/10.1126/science.aaj2239)
was not presented in an analysis-friendly manner, we extracted the necessary
information from the PDF document with figures and gathered it here.
(See [README in the SCRIPTS folder](../SCRIPTS/README.md) for details).

- `extracted_lines_clean.txt` - raw extracted text
- `extracted_lines_clean.txt` - cleaned raw data
- `extracted_lines_clean_all.dat` - final dataset in the format:

|  _variable name_       |   _format_     |    _description_       |
|:-----------------------|:--------------:|:-----------------------|
| `TF_name`              | character      | TF name (gene name)    |
| `part_tested`          | character      | which part of the protein was tested? (check Yin et al. for details) |
| `class`                | character      | classification of TF according to sensitivity to methylation of DNA binding site(s) |
| `note`                 | character      | comments (if any)       |

- `yin_science_2017_TF_types.txt` gives a short overview of the classification
of TFs made in Yin et al. study

### Checking significant CpGs vs TFs that bind methylated DNA

- `TFs_binding_to_CpGs_JASPAR_ensembl.dat` - data collected from JASPAR 2022 TFBS hg19 via visualization in ensembl; format:

|  _variable name_       |   _format_     |    _description_       |
|:-----------------------|:--------------:|:-----------------------|
| `CpG`                  | character      | CpG name               |
| `TF_name`              | character      | TF name (gene name)    |
| `in_MeDReader`         | boolean        | included in MeDReader DB or not? |

- `TFs_binding_signif_CpGs_cleaned_all.dat` - all the TFs that bind to the significant CpGs and were checked in Yin et al. study; format:

|  _variable name_       |   _format_     |    _description_       |
|:-----------------------|:--------------:|:-----------------------|
| `CpG`                  | character      | CpG name               |
| `TF_name`              | character      | TF name (gene name)    |
| `in_MeDReader`         | boolean        | included in MeDReader DB or not? |
| `part_tested`          | character      | which part of the protein was tested? (check Yin et al. for details) |
| `class`                | character      | classification of TF according to sensitivity to methylation of DNA binding site(s) |
| `note`                 | character      | comments (if any)       |

- `TFs_methyl_sensitive_signif_CpGs.dat` - only methylation-sensitive TFs from the file above (the same format)

## Results from the external cohort

The files starting with `CHART_` contain XWAS results for the replication cohort
(CHART cohort, Australia), using `limma` R package. Format:

|  _variable name_       |   _format_     |    _description_       |
|:-----------------------|:--------------:|:-----------------------|
|                        | character      | CpG name (rowname)     |  
|   logFC                | numeric        | log fold change        |
|       t                | numeric        | t statistic            |
| P.Value                | numeric        | raw p-value            |
| adj.P.Val              | numeric        | FDR-adjusted p-value   |
|       B                | numeric        | |
| Male_ART_mean          | numeric        | mean value of beta (DNAm) for males conceived through ART |
| Female_ART_mean        | numeric        | mean value of beta (DNAm) for females conceived through ART |
| Male_nonART_mean       | numeric        | mean value of beta (DNAm) for males conceived naturally |
| Female_nonART_mean     | numeric        | mean value of beta (DNAm) for females conceived naturally |
| bacon_pvalue           | numeric        | p-value after adjustment with BACON algortithm |

Moreover, there are also `.rds` files, where the genes and regulatory regions
co-localized with the significant results were collected. _(currently, empty)_

## Bootstrapping results

Rank tables after bootstrapping of the XWAS results on MoBa data are available
in `bootstrap_results_boys_final.csv` and `bootstrap_results_girls_final.csv`.
