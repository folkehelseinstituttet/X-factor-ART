## DATA

All the data and results that can be made public.

### XchromosomeResultsSexStratified.rds

All the results of EWAS on the CpG probes on X-chromosome.

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

All the results of running DMRff (finding differentially methylated regions).

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
on the net: https://support.illumina.com/array/array_kits/infinium-methylationepic-beadchip-kit/downloads.html).

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
