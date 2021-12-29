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

