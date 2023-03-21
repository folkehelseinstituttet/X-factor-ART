# DESCRIPTION: Results from XWAS on CHART data - plot positions and p-values
# AUTHOR: Julia Romanowska
# DATE CREATED: 2023-03-15
# DATE LAST MODIFIED: 2023-03-21

# SETUP ----
library(karyoploteR)
library(tidyverse)
library(here)
library(biomaRt)
library(conflicted)
conflicts_prefer(
  dplyr::filter,
  dplyr::select,
  dplyr::rename
)

source("grabRegulRegions.R")
source("grabGenes.R")

threshold_signif <- 0.01

# READ DATA ----
xwas_res_girls <- read_csv(
  here("CHART_results",
  "CHART_Female_limma_output_maternal_smoking_and_sample_plate_with_bacon_pvalue.csv"
  )
) %>%
  rename(CpG = `...1`)
xwas_res_girls_all_vars <- read_csv(
  here("CHART_results", 
 "CHART_Xchromosome_FEMALES_nonART_v_ART_LIMMA_output_maternal_smoking_and_sample_plate_variables.csv"
 ), col_select = -`...1`
 ) %>%
  rename(CpG = Row.names) %>%
  select(-starts_with("X"))

xwas_res_girls_all_vars <- xwas_res_girls_all_vars %>%
  select(CpG:pos) %>%
  left_join(
    xwas_res_girls %>%
      select(CpG, t, bacon_pvalue)
  ) %>%
  mutate(adj_bacon_pvalue = p.adjust(bacon_pvalue))
xwas_res_girls_all_vars
xwas_res_girls_all_vars %>% filter(is.na(bacon_pvalue))

xwas_res_boys <- read_csv(
  here("CHART_results",
  "CHART_Male_limma_output_maternal_smoking_and_sample_plate_with_bacon_pvalue.csv"
  )
) %>%
  rename(CpG = `...1`)
xwas_res_boys_all_vars <- read_csv(
  here("CHART_results", 
 "CHART_Xchromosome_Males_nonART_v_ART_LIMMA_output_maternal_smoking_and_sample_plate_variables.csv"
 ), col_select = -`...1`
 ) %>%
  rename(CpG = Row.names) %>%
  select(-starts_with("X"))

xwas_res_boys_all_vars <- xwas_res_boys_all_vars %>%
  select(CpG:pos) %>%
  left_join(
    xwas_res_boys %>%
      select(CpG, t, bacon_pvalue)
  ) %>%
  mutate(adj_bacon_pvalue = p.adjust(bacon_pvalue))
xwas_res_boys_all_vars
xwas_res_boys_all_vars %>% filter(is.na(bacon_pvalue))

all_results_list <- list(boys = xwas_res_boys_all_vars,
                         girls = xwas_res_girls_all_vars)
all_results_signif_list <- list(
  boys = xwas_res_boys_all_vars %>%
    filter(adj_bacon_pvalue < threshold_signif),
  girls = xwas_res_girls_all_vars %>%
    filter(adj_bacon_pvalue < threshold_signif)
)
all_results_signif_list

## Fetch regulatory region annotation from the ensembl ----

regul_regs_chart_file <- here("CHART_results", "CHART_regul_regs.rds")

if(!file.exists(regul_regs_chart_file)){
  ensembl_reg <- useEnsembl("regulation",
    dataset = "hsapiens_regulatory_feature",
    GRCh = 37)
  
  my_attribs <- c(
    "regulatory_stable_id",
    "chromosome_name",
    "chromosome_start",
    "chromosome_end",
    "feature_type_name",
    "feature_type_description",
    "so_accession"
  )
  my_filters <- c(
    "chromosome_name",
    "start", "end"
  )
  
  regul_regs_chart <- map(all_results_signif_list, function(cur_res){
    data_in <- cur_res %>%
      dplyr::select(chr, start = pos) %>%
      mutate(end = start + 1)
    colnames(data_in) <- my_filters
  
    cur_regul_regs <- grabRegulRegions(
      data_in, ensembl_reg, my_attribs, my_filters, asGRanges = TRUE
    )
    return(cur_regul_regs)
  })
  names(regul_regs_chart) <- names(all_results_signif_list)
  
  saveRDS(regul_regs_chart, regul_regs_chart_file)
} else {
  regul_regs_chart <- readRDS(regul_regs_chart_file)
}
sapply(regul_regs_chart, length)

## Fetch genes ----

genes_chart_file <- here("CHART_results", "CHART_genes_regs.rds")

if(!file.exists(genes_chart_file)){
  genes_chart <- map(all_results_signif_list, function(cur_res){
    data_in <- cur_res %>%
      dplyr::select(chr, start = pos) %>%
      mutate(end = start + 1)
    colnames(data_in) <- c("chromosome_name", "start", "end")
  
    cur_genes <- grabGenes(
      data_in, asGRanges = TRUE
    )
    return(cur_genes)
  })
  names(genes_chart) <- names(all_results_signif_list)

  saveRDS(genes_chart, genes_chart_file)
} else {
  genes_chart <- readRDS(genes_chart_file)
}
sapply(genes_chart, function(x){
  nrow(x[[1]]$genes)
})

# grab also ensembl information
genes_ensembl_file <- here("CHART_results", "CHART_genes_ensembl.rds")
if(!file.exists(genes_ensembl_file)){
  ensembl_genes <- useEnsembl("genes",
    dataset = "hsapiens_gene_ensembl",
    GRCh = 37)
  
  all_attribs <- listAttributes(ensembl_genes)
  my_attribs <- c(
    "ensembl_gene_id",
    "description",
    "external_gene_name",
    "start_position", "end_position"
  )
  all_filters <- listFilters(ensembl_genes)
  my_filters <- c(
    "chromosome_name",
    "start", "end"
  )
  
  genes_chart_ensembl <- map(all_results_signif_list, function(cur_res){
    data_in <- cur_res %>%
      dplyr::select(chr, start = pos) %>%
      mutate(end = start + 1)
    colnames(data_in) <- my_filters
  
    cur_regul_regs <- grabRegulRegions(
      data_in, ensembl_genes, my_attribs, my_filters, asGRanges = FALSE
    )
    return(cur_regul_regs)
  })
  names(genes_chart_ensembl) <- names(all_results_signif_list)

  saveRDS(genes_chart_ensembl, file = genes_ensembl_file)
} else {
  genes_chart_ensembl <- readRDS(genes_ensembl_file)
}
sapply(genes_chart_ensembl, length)

# PLOT ----
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene
all.genes <- makeGenesDataFromTxDb(
  txdb, plot.transcripts = FALSE, plot.transcripts.structure = FALSE)$genes
all.genes.tbl <- as_tibble(as.data.frame(all.genes))

walk(names(all_results_list), function(x){
  cur_data_all <- all_results_list[[x]]
  cur_data_signif <- all_results_signif_list[[x]]
  # cur_data_signif_more <- all_results_list[[x]] %>%
    # filter(adj_bacon_pvalue < 0.01)

  cur_genes_regs <- genes_chart[[x]][[1]]
  cur_reg_regs <- regul_regs_chart[[x]]

  cur_y_lim <- range(-log10(cur_data_all$adj_bacon_pvalue))
  
  # tiff(
  #   filename = here("FIGURES", paste0("xwas_res_pval_CHART_", x, ".tiff")),
  #   # width = 1200,
  #   # height = 450,
  #   width = 6.7,
  #   height = 2.5,
  #   units = "in",
  #   pointsize = 18,
  #   res = 300
  # )
  png(
    filename = here("FIGURES", paste0("xwas_res_pval_CHART_", x, ".png")),
    width = 1200,
    height = 450,
    pointsize = 18,
    res = 300
  )
  kp <- plotKaryotype(
    chromosomes = "chrX",
    plot.type = 1
  )
  # plotting gene density
  kpPlotDensity(
    kp, all.genes,
    data.panel = 1, r0 = 0, r1 = 0.4
  )
  kpAddLabels(
    kp, labels = "gene density",
    data.panel = 1, r0 = 0, r1 = 0.3,
    label.margin = -0.02
  )
  # plotting CpG p-values
  #   - first - significant (if any)
  # if(nrow(cur_data_signif_more) != 0){
  #   kpPoints(
  #     kp,
  #     x = cur_data_signif_more$pos,
  #     y = -log10(cur_data_signif_more$adj_bacon_pvalue),
  #     chr = "chrX",
  #     data.panel = 1, cex = 2,
  #     col = "navy",
  #     r0 = 0.5, r1 = 1,
  #     ymin = 0,
  #     ymax = cur_y_lim[2]
  #   )
  # }
  if(nrow(cur_data_signif) != 0){
    kpPoints(
      kp,
      x = cur_data_signif$pos,
      y = -log10(cur_data_signif$adj_bacon_pvalue),
      chr = "chrX",
      data.panel = 1, cex = 1.5,
      col = "cyan",
      r0 = 0.5, r1 = 1,
      ymin = 0,
      ymax = cur_y_lim[2]
    )
  }
  kpPoints(
    kp,
    x = cur_data_all$pos,
    y = -log10(cur_data_all$adj_bacon_pvalue),
    chr = "chrX",
    data.panel = 1,
    r0 = 0.5, r1 = 1,
    ymin = 0,
    ymax = cur_y_lim[2]
  )
  kpAxis(
    kp, data.panel = 1,
    r0 = 0.5, r1 = 1,
    ymin = 0,
    ymax = cur_y_lim[2]
  )
  kpAddLabels(
    kp, labels = "-log10(FDR)",
    srt = 90, data.panel = 1,
    r0 = 1, r1 = 0.9,
    label.margin = 0.07
  )
  if(length(cur_reg_regs) != 0){
    # plotting regulatory regs
    karyoploteR::kpDataBackground(
      kp, data.panel = 2,
      r0 = 0.55, r1 = 0.8,
      color = "#D3E1FF"
    )
    karyoploteR::kpPlotRegions(
      kp, data = cur_reg_regs,
      data.panel = 2, r0 = 0.55, r1 = 0.8,
      col = "navyblue",
      border = "navyblue"
    )
    karyoploteR::kpPlotMarkers(
      kp, data = cur_reg_regs,
      labels = cur_reg_regs$feature_type_name,
      data.panel = 2, r0 = 0, r1 = 1,
      text.orientation = "horizontal",
      label.dist = 0.01,
      label.margin = 5,
      label.color = "navyblue",
      line.color = "navyblue",
      ignore.chromosome.ends = TRUE
    )
    karyoploteR::kpAddLabels(
      kp, labels = "regulatory\n regions",
      data.panel = 2, r0 = 0.6, r1 = 1,
      col = "navyblue"
    )
  }
  if(!is.null(nrow(cur_genes_regs$genes))){
  # plotting genes
  karyoploteR::kpPlotRegions(
    kp, data = cur_genes_regs,
    data.panel = 2, r0 = 0, r1 = 0.2
  )
  karyoploteR::kpPlotMarkers(
    kp, data = cur_genes_regs,
    labels = cur_genes_regs$gene_name,
    data.panel = 2, r0 = 0, r1 = 0.4,
    text.orientation = "horizontal",
    label.dist = 0.01,
    label.margin = 5,
    ignore.chromosome.ends = TRUE
  )
  karyoploteR::kpAddLabels(
    kp, labels = "genes",
    data.panel = 2, r0 = 0.2, r1 = 0.45
  )
  }
  dev.off()
})

## Save for sharing ----
write_csv(
  xwas_res_girls_all_vars %>%
    filter(P.Value < 0.05),
  file = here("DATA", "pval_0.05_from_Female_ART_with_smoke_model.csv")
)

write_csv(
  xwas_res_boys_all_vars %>%
    filter(P.Value < 0.05),
  file = here("DATA", "pval_0.05_from_Male_ART_with_smoke_model.csv")
)
