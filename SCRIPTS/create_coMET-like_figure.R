# DESCRIPTION: manipulating the correlation matrix and adding chromosome
#   location of the SNPs and CpGs (using karyoploter); to create a figure
#   that looks like from the coMET package
# AUTHOR: Julia Romanowska
# DATE CREATED: 2021.10.15
# DATE UPDATED: 2021.12.27
# ----------------------

# SETUP ----
library(karyoploteR)
library(here)
library(tidyverse)
library(ggcorrplot)
library(biomaRt)
library(grid)

which_gender <- "girls"
which_model <- 1
which_signif <- 5
all_model_names <- c("1c", "parents.1c", "1d.bw", "parents.1d.bw")
which_model_name <- all_model_names[which_model]
which_model_name

p_threshold <- 0.05
min_cpgs_in_dmr <- 3

source("grabRegulRegions.R")
source("grabGenes.R")

# READ DATA ----
load(
  here("DATA", "XchromosomeResults.RData")
)
all_res_names <- ls(pattern = "^Both.|^Boys.|^Girls.")
rm(list = all_res_names)

# EWAS results
ewas_results <- readRDS(here("DATA", "Xchromosome_results_TEXT.rds"))
names(ewas_results) <- all_res_names

# CpG data
cpg_info <- readRDS(here("DATA", "CpG_info_Xchrom_manifest.rds"))

# DMR data
dmrff_res_all <- readRDS(here("DATA", "DMRFF_results_TEXT.rds"))
names(dmrff_res_all) <- all_res_names

# ADAPT DATA ----
# how many base pairs to show around the chosen EWAS top hit
# margin_zoom <- rep(50000, 2)

# !NB! lower range for cg17479100 due to high density of CpGs there
margin_zoom <- rep(10000, 2)

# !NB! lower range for cg00243584 due to high density of CpGs there
# margin_zoom <- rep(20000, 2)

# get all the CpGs in the DMR
cur_chosen_ewas_res <- ewas_results[[
      paste(
        stringr::str_to_sentence(which_gender),
        which_model_name,
        sep = "."
      )
    ]]$results %>%
      mutate(ps_adj_BH = stats::p.adjust(ps, method = "BH")) %>% 
      dplyr::select(-z_stat, -starts_with("old_"))

# get info about the chosen CPG
best_ewas <- cur_chosen_ewas_res %>%
  arrange(ps, ps_adj_BH) %>%
  slice(which_signif) %>%
  left_join(
    cpg_info %>% 
      dplyr::select(cpg_id = Name, CHR, MAPINFO, Strand)
  )
best_ewas

# read correlation data
corr_plot <- readRDS(
  here("DATA", paste0("correlation_plot_", which_gender, "_model",
                      which_model, "_", best_ewas$cpg_id, ".rds.txt"))
)
# correlation matrix
corr_matrix <- readRDS(
  here("DATA", paste0("correlation_matrix_", which_gender, "_model",
                      which_model, "_", best_ewas$cpg_id, ".rds.txt"))
)
dim(corr_matrix)
colnames(corr_matrix)
corr_matrix

# find the CpGs within the margin
cur_chosen_cpgs <- cpg_info %>%
  filter(MAPINFO <= (best_ewas$MAPINFO + margin_zoom[2]) &
           MAPINFO >= (best_ewas$MAPINFO - margin_zoom[1])) %>%
  dplyr::select(
    CpG = Name, CHR, MAPINFO, Strand, contains("UCSC_")
  ) %>%
  left_join(
    cur_chosen_ewas_res,
    by = c("CpG" = "cpg_id")
  ) %>%
  filter(!is.na(effect_size)) %>%
  arrange(MAPINFO)
cur_chosen_cpgs

# extract the current DMRs (if there are any)
cur_dmrff_signif <- as_tibble(dmrff_res_all[[
    paste(
      stringr::str_to_sentence(which_gender),
      which_model_name,
      sep = "."
    )
  ]]) %>%
  filter(p.adjust < p_threshold & n >= min_cpgs_in_dmr)

cur_region <- c(best_ewas$MAPINFO - margin_zoom[1], best_ewas$MAPINFO + margin_zoom[2])
cur_chosen_dmr <- cur_dmrff_signif %>%
  filter(
    start <= cur_region[2] & (start >= cur_region[1] | end >= cur_region[1])
  )
cur_chosen_dmr

if(nrow(cur_chosen_dmr) > 0){
# this object is used later to plot the DMR as a rectangle
  cur_dmrff_signif_regs <- regioneR::toGRanges(
    data.frame(
      chr = "chrX",
      start = cur_chosen_dmr$start,
      end = cur_chosen_dmr$end,
      n_cpgs = cur_chosen_dmr$n
    )
  )
}

# fetch annotations for this region
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
data_in <- tibble(chr = "X", start = cur_region[1], end = cur_region[2])
colnames(data_in) <- my_filters
data_in

cur_regul_regs <- grabRegulRegions(
  data_in, ensembl_reg, my_attribs, my_filters, asGRanges = TRUE
)
cur_regul_regs

saveRDS(
  cur_regul_regs,
  here("DATA", paste0("cur_regul_regs_", which_gender, "_model", which_model,
                      "_", best_ewas$cpg_id, ".rds"))
)

cur_genes <- grabGenes(data_in, "chrX")
length(cur_genes[[1]]$genes)

saveRDS(
  cur_genes,
  here("DATA", paste0("cur_genes_", which_gender, "_model", which_model, "_",
                      best_ewas$cpg_id, ".rds"))
)

# PLOT ----

# here's the palette
scico::scico_palette_show("berlin")

# karyotype + genomic location of DMRs
png(
  here("FIGURES", paste0("coMET-like_fig_chrom_pos_", which_gender,
                          "_model", which_model, "_", best_ewas$cpg_id, ".png")),
  width = 12,
  height = 8.5,
  units = "in",
  res = 300
)

zoom_region <- regioneR::toGRanges(
    data.frame(
      chr = "chrX",
      start = cur_region[1],
      end = cur_region[2]
    )
  )

pp <- getDefaultPlotParams(plot.type = 2)
pp$topmargin <- 8
pp$data1height <- 100

kp_tmp <- karyoploteR::plotKaryotype(
  chromosome = "chrX",
  zoom = zoom_region,
  plot.type = 2,
  plot.params = pp
)
karyoploteR::kpAddBaseNumbers(
  kp_tmp,
  tick.dist = 10000,
  tick.len = 10,
  minor.ticks = TRUE,
  minor.tick.dist = 5000,
  minor.tick.len = 5,
  cex = 1,
  add.units = TRUE,
  digits = 4
)

# plot EWAS results
# plotting axis on the left side, setting range from data
ewas_max <- max(cur_chosen_cpgs$effect_size, na.rm = TRUE) + 0.05
ewas_min <- min(cur_chosen_cpgs$effect_size, na.rm = TRUE) - 0.1

kpAxis(
  kp_tmp,
  numticks = 6,
  r0 = 0, r1 = 1,
  ymin = ewas_min,
  ymax = ewas_max
)

# setting a nice color scale
colors_acc_pval <- scales::col_numeric(
  palette = scico::scico(100, palette = 'berlin'),
  domain = range(cur_chosen_ewas_res$ps_adj_BH)
)

# plotting background for various layers of data
kpDataBackground(kp_tmp, data.panel = 1)
if(nrow(cur_chosen_dmr) > 0){
  kpDataBackground(
    kp_tmp, data.panel = 2,
    r0 = 0.75, r1 = 0.85,
    color = "#EFFFEF"
  )
}
if(!is.null(cur_regul_regs)){
  kpDataBackground(
    kp_tmp, data.panel = 2,
    r0 = 0.4, r1 = 0.5,
    color = "#D3E1FF"
  )
}

# plotting the EWAS results
karyoploteR::kpPoints(
  kp_tmp,
  chr = "chrX",
  x = cur_chosen_cpgs$MAPINFO,
  y = cur_chosen_cpgs$effect_size,
  ymax = ewas_max,
  ymin = ewas_min,
  data.panel = 1,
  r0 = 0, r1 = 1,
  cex = 1.5,
  col = colors_acc_pval(cur_chosen_cpgs$ps_adj_BH)
)
# adding the y axis title
karyoploteR::kpAddLabels(
  kp_tmp,
  labels = "EWAS effect size",
  data.panel = 1, r0 = 0.8, r1 = 1,
  srt = 90,
  label.margin = 0.08
)

# plot CpGs
karyoploteR::kpPlotMarkers(
  kp_tmp, x = cur_chosen_cpgs$MAPINFO,
  chr = "chrX",
  labels = cur_chosen_cpgs$CpG,
  data.panel = 'all',
  r0 = 1, r1 = 0.05,
  ymin = 0, ymax = 0.65,
  # text.orientation = "horizontal",
  # label.dist = 0.05,
  clipping = FALSE,
  line.color = "grey60",
  pos = 2,
  max.iter = 300
)

# plot genes
if(!is.null(cur_genes)){
  karyoploteR::kpPlotGenes(
    kp_tmp,
    data = cur_genes[[1]],
    data.panel = 2,
    r0 = 0.1,
    r1 = 0.4,
    gene.name.position = "right"
  )
  karyoploteR::kpAddLabels(
    kp_tmp, labels = "genes",
    data.panel = 2, r0 = 0.1, r1 = 0.3,
    label.margin = 0.02
  )
}

# plot the regulatory regions, if any
if(!is.null(cur_regul_regs)){
  karyoploteR::kpPlotRegions(
    kp_tmp,
    data = cur_regul_regs,
    data.panel = 2,
    r0 = 0.4, r1 = 0.5,
    col = "navyblue",
    border = "navyblue"
  )
  karyoploteR::kpPlotMarkers(
    kp_tmp,
    data = cur_regul_regs,
    labels = cur_regul_regs$type,
    data.panel = 2,
    r0 = 0.4, r1 = 0.7,
    # text.orientation = "horizontal",
    # label.dist = 0.01,
    label.color = "navyblue", line.color = "navyblue",
    ignore.chromosome.ends = TRUE,
    clipping = FALSE,
    cex = 0.7
  )
  karyoploteR::kpAddLabels(
    kp_tmp, labels = "and\n regulatory\n regions",
    data.panel = 2,
    r0 = 0.3, r1 = 0.8,
    col = "navyblue",
    label.margin = 0.02
  )
}

if(nrow(cur_chosen_dmr) > 0){
  # plot DMR
  karyoploteR::kpPlotRegions(
    kp_tmp,
    data = cur_dmrff_signif_regs,
    data.panel = 2,
    r0 = 0.75, r1 = 0.85,
    col = "darkgreen", border = "darkgreen"
  )
  karyoploteR::kpAddLabels(
    kp_tmp, labels = "DMR",
    data.panel = 2,
    r0 = 0.8, r1 = 0.85,
    col = "darkgreen",
    label.margin = 0.02
  )
}

dev.off()

# CORRELATION MATRIX, ROTATED ----
png(
  filename = here("FIGURES", paste0("rotated_corr_matrix_", which_gender,
                        "_model", which_model, "_", best_ewas$cpg_id, ".png")),
  width = 17, # for girls - larger matrix
  height = 12, # for girls - larger matrix
  # width = 13, # for girls - larger matrix
  # height = 9, # for girls - larger matrix
  # width = 9,
  # height = 6,
  units = "in",
  res = 300
)

# corr_plot
grid.newpage()
corr_plot +
  theme(
    legend.position = c(1.12, 0.2)
  )
grid.rect(
  x = unit(0.4, "npc"),
  y = unit(0.5, "npc"),
  width = unit(0.91, "npc"),
  # width = unit(0.88, "npc"), # for boys
  height = unit(1, "npc"),
  gp = gpar(fill = "white", col = "white")
)
print(
  corr_plot +
    annotate(
      "text",
      x = 0:(nrow(cur_chosen_cpgs) - 2),
      y = 1:(nrow(cur_chosen_cpgs) - 1),
      label = paste(unique(corr_matrix$parameter1), "-"),
      hjust = 0.8
    ) +
    coord_equal(clip = "off") +
    theme(
      legend.position = "none",
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank()
    ),
  vp = viewport(
    x = unit(0.5, "npc"),
    y = unit(0.7, "npc"),
    width = unit(1, "npc"),
    height = unit(1, "npc"),
    angle = -45
  )
)
grid.text(
  label = "x marks unsignificant correlation\n (adjustment method: BY)",
  x = unit(0.8, "npc"),
  y = unit(0.1, "npc")
)

dev.off()

# SAVE ----
