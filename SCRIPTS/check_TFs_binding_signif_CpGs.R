# DESCRIPTION: ------------
#  I will check which TFs are actually binding _preferrably_ to the methylated
#  sequence, based on data from Yin, Y., et al. (2017). Science, 356(6337),
#  eaaj2239. https://doi.org/10.1126/science.aaj2239 (prepared in 
#  'create_TF_methyl_binding_dataset.R') and data collected from JASPAR 2022
#  TFBS hg19 via visualization in ensembl
#
# AUTHOR: Julia Romanowska
# DATE CREATED: 2022-01-06
# DATE MODIFIED:

# SETUP ------------
library(tidyverse)
library(here)

# READ DATA ----------------
tfs_binding_signif_cpgs <- read_delim(
  here("DATA", "TFs_binding_to_CpGs_JASPAR_ensembl.dat"),
  delim = "\t"
)
tfs_binding_signif_cpgs

binding_tf_methyl <- read_delim(
  here("DATA", "extracted_lines_clean_all.dat"),
  delim = "\t"
)
binding_tf_methyl

# ANALYSE ---------------
tfs_binding_signif_cpgs <- tfs_binding_signif_cpgs %>%
  left_join(binding_tf_methyl)
tfs_binding_signif_cpgs

# check how many were not matched
tfs_binding_signif_cpgs %>%
  filter(is.na(part_tested) & is.na(class))
# some of those are not matched because the name of TF is just an alias for
#   another name, which was matched
tfs_binding_signif_cpgs %>%
  filter(is.na(part_tested) & is.na(class) & in_MeDReader)

# how many of those matched have a note?
tfs_binding_signif_cpgs %>%
  filter(!is.na(part_tested) & !is.na(class)) %>%
  count(note)

tfs_binding_signif_cpgs %>%
  filter(!is.na(part_tested) & !is.na(class)) %>%
  count(note) %>%
  filter(!is.na(note)) %>%
  summarise(sum(n))

# to proceed with a clean dataset, I'll remove those that were not matched
#  as well as those that have unreliable measures (i.e., contain something
#  in the 'note' column)
tfs_binding_signif_cpgs_clean <- tfs_binding_signif_cpgs %>%
  filter((!is.na(part_tested) & !is.na(class)) &
          is.na(note))
DT::datatable(tfs_binding_signif_cpgs_clean)

# ANALYSE - how many duplicated entries? ------------
duplicates <- tfs_binding_signif_cpgs_clean %>%
  count(TF_name) %>%
  filter(n > 1)
duplicates

tfs_binding_signif_cpgs_clean %>%
  filter(TF_name %in% duplicates$TF_name) %>%
  arrange(TF_name)

# how many of those belong to various classes?
tfs_binding_signif_cpgs_clean %>%
  filter(TF_name %in% duplicates$TF_name) %>%
  arrange(TF_name) %>%
  group_by(TF_name) %>%
  count(class) %>%
  filter(n == 1)

# THAT'S GOOD - all the duplicates were basically replicates :)
# let's get rid of the duplicates

#  - first: remove the plain duplicates
nrow(tfs_binding_signif_cpgs_clean)
tfs_binding_signif_cpgs_clean <- tfs_binding_signif_cpgs_clean %>%
  distinct()
nrow(tfs_binding_signif_cpgs_clean)

duplicates <- tfs_binding_signif_cpgs_clean %>%
  count(TF_name) %>%
  filter(n > 1)
duplicates
#  - second: extract those that are not duplicated
tfs_binding_signif_cpgs_clean_part1 <- tfs_binding_signif_cpgs_clean %>%
  filter(!(TF_name %in% duplicates$TF_name))
#  - third: extract those that were duplicated and get only one copy
#           out of there: if there is 'FL' part tested, take it, if not,
#           take any
tfs_binding_signif_cpgs_clean_part2 <- tfs_binding_signif_cpgs_clean %>%
  filter(TF_name %in% duplicates$TF_name)
tfs_binding_signif_cpgs_clean_part2_fl <- tfs_binding_signif_cpgs_clean_part2 %>%
  filter(part_tested == "FL")
tfs_binding_signif_cpgs_clean_part2_fl

tfs_binding_signif_cpgs_clean_part2 %>%
  filter(!(TF_name %in% tfs_binding_signif_cpgs_clean_part2_fl$TF_name)) %>%
  arrange(TF_name)
# these are binding in various places, so that's OK

#  - finally: merge the two parts
tfs_binding_signif_cpgs_clean <- bind_rows(
  tfs_binding_signif_cpgs_clean_part1,
  tfs_binding_signif_cpgs_clean_part2_fl,
  tfs_binding_signif_cpgs_clean_part2 %>%
    filter(!(TF_name %in% tfs_binding_signif_cpgs_clean_part2_fl$TF_name)) %>%
    arrange(TF_name),
  tfs_binding_signif_cpgs %>%
    filter(is.na(part_tested) & is.na(class) & !in_MeDReader)
)
DT::datatable(tfs_binding_signif_cpgs_clean)

# ANALYSE - how many of those TFs can bind to a methylated seq.? ------
# How many various TFs per CpG?
tfs_binding_signif_cpgs_clean %>%
  count(CpG)

# I will do this separately, per CpG
knitr::kable(
  tfs_binding_signif_cpgs_clean %>%
    group_by(CpG) %>%
    count(class),
  caption = "Types of TF binding to a motif containing one of the significant CpGs"
)

# It's the TFs that are in class 'MethylPlus' or 'MethylMinus' that are of
#  interest - those react on methylation status of the binding sequence
knitr::kable(
  tfs_binding_signif_cpgs_clean %>%
    filter(class %in% c("MethylMinus", "MethylPlus")) %>%
    arrange(CpG, TF_name)
)

# SAVE ------------
write_delim(
  tfs_binding_signif_cpgs_clean %>%
    filter(class %in% c("MethylMinus", "MethylPlus")) %>%
    arrange(CpG, TF_name),
  file = here("DATA", "TFs_methyl_sensitive_signif_CpGs.dat"),
  delim = "\t"
)
