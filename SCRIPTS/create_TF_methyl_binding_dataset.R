# DESCRIPTION ---------------
#   Get the data extracted from Yin, Y., et al. (2017). 'Impact of cytosine
#   methylation on DNA binding specificities of human transcription factors.'
#   Science, 356(6337), eaaj2239. https://doi.org/10.1126/science.aaj2239 .
#   The supplementary PDF file (Databases_S1_for_submission.pdf) includes
#   all the measured TFs with their classification based on whether their
#   binding sequence includes a CpG or not and based on whether they bind
#   better to the methylated CpG or not.
#   PDF is not a good way of keeping that information; I've exported the text
#   only into DATA/yin_science_2017_TF_types.txt and then ran the script
#   'extract_data_from_yin.sh' to create a file with one TF per line
#   (DATA/extracted_lines_clean.txt).
#   I can analyse this data from there.
#
#  From the supplementary info - some details on special marking:
#    *A*      Signal from Another TF family: Raw sequence data contains reads
#             belonging to another TF family. These are not included to the
#             model presented as the seed excludes them. Information provided
#             to aid bioinformatic analyses of full data.
#    *B*      CG frequency Below cut-off: CG containing subsequences were
#             detected, but frequency of the CG was below the 10% cut-off for
#             calling the effect of methylation on binding.
#    *C*      Call based on bisulfite SELEX: kmer analysis does not detect
#             effect, call is based on bisulfite SELEX. This is caused mostly
#             by absolute requirement for CG in site, in which case the small
#             fraction of unmethylated CG containing ligands that is left after
#             MSssI treatment is enriched.
#    *M*      Motif-dependent effect: TF has multiple motifs, with different
#             effect of CG methylation.
#    *N*      Weak enrichment with Non-specific signal: Enrichment for the motif
#             was weak, there is also signal for enrichment of subsequences
#             that enrich weakly in many weak or failed experiments.
#    *S*      Suboptimal k-mer length: kmer is too short or long compared to
#             the corresponding motif. 8mer analysis does not detect the CG or
#             doesn't detect the CG at the correct position.
#    *U*      Unique reads used to detect signal: Selection library has low
#             diversity after selection, unique reads used to detect signal.
#    *W*      Weak: Obtained motif has relatively low number of reads or weak
#             enrichment.
#
#   Additionally, the abbreviations after the name of the TF mean:
#    *FL*    full length of the protein was tested
#    *DBD*   DNA-binding domain
#    *eDBD*  extended DNA-binding domain
#
# AUTHOR: Julia Romanowska
# DATE CREATED: 2022-01-06
# DATE MODIFIED:

# SETUP -----------------
library(tidyverse)
library(here)

# READ DATA -----------------------
extracted_lines <- read_delim(
  here("DATA", "extracted_lines_clean.txt"),
  delim = ":",
  col_names = FALSE
)
extracted_lines
colnames(extracted_lines) <- c("TF_full_name", "class_full")
extracted_lines

# TIDY ------------------
# remove all the whitespace characters and divide into extra columns
extracted_lines_clean <- extracted_lines %>%
  separate(
    col = TF_full_name,
    into = c("TF_name", "part_tested"),
    sep = "_"
  ) %>%
  separate(
    col = class_full,
    into = c("class", "note"),
    sep = "\\*",
    extra = "merge",
    fill = "right"
  ) %>%
  mutate(across(everything(), str_trim)) %>%
  mutate(across(c(class, part_tested), as.factor))
extracted_lines_clean

skimr::skim(extracted_lines_clean)

# check the various classes
#  - TF name
extracted_lines_clean %>%
  count(TF_name) %>%
  count(n)

# why are there more than one entry for some TFs?
surplus_TF_entries <- extracted_lines_clean %>%
  count(TF_name) %>%
  filter(n > 1) %>%
  pull(TF_name)

extracted_lines_clean %>%
  filter(TF_name %in% surplus_TF_entries)
# looks OK

#  - check part_tested
extracted_lines_clean %>%
  count(part_tested)

#  - check class
extracted_lines_clean %>%
  count(class)

# need to fix wrong typing
extracted_lines_clean <- extracted_lines_clean %>%
  mutate(class = as.character(class)) %>%
  mutate(class = case_when(
    class == "inconclusive" ~ "Inconclusive",
    class == "MethylMInus" ~ "MethylMinus",
    TRUE ~ class
  )) %>%
  mutate(class = as.factor(class))
extracted_lines_clean %>%
  count(class)

#  - check note
extracted_lines_clean %>%
  count(note)

# need to remove the asterisk
extracted_lines_clean <- extracted_lines_clean %>%
  mutate(note = str_replace(note, fixed("*"), ""))

extracted_lines_clean %>%
  count(note)

# there are two types of 'note' that are the same: NS and SN
extracted_lines_clean <- extracted_lines_clean %>%
  mutate(note = ifelse(
    note == "SN",
    yes = "NS",
    no = note
  )) %>%
  mutate(note = as.factor(note))
extracted_lines_clean %>%
  count(note)

# check everything once more
skimr::skim(extracted_lines_clean)

# SAVE ---------------------
write_delim(
  extracted_lines_clean,
  file = here("DATA", "extracted_lines_clean_all.dat"),
  delim = "\t"
)
