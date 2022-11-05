## This script does the following:
## 1) Load curated ChIP-seq metadata from Gsheets
## 2) Subsets metadata for completed runs and unique Run title/experiment IDs
## 3) Associates the cryptically named ENCODE pipeline output dirs to each run
## 4) Ensures that each run has expected output and pipeline meta matches curated meta
## 5) Saves out table of completed runs for curating GEO groups in Gsheets, reads curated back in
## 6) Save all tables locally
## -----------------------------------------------------------------------------

library(tidyverse)
library(googlesheets4)
library(rjson)
source("R/setup-01_config.R")

# the root dir where the output of the pipeline lives
chip_dir <- paste0(pipeout_dir, "chip/")

# Output paths of the various tables that are generated:

# table that links qc.jsons to experiment IDs for input to the qc2tsv tool
qc_out <- paste0(pipeout_dir, "qc_reports/batch1_qcfiles_", date, ".txt")

# table associating experiment IDs to its corresponding output directory
run_dir_output <- paste0(meta_dir, "Chipseq/batch1_run_dirs_", date, ".tsv")

# output metadata of all processed ChIP-seq samples
meta_all_output <- paste0(meta_dir, "Chipseq/batch1_chip_meta_all_", date, ".tsv")

# output metadata of all processed ChIP-seq experiments
meta_distinct_output <- paste0(meta_dir, "Chipseq/batch1_chip_meta_completed_runs_", date, ".tsv")

# Output of curated GEO/lab groups for each ChIP-seq experiment
geo_group_output <- paste0(meta_dir, "Chipseq/batch1_chip_geo_groups_", date, ".tsv")


# Load and clean metadata
# MECP2 samples were run on both the TF and the histone pipelines. This involved
# separating their metadata tracking sheets. Ultimately used MECP2 histone
# for analysis, so combine and only keep MECP2 histone runs.
# ------------------------------------------------------------------------------


# Load gsheets of the curated metadata
meta <- read_sheet(
  ss = gsheets_chip,
  sheet = "Master_batch1",
  trim_ws = TRUE,
  col_types = "c"
)

# mecp2 samples that were run in the histone mode of the pipeline
meta_mh <- read_sheet(
  ss = gsheets_chip,
  sheet = "Mecp2_histone",
  trim_ws = TRUE,
  col_types = "c"
)

# Coerce chr "NA" to actual NAs (needed for NA detection)
meta[meta == "NA"] <- NA
meta_mh[meta_mh == "NA"] <- NA


# Combine and create subsets: All experiments; Keeping only completed experiments
# with Mecp2 from the histone pipeline (used in analysis); Distinct meta so 
# that each row is a unique experiment - collapse sample SRX/GSM IDs 


meta_all <- rbind(meta, meta_mh)


meta_comp_all <- rbind(
  filter(meta, str_to_title(Symbol) != "Mecp2", Complete == "1"), 
  filter(meta_mh, Complete == "1")) %>% 
  arrange(Symbol)


collapse_samples <- meta_comp_all %>% 
  group_by(Experiment_ID) %>% 
  mutate(ID = paste(ID, collapse = ", "),
         GSM = paste(GSM, collapse = ", ")) %>% 
  dplyr::select(Experiment_ID, ID, GSM) %>% 
  distinct(Experiment_ID, .keep_all = TRUE)


meta_distinct <- meta_comp_all %>% 
  dplyr::select(-c(ID, GSM)) %>%
  distinct(Experiment_ID, .keep_all = TRUE) %>% 
  left_join(x = ., y = collapse_samples, by = "Experiment_ID") %>% 
  dplyr::select(-Complete)


meta_distinct <- meta_distinct[, order(match(colnames(meta_distinct), colnames(meta_comp_all)))]


stopifnot(identical(
  filter(meta_distinct, Experiment_ID == "GSE114172_Ascl1_Mouse_Dox-OE-12hr")$ID,
  paste(filter(meta_all, Experiment_ID == "GSE114172_Ascl1_Mouse_Dox-OE-12hr")$ID, collapse = ", ")
))



# Each completed ChIP-seq run/experiment is in a cryptically named directory.
# The following gets all dirs, loads the metadata json to get the associated
# experiment ID/run title and samples, ensures that the sample/input IDs
# in that experiment/run match those in the corresponding curated metadata sheet,
# and checks for presence of the peak (IDR/Overlap) subdirs and the QC json.
# A df of the run info and associated directory are returned
# ------------------------------------------------------------------------------


check_meta <- function(run_meta, chip_meta = meta_all) {
  
  # Ensures that the SRX IDs associated in the pathing of the ENCODE
  # run meta json match the SRX IDs in curated ChIP-seq meta. Returns T/F
  
  ctl <- unlist(run_meta$inputs[str_detect(names(run_meta$inputs), "chip.ctl")])
  ctl <- unique(unlist(str_extract_all(ctl, "SRX[:digit:]+")))
  
  sample <- unlist(run_meta$inputs[str_detect(names(run_meta$inputs), "chip.fastqs|chip.tas")])
  sample <- unique(unlist(str_extract_all(sample, "SRX[:digit:]+")))
  
  sub_meta <- chip_meta[chip_meta$Experiment_ID == run_meta$inputs$chip.title, ]
  meta_ctl <- unlist(unique(str_split(sub_meta$Input_ID, ", ")))
  
  # are all input and sample IDs accounted for?
  suppressWarnings(all((ctl %in% meta_ctl) & (sample %in% sub_meta$ID)))
}



get_run_info <- function(dir) {
  
  # Given a directory output from the ENCODE pipeline, load the metadata json
  # file to extract the run title/experiment ID, then check for presence of 
  # peak dir and QC json
  
  df <- data.frame(Experiment_ID = NA,
                   Dir = dir, 
                   IDR = FALSE,
                   QC = FALSE,
                   ID_match = FALSE)
  
  # End if no metadata found
  if (!file.exists(paste0(dir, "/metadata.json"))) {
    return(df)
  }
  
  # load metadata and extract run title
  run_meta <- fromJSON(file = paste0(dir, "/metadata.json"))
  df$Experiment_ID <- run_meta$inputs$chip.title
  
  # check for IDR/overlap presence
  peak_dirs <- list.files(paste0(dir, "/peak"))
  
  if ("idr_reproducibility" %in% peak_dirs) {
    df$IDR <-  TRUE
  } else if ("overlap_reproducibility" %in% peak_dirs) {
    df$IDR <- "Overlap"
  }
  
  # Check for qc.json presence
  df$QC <- file.exists(paste0(dir, "/qc/qc.json"))
  
  # Does the curated meta SRX ID names match meta json pathing?
  df$ID_match <- check_meta(run_meta)
  
  return(df)
}



# Get the dirs from the ChIP-seq output dir
chip_dirs <- list.dirs(chip_dir, recursive = FALSE)

# Get the df of runs, dir, and status checks
run_info <- do.call(rbind, lapply(chip_dirs, get_run_info))


# Tack on symbol/species information (use all meta as want to inspect failures)

run_info <- run_info %>%
  left_join(.,
            meta_all[, c("Symbol", "Species", "Experiment_ID")],
            by = "Experiment_ID") %>%
  distinct(.keep_all = TRUE)


# Inspect/check the output to ensure everything is accounted for
# ------------------------------------------------------------------------------


# check for when no metadata was found - suggests failed run to be removed
stopifnot(nrow(run_info[which(is.na(run_info$Experiment_ID)), ]) == 0)

# Have any runs been duplicated?
stopifnot(nrow(run_info[which(duplicated(run_info$Experiment_ID)), ]) == 0)

# Ensure that IDs match
stopifnot(sum(!run_info$ID_match) == 0)


# Inspect run/experiment ID in meta. Note that some 'additional' runs may
# be present - eg, tests/data outside of Batch 1. Exclude these and experiments
# with Fastq failures as they would not be submitted to the pipeline


missing_exps <- meta_all %>% 
  filter(!str_detect(Complete, "Fastq failure")) %>% 
  filter(!Experiment_ID %in% run_info$Experiment_ID)

stopifnot(nrow(missing_exps) == 0)

run_info <- filter(run_info, (Experiment_ID %in% meta_all$Experiment_ID))


# Are there are any runs where the inputs have not been shared across reps?
# This checks for curation errors on the metadata sheet but also make sure to 
# check corresponding run to see if all inputs were included for pooling

input_per_run <- meta_all %>% 
  group_by(Experiment_ID) %>% 
  summarise(Count_input = n_distinct(Input_ID))

stopifnot(nrow(meta_all[which(input_per_run$Count_input > 1), ]) == 0)


# Is QC missing? NOTE: A qc.json will not be generated if IDR fails, even if
# peak overlap is successful, so turns out this isn't useful. Keep as reminder
# view(run_info[!run_info$QC, ])


# Ensure absent/failed IDR is reflected in metadata
# ------------------------------------------------------------------------------


idr_fail <- run_info %>%
  filter(Experiment_ID %in% meta$Experiment_ID) %>%
  filter(IDR %in% c(FALSE, "Overlap")) %>% 
  pull(Experiment_ID)

meta_idr_fail <- meta %>% 
  filter(!str_detect(Complete, "Fastq failure")) %>% 
  filter(str_detect(str_to_lower(Complete), "idr.*"))

stopifnot(meta_idr_fail$Experiment_ID %in% idr_fail)



# Organize QC file locations for input into the ENCODE qc2tsv tool
# ------------------------------------------------------------------------------


qc_df <- run_info %>% 
  filter(Experiment_ID %in% meta_all$Experiment_ID & QC) %>% 
  mutate(QC_file = paste0(Dir, "/qc/qc.json")) %>% 
  dplyr::select(Experiment_ID, QC_file)



# Save out Gsheets for curating GEO groups (used as blocking variable for bind
# scores). Load curated sheet back to save locally 


geo_dup <- meta_distinct %>% 
  dplyr::select(Experiment_ID, Symbol, GSE) %>% 
  arrange(GSE) %>% 
  dplyr::mutate(Duplicated = duplicated(GSE),
                GEO_Group = NA)

write_sheet(geo_dup, ss = gsheets_chip, sheet = paste0("GEO_groups_", date))


geo_dup_cur <- read_sheet(
  ss = gsheets_chip,
  sheet = paste0("GEO_groups_curated_", date),
  trim_ws = TRUE,
  col_types = "c"
)


# Save out tables
# ------------------------------------------------------------------------------


write.table(
  run_info,
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE,
  sep = "\t",
  file = run_dir_output
)


write.table(
  meta_all,
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE,
  sep = "\t",
  file = meta_all_output
)


write.table(
  meta_distinct,
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE,
  sep = "\t",
  file = meta_distinct_output
)


write.table(
  geo_dup_cur,
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE,
  sep = "\t",
  file = geo_group_output
)


write.table(
  qc_df,
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE,
  sep = "\t",
  file = qc_out
)
