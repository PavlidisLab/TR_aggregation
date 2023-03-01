## Download and format Chip Atlas's ChIP-seq metadata.
## Only want to keep TF data sets in mouse and human. The downloaded table
## then undergoes additional manual curation to verify details, match inputs, 
## and add additional experiments from other databases.
## -----------------------------------------------------------------------------

library(tidyverse)
library(parallel)
source("R/setup-01_config.R")

# Be aware of date of latest update! https://github.com/inutano/chip-atlas/wiki
# Last update: August 2020 (checked May2022, still same)
ca_date <- "aug2020"

# Chip Atlas download URL
url <- "http://dbarchive.biosciencedbc.jp/kyushu-u/metadata/experimentList.tab"

# Output
out_dir <- file.path(meta_dir, "Chipseq", "Chipatlas")
meta_path <- file.path(out_dir, paste0("chipatlas_all_meta_", ca_date, ".tsv"))

# List of TFs to subset meta
tf_mm <- read.delim(tf_path_mm, stringsAsFactors = FALSE)
tf_hg <- read.delim(tf_path_hg, stringsAsFactors = FALSE)


# Download all chip atlas metadata
# -----------------------------------------------------------------------------


if (!file.exists(meta_path)) {
  download.file(url, meta_path)
}

assertthat::is.readable(meta_path)


# The table has a lot of empty elements - need to get field count. 

nfields <- max(count.fields(meta_path, sep = "\t"), na.rm = TRUE)

meta_all <- read.delim(
  meta_path,
  col.names = 1:nfields,
  quote = "",
  header = FALSE,
  stringsAsFactors = FALSE
)

dim(meta_all) # aug2020: 439593    131


# Format. Note that the following assumes fixed columns, so inspect!
# -----------------------------------------------------------------------------


# Only want mouse and human for input and TF ChIP experiments

meta_all_clean <- meta_all %>%
  filter(X2 %in% c("mm9", "hg19") &
           X3 %in% c("Input control", "TFs and others"))


# Columns 10 and beyond correspond to submitted user information on the GEO 
# page - variable in presence, naming, and position. Looking for any mention
# of antibody for later curation.

get_antibody <- function(meta, ncores = cores) {
  
  antibody <- mclapply(1:nrow(meta), function(x) {
    match_text <- str_extract(meta[x, 10:ncol(meta)], ".*antibody.*")
    match_text <- paste(na.omit(match_text), collapse = " | ")
    ifelse(match_text == "", NA, match_text)
  }, mc.cores = ncores)

  return(unlist(antibody)) 
}


antibody <- get_antibody(meta_all_clean)


# Split the compressed QC columns (nReads, % mapped", % duplicated, nPeaks)

qc_cols <- str_split_fixed(meta_all_clean[, "X8"], ",", 4)

# Keep main columns and add antibody/QC

meta_all_clean_subset <- meta_all_clean[ , c(1:7, 9)]
meta_all_clean_subset <- cbind(meta_all_clean_subset, unlist(antibody), qc_cols)

# Set names

if (!all(unique(meta_all_clean$X2) %in% c("hg19", "mm9"))) {
  warning("double check column names are in order!")
}


colnames(meta_all_clean_subset) <-
  c(
    "ID",
    "Species",
    "Antigen_Class",
    "Symbol",
    "Cell_Type_Class",
    "Cell_Type",
    "Cell_Type_Description",
    "GEO_Title",
    "Antibody",
    "CA_Read_Count",
    "CA_Perc_Mapped",
    "CA_Perc_Duplicated",
    "CA_Peak_Count"
  )


# Chip-Atlas uses mm9/hg19. Pipeline uses mm10/hg38. Generalizing to species

meta_all_mm <- filter(meta_all_clean_subset, Species == "mm9")
meta_all_mm$Species <- "Mouse"

meta_all_hg <- filter(meta_all_clean_subset, Species == "hg19")
meta_all_hg$Species <- "Human"


# Separate input and ChIP meta

meta_input_mm <- filter(meta_all_mm, Antigen_Class == "Input control")
meta_tf_mm <- filter(meta_all_mm, Antigen_Class == "TFs and others")

meta_input_hg <- filter(meta_all_hg, Antigen_Class == "Input control")
meta_tf_hg <- filter(meta_all_hg, Antigen_Class == "TFs and others")


# Restrict to TFs 

meta_tf_mm <- filter(meta_tf_mm, Symbol %in% tf_mm$Symbol)
meta_tf_hg <- filter(meta_tf_hg, Symbol %in% tf_hg$Symbol)


# Save out metadata

table_list <- list(
  meta_tf_hg,
  meta_input_hg,
  meta_tf_mm,
  meta_input_mm
)

table_names <-  c(
  "chipatlas_hg_tf_meta_",
  "chipatlas_hg_input_meta_",
  "chipatlas_mm_tf_meta_",
  "chipatlas_mm_input_meta_"
)


for (i in 1:length(table_list)) {
  write.table(
    x = table_list[[i]],
    sep = "\t",
    row.names = FALSE,
    na = "NA",
    file = file.path(out_dir, paste0(table_names[i], ca_date, ".tsv"))
  )
}
