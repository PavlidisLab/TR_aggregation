## Download and organize relevant information from Chip Atlas's ChIP-seq meta
## sheet. Only want to keep TF data sets in mouse and human. The downloaded table
## then undergoes additional manual curation to verify details, match inputs, and
## add additional experiments from other databases.
## -----------------------------------------------------------------------------

library(tidyverse)
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

if (!file.exists(meta_path)) {
  download.file(url, meta_path)
}

assertthat::is.readable(meta_path)

# The table has a lot of empty elements - need to get field count. After loading
# inspect for column names.

nfields <- max(count.fields(meta_path, sep = "\t"), na.rm = TRUE)

all_meta <- read.delim(
  meta_path,
  col.names = 1:nfields,
  quote = "",
  header = FALSE,
  stringsAsFactors = FALSE
)

dim(all_meta) # aug2020: 439593    131

# only want mouse and human for input and TF ChIP experiments

all_meta_clean <- all_meta %>%
  filter(X2 %in% c("mm9", "hg19") &
           X3 %in% c("Input control", "TFs and others"))


# Columns 10 and beyond correspond to submitted user information on the GEO 
# page - variable in presence, naming, and position. Looking for any mention
# of antibody for later curation.


antibody <- parallel::mclapply(1:nrow(all_meta_clean), function(x) {
  match_text <- str_extract(all_meta_clean[x, 10:nfields], ".*antibody.*")
  match_text <- paste(na.omit(match_text), collapse = " | ")
  ifelse(match_text == "", NA, match_text)
}, mc.cores = cores)


# Split the compressed QC columns

qc_cols <- str_split_fixed(all_meta_clean[, "X8"], ",", 4)

# Keep main columns and add antibody/QC

all_meta_clean_subset <- all_meta_clean[ , c(1:7, 9)]
all_meta_clean_subset <- cbind(all_meta_clean_subset, unlist(antibody), qc_cols)

# Set names

if (!all(unique(all_meta_clean$X2) %in% c("hg19", "mm9"))) {
  warning("double check column names are in order!")
}


colnames(all_meta_clean_subset) <-
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

all_meta_mm <- filter(all_meta_clean_subset, Species == "mm9")
all_meta_mm$Species <- "Mouse"

all_meta_hg <- filter(all_meta_clean_subset, Species == "hg19")
all_meta_hg$Species <- "Human"


# Separate input and ChIP meta

input_meta_mm <- filter(all_meta_mm, Antigen_Class == "Input control")
tf_meta_mm <- filter(all_meta_mm, Antigen_Class == "TFs and others")

input_meta_hg <- filter(all_meta_hg, Antigen_Class == "Input control")
tf_meta_hg <- filter(all_meta_hg, Antigen_Class == "TFs and others")


# Restrict to TFs 

tf_meta_mm <- filter(tf_meta_mm, Symbol %in% tf_mm$Symbol)
tf_meta_hg <- filter(tf_meta_hg, Symbol %in% tf_hg$Symbol)


# Save out metadata

table_list <- list(
  tf_meta_hg,
  input_meta_hg,
  tf_meta_mm,
  input_meta_mm
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
