## Download and organize relevant information from Chip Atlas's ChIP-seq meta
## sheet. Only want to keep TF data sets in mouse and human. The downloaded table
## then undergoes additional manual curation to verify details, match inputs, and
## add additional experiments from other databases.

## Be aware of date of latest update! https://github.com/inutano/chip-atlas/wiki
## Last update: August 2020 (checked May2022, still same)
## -----------------------------------------------------------------------------

library(tidyverse)

date <- "aug2020"

meta_name <- paste0("chipatlas_all_meta_", date, ".tsv")
outdir <- "~/Data/Metadata/Chipseq/Chipatlas/"
meta_path <- paste0(outdir, meta_name)
url <- "http://dbarchive.biosciencedbc.jp/kyushu-u/metadata/experimentList.tab"

mouse_tfs <- read.delim("~/Data/Metadata/mouse_tfs.tsv", stringsAsFactors = FALSE)
human_tfs <- read.delim("~/Data/Metadata/human_tfs.tsv", stringsAsFactors = FALSE)

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

dim (all_meta) # aug2020: 439593    131

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
}, mc.cores = 8)


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
mouse_all_meta <- filter(all_meta_clean_subset, Species == "mm9")
mouse_all_meta$Species <- "Mouse"

human_all_meta <- filter(all_meta_clean_subset, Species == "hg19")
human_all_meta$Species <- "Human"


# Separate input and ChIP meta
mouse_input_meta <- filter(mouse_all_meta, Antigen_Class == "Input control")
mouse_tf_meta <- filter(mouse_all_meta, Antigen_Class == "TFs and others")

human_input_meta <- filter(human_all_meta, Antigen_Class == "Input control")
human_tf_meta <- filter(human_all_meta, Antigen_Class == "TFs and others")


# Restrict to known TFs
mouse_tf_meta <- filter(mouse_tf_meta, Symbol %in% mouse_tfs$Symbol)
human_tf_meta <- filter(human_tf_meta, Symbol %in% human_tfs$Symbol)

# save out metadata

table_list <- list(
  human_tf_meta,
  human_input_meta,
  mouse_tf_meta,
  mouse_input_meta
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
    file = paste0(outdir, table_names[i], date, ".tsv")
  )
}
