## Subset Chip Atlas metadata to TF of interest, format and export for curation.
## This was structured around 'batch 1' TFs
## Ascl1, Hes1, Mecp2, Mef2c, Neurod1, Pax6, Runx1, Tcf4
## -----------------------------------------------------------------------------

library(tidyverse)
library(RCurl)
library(googlesheets4)
source("R/setup-01_config.R")

# tfs of interest
tfs <- c("pax6", "ascl1", "neurod1", "mecp2", "hes1", "runx1", "mef2c", "tcf4")

# Current Chip atlas meta update and previous, to isolate only new additions
ca_date <- "aug2020"
ca_date_old <- "dec2019"

# Output
meta_dir <- file.path(meta_dir, "Chipseq", "Chipatlas")
outfile_all <- file.path(meta_dir, paste0("batch1_chipmeta_tobecurated_all_", ca_date, ".tsv"))
outfile_new <- file.path(meta_dir, paste0("batch1_chipmeta_tobecurated_new_", ca_date, ".tsv"))

# Loading input metadata to be subset
meta_hg <- read.delim(file.path(meta_dir, paste0("chipatlas_hg_tf_meta_", ca_date, ".tsv")), stringsAsFactors = FALSE)
meta_mm <- read.delim(file.path(meta_dir, paste0("chipatlas_mm_tf_meta_", ca_date, ".tsv")), stringsAsFactors = FALSE)
meta_input_hg <- read.delim(file.path(meta_dir, paste0("chipatlas_hg_input_meta_", ca_date, ".tsv")), stringsAsFactors = FALSE)
meta_input_mm <- read.delim(file.path(meta_dir, paste0("chipatlas_mm_input_meta_", ca_date, ".tsv")), stringsAsFactors = FALSE)

# load older metadata to filter for 'new' experiments only for curation
meta_old <- read.delim(file.path(meta_dir, paste0("batch1_chipmeta_tobecurated_all_", ca_date_old, ".tsv")), stringsAsFactors = FALSE)


# Function to search GEO to match a GSM sample ID to the GSE series ID

get_gse <- function(gsm) {
  page <-  RCurl::getURL(paste0("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=", gsm))
  paste(stringr::str_extract_all(page, "GSE[0-9]*?(?=<)")[[1]], collapse = ", ")
}


# Subset to batch 1 TFs, add GEO accession information, and keep cols of interest
# ------------------------------------------------------------------------------


stopifnot(identical(colnames(meta_hg), colnames(meta_mm)))
stopifnot(identical(colnames(meta_input_hg), colnames(meta_input_mm)))

meta_filter <- filter(rbind(meta_hg, meta_mm), str_to_lower(Symbol) %in% tfs)

input_meta <- rbind(meta_input_hg, meta_input_mm)

# Get the GSM from the GEO_Title field. Some entries have no GSM, or multiple 
# due to curation errors, so coerce to NA

gsm <- str_extract_all(meta_filter$GEO_Title, "^GSM[0-9]+")

gsm <- unlist(lapply(gsm, function(x) {
  if (length(x) != 1) x <- NA 
  return(x)
}))

gse <- unlist(lapply(gsm, get_gse))

# only keeping minimal columns/details for downstream manual curation, and
# create separate columns for GEO identifiers

keep_cols <- c("ID", "Species", "Symbol", "Cell_Type", "GEO_Title", "Antibody")

meta_filter <- meta_filter %>%
  dplyr::select(all_of(keep_cols)) %>%
  mutate(
    GEO_Title = str_replace(GEO_Title, "^GSM[0-9]+: ", ""),
    GSM = gsm,
    GSE = gse
  )


# Look to match input SRX IDs with the TF Chip experiment. Heuristic: extract
# numerical part of SRX ID, add a 'window' and search input SRX IDs for matches
# within this window. this assumes that the deposited input IDs will be in
# 'neighborhood' of the experimental sample.
# ------------------------------------------------------------------------------


# Returns the experiment meta (exp_meta) with cols for SRX and GSM IDs

match_input <- function(exp_meta, input_meta, window_size = 500) {
  
  match_input <- lapply(1:nrow(exp_meta), function(x) {
    
    # Restrict the input to common cell type as query and get integers from the IDs
    filter_input <- filter(input_meta, Cell_Type == exp_meta$Cell_Type[x])
    filter_srx_int <- as.integer(str_extract_all(filter_input$ID, "[:digit:]+", simplify = TRUE))
    query_srx_int <- as.integer(str_extract(exp_meta$ID[x], "[:digit:]+"))
    
    # look for input IDS within integer window of the query ID
    match_ids <- which(
      filter_srx_int <= (query_srx_int + window_size) &
        filter_srx_int >= (query_srx_int - window_size))
    
    filter_input <- filter_input[match_ids, ]
    
    # Return NA if no match, otherwise return all candidates input SRXs as well
    # as the corresponding GSM IDs
    if (length(match_ids) == 0) {
      return(NA)
    } else {
      list(
        gsm = paste(str_extract_all(filter_input$GEO_Title, "GSM[:digit:]+"), collapse = ", "),
        srx = paste(filter_input$ID, collapse = ", ")
      )
    }
  })
  
  exp_meta$Input_GSM <- unlist(lapply(match_input, `[`, "gsm"))
  exp_meta$Input_ID <- unlist(lapply(match_input, `[`, "srx"))
  return(exp_meta)
}


meta_filter <- match_input(meta_filter, input_meta)


# Adding columns to be curated, remove redundant antibody text, reorder cols.
# Replicate: [0/1] Does the sample belong to a set of replicates
# Condition [0/1/2] 0: WT; 1: The TF is directly perturbed; 2: A treatment that
# may not directly perturb the TF but could affect its function
# -----------------------------------------------------------------------------


meta_filter <- meta_filter %>% 
 mutate(
   Replicate = NA,
   Condition = NA,
   Antibody = str_replace(Antibody, "chip antibody=|antibody=", "")
 ) %>% 
  relocate(Antibody, .after = last_col())


# Split for experiments only found in latest update

meta_new <- filter(meta_filter, !(ID %in% meta_old$ID))


# Exporting local copies and for manual curation in google sheets
# ------------------------------------------------------------------------------


sheet_write(
  data = meta_filter,
  ss = gsheets_chip,
  sheet = paste0("Raw_Batch1_all_", ca_date)
)


write.table(
  x = meta_filter,
  sep = "\t",
  row.names = FALSE,
  quote = FALSE,
  file = outfile_all
)


sheet_write(
  data = meta_new,
  ss = gsheets_chip,
  sheet = paste0("Raw_Batch1_newonly_", ca_date)
)


write.table(
  x = meta_new,
  sep = "\t",
  row.names = FALSE,
  quote = FALSE,
  file = outfile_new
)
