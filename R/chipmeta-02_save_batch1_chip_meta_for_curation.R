## Subset the list of experiments from Chip Atlas metadata to only include TFs
## of interest. At time of writing this was structured around 'batch 1' TFs - 
## Ascl1, Hes1, Mecp2, Mef2c, Neurod1, Pax6, Runx1, Tcf4
## This table is then exported to googlesheets for curation (verify context,
## match input controls, tidy). 
## -----------------------------------------------------------------------------

library(tidyverse)
library(RCurl)
library(googlesheets4)

date <- "aug2020"  # latest chip atlas metadata update
date_old <- "dec2019"  # previous update, to isolate only new additions
metadir <- "~/Data/Metadata/Chipseq/Chipatlas/"
gsheets_id <- "1rGVnLL0eXHqr97GM1tloiWwwrJUUmj_ZjW5UOHFN1cc"
outfile_all <- paste0(metadir, "batch1_chipmeta_tobecurated_all_", date, ".tsv")
outfile_new <- paste0(metadir, "batch1_chipmeta_tobecurated_new_", date, ".tsv")

human_meta <- read.delim(paste0(metadir, "chipatlas_hg_tf_meta_", date, ".tsv"), stringsAsFactors = FALSE)
mouse_meta <- read.delim(paste0(metadir, "chipatlas_mm_tf_meta_", date, ".tsv"), stringsAsFactors = FALSE)
human_input_meta <- read.delim(paste0(metadir, "chipatlas_hg_input_meta_", date, ".tsv"), stringsAsFactors = FALSE)
mouse_input_meta <- read.delim(paste0(metadir, "chipatlas_mm_input_meta_", date, ".tsv"), stringsAsFactors = FALSE)

# load older metadata to filter for 'new' experiments only for curation
old_meta <- read.delim(paste0(metadir, "batch1_chipmeta_tobecurated_all_", date_old, ".tsv"))

get_gse <- function(gsm) {
  # Searches GEO to match a GSM sample ID to the GSE series
  page <-  RCurl::getURL(paste0("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=", gsm))
  paste(stringr::str_extract_all(page, "GSE[0-9]*?(?=<)")[[1]], collapse = ", ")
}


# Subset to batch 1 TFs and add GEO accession information
# ------------------------------------------------------------------------------


tfs <- c("pax6", "ascl1", "neurod1", "mecp2", "hes1", "runx1", "mef2c", "tcf4")

stopifnot(identical(colnames(human_meta), colnames(mouse_meta)))
stopifnot(identical(colnames(human_input_meta), colnames(mouse_input_meta)))

filter_meta <- filter(rbind(human_meta, mouse_meta), tolower(Symbol) %in% tfs)

input_meta <- rbind(human_input_meta, mouse_input_meta)

# get the GSM from the GEO_Title field
gsm <- str_extract_all(filter_meta$GEO_Title, "^GSM[0-9]+")

# some entries have no GSM, or multiple due to curation errors. coerce to NA
gsm <- unlist(lapply(gsm, function(x) {
  if (length(x) != 1) x <- NA 
  return (x)
}))

# get the GSE identifier from the GSM id. Slow!
gse <- unlist(lapply(gsm, get_gse))

# only keeping minimal columns/details for downstream manual curation, and
# create seperate columns for GEO identifiers
keep_cols <- c("ID", "Species", "Symbol", "Cell_Type", "GEO_Title", "Antibody")
filter_meta_raw <- filter_meta[, keep_cols]
filter_meta_raw$GEO_Title <- str_replace(filter_meta_raw$GEO_Title, "GSM[0-9]+: ", "")
filter_meta_raw$GSM <- gsm
filter_meta_raw$GSE <- gse


# Look to match input SRX IDs with the TF Chip experiment. Heuristic: extract
# numerical part of SRX ID, add a 'window' and search input SRX IDs for matches
# within this window. this assumes that the deposited input IDs will be in
# 'neighbourhood' of the experimental sample.
# ------------------------------------------------------------------------------


window_size <-  500

match_input <- lapply(1:nrow(filter_meta_raw), function(x) {
  
  # restrict the input to common cell type as query and get integers from the IDs
  filter_input <- filter(input_meta, Cell_Type == filter_meta_raw$Cell_Type[x])
  filter_srx_int <- as.integer(str_extract_all(filter_input$ID, "[:digit:]+", simplify = TRUE))
  query_srx_int <- as.integer(str_extract(filter_meta_raw$ID[x], "[:digit:]+"))

  # look for input IDS within arbitrary window of the query ID
  match_ids <- which(
    filter_srx_int < (query_srx_int + window_size) &
      filter_srx_int > (query_srx_int - window_size))
  
  filter_input <- filter_input[match_ids, ]

  # return NA if no match, otherwise return all candidates input SRXs as well
  # as the corresponding GSM IDs
  if (length(match_ids) == 0) {
    return (NA)
  } else {
    list(
      gsm = paste(str_extract_all(filter_input$GEO_Title, "GSM[:digit:]+"), collapse = ", "),
      srx = paste(filter_input$ID, collapse = ", ")
    )
  }
}
)
names(match_input) <- filter_meta_raw$ID

filter_meta_raw$Input_GSM <- unlist(lapply(match_input, `[`, "gsm"))
filter_meta_raw$Input_ID <- unlist(lapply(match_input, `[`, "srx"))


# Replicate: [0/1] Does the sample belong to a set of replicates
# Condition [0/1/2] 0: WT; 1: The TF is directly perturbed; 2: A treatment that
# may not directly perturb the TF but could affect its function
# -----------------------------------------------------------------------------


filter_meta_raw$Replicate <- rep(NA, nrow(filter_meta_raw))
filter_meta_raw$Condition <- rep(NA, nrow(filter_meta_raw))

# column ordering and remove redundant antibody text

filter_meta_raw <- filter_meta_raw[, !names(filter_meta_raw) == "Antibody"]
filter_meta_raw$Antibody <- filter_meta$Antibody

filter_meta_raw$Antibody <- str_replace(filter_meta_raw$Antibody,
                                        "chip antibody=|antibody=", "")


# Split for experiments only found in latest update

meta_new_only <- filter(filter_meta_raw, !(ID %in% old_meta$ID))


# Exporting for manual curation in google sheets
# ------------------------------------------------------------------------------


sheet_write(
  data = filter_meta_raw,
  ss = gsheets_id,
  sheet = paste0("Raw_Batch1_all_", date)
)


write.table(
  x = filter_meta_raw,
  sep = "\t",
  row.names = FALSE,
  quote = FALSE,
  file = outfile_all
)


sheet_write(
  data = meta_new_only,
  ss = gsheets_id,
  sheet = paste0("Raw_Batch1_newonly_", date)
)


write.table(
  x = meta_new_only,
  sep = "\t",
  row.names = FALSE,
  quote = FALSE,
  file = outfile_new
)
