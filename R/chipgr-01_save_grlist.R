## Save GRlist objects for the Batch 1 ChIP-seq experiments
## -----------------------------------------------------------------------------

library(tidyverse)
library(GenomicRanges)
source("~/regnetR/R/utils/range_table_functions.R")

date <- "Apr2022"
peakset <- "idr"  # idr for TF (but overlap still for mecp2) or all overlap
min_peaks <- 100
pipeline_dir <- "/cosmos/data/pipeline-output/chipseq-encode-pipeline/chip"
out_dir <- "~/Data/Annotated_objects/GRanges/"

# batch 1 ChIP-seq meta and directories of peak files
run_ids <- read.delim(paste0("~/Data/Metadata/Chipseq/batch1_run_dirs_", date, ".tsv"), stringsAsFactors = FALSE)
meta <- read.delim(paste0("~/Data/Metadata/Chipseq/batch1_chip_meta_final_", date, ".tsv"), stringsAsFactors = FALSE)
stopifnot(all(meta$Experiment_ID %in% run_ids$Experiment_ID))

# ENCODE blacklists
bl_hg <- read.delim("~/Data/Chromosome_info/blacklist_hg38.tsv", stringsAsFactors = FALSE)
bl_mm <- read.delim("~/Data/Chromosome_info/blacklist_mm10.tsv", stringsAsFactors = FALSE)

ids <- left_join(meta, run_ids[, c("Experiment_ID", "Dir")], by = "Experiment_ID")

stopifnot(!any(ids$Dir == "" | is.na(ids$Dir)))

ids_hg <- ids[ids$Species == "Human", ]
ids_mm <- ids[ids$Species == "Mouse", ]


# Load each complete ChIP-seq experiment into a list of GRanges objects
#-------------------------------------------------------------------------------


# Human


bl_gr_hg <- bl_to_gr(bl_hg)


gr_list_hg <- lapply(1:nrow(ids_hg), function(i) {
  
  pipeline_type <- 
    ifelse(str_detect(ids_hg$Experiment_ID[i], "HISTONE"), "overlap", peakset)
  
  peak_table <- read_encpeak(ids_hg$Dir[i], peakset = pipeline_type) %>% 
    get_standard_chr()  # consistent chr names with ensembl table
  
  # Convert to a GRanges object and remove blacklist regions.
  # Format = TRUE will strip start and stop and only return summit
  peak_gr <- peak_to_gr(peak_table, format = FALSE) %>% 
    filter_blacklist(bl_gr = bl_gr_hg)

  message(ids_hg$Experiment_ID[i], " complete ", Sys.time())
  
  return(peak_gr)
  
})
  
names(gr_list_hg) <- ids_hg$Experiment_ID

# Convert to GRList object
gr_list_hg <- GRangesList(gr_list_hg)


# Mouse


bl_gr_mm <- bl_to_gr(bl_mm)


gr_list_mm <- lapply(1:nrow(ids_mm), function(i) {
  
  pipeline_type <- 
    ifelse(str_detect(ids_mm$Experiment_ID[i], "HISTONE"), "overlap", peakset)
  
  peak_table <- read_encpeak(ids_mm$Dir[i], peakset = pipeline_type) %>% 
    get_standard_chr()  # consistent chr names with ensembl table
  
  # Convert to a GRanges object and remove blacklist regions.
  # Format = TRUE will strip start and stop and only return summit
  peak_gr <- peak_to_gr(peak_table, format = FALSE) %>% 
    filter_blacklist(bl_gr = bl_gr_mm)
  
  message(ids_mm$Experiment_ID[i], " complete ", Sys.time())
  
  return(peak_gr)
  
})

names(gr_list_mm) <- ids_mm$Experiment_ID


# Convert to GRList object
gr_list_mm <- GRangesList(gr_list_mm)


# Save 
#-------------------------------------------------------------------------------


saveRDS(
  object = gr_list_hg,
  file = paste0(out_dir, "human_batch1_grlist_peakset=", peakset, "_", date, ".RDS")
)


saveRDS(
  object = gr_list_mm,
  file = paste0(out_dir, "mouse_batch1_grlist_peakset=", peakset, "_", date, ".RDS")
)
