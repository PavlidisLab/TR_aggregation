## Save ChIP-seq experiments as list of GenomicRange objects
## -----------------------------------------------------------------------------

library(tidyverse)
library(GenomicRanges)
source("R/setup-01_config.R")
source("R/utils/range_table_functions.R")

# batch 1 ChIP-seq meta and directories of peak files
run_ids <- read.delim(chip_run_path, stringsAsFactors = FALSE)
meta <- read.delim(chip_meta_path, stringsAsFactors = FALSE)
stopifnot(all(meta$Experiment_ID %in% run_ids$Experiment_ID))

# Load ENCODE blacklists and convert to GRanges
bl_hg <- bl_to_gr(read.delim(bl_path_hg, stringsAsFactors = FALSE))
bl_mm <- bl_to_gr(read.delim(bl_path_mm, stringsAsFactors = FALSE))

ids <- left_join(meta, run_ids[, c("Experiment_ID", "Dir")], by = "Experiment_ID")

stopifnot(!any(ids$Dir == "" | is.na(ids$Dir)))

ids_hg <- ids[ids$Species == "Human", ]
ids_mm <- ids[ids$Species == "Mouse", ]


# Load each complete ChIP-seq experiment into a list of GRanges objects
#-------------------------------------------------------------------------------


load_gr_list <- function(id_df, bl_gr, peakset = "idr") {
  
  stopifnot(peakset %in% c("idr", "overlap"))
  
  gr_list <- lapply(1:nrow(id_df), function(i) {
    
    # Mecp2 ran on histone pipeline which produces overlap, not idr, peak sets
    pipeline_type <- 
      ifelse(str_detect(id_df$Experiment_ID[i], "HISTONE"), "overlap", peakset)
    
    # load peak table and coerce consistent chr names with protein coding table
    peak_table <- read_encpeak(id_df$Dir[i], peakset = pipeline_type) %>% 
      get_standard_chr()  
    
    # Convert to a GRanges object and remove blacklist regions.
    # Format = TRUE will strip start and stop and only return summit
    peak_gr <- peak_to_gr(peak_table, format = FALSE) %>% 
      filter_blacklist(bl_gr = bl_gr)
    
    message(id_df$Experiment_ID[i], " complete ", Sys.time())
    
    return(peak_gr)
    
  })
  
  names(gr_list) <- id_df$Experiment_ID
  
  # Convert to GRList object
  gr_list <- GRangesList(gr_list)
  
  return(gr_list)
  
}


# Human
gr_hg <- load_gr_list(id_df = ids_hg, bl_gr = bl_hg)

# Mouse
gr_mm <- load_gr_list(id_df = ids_mm, bl_gr = bl_mm)


# Save 
#-------------------------------------------------------------------------------


saveRDS(
  object = gr_hg,
  file = grlist_hg_path
)


saveRDS(
  object = gr_mm,
  file = grlist_mm_path
)
