## Save out a list of all curated/successfully loaded batch 1 TF perturbation
## resultsets (after processing) into an RDS object. 
## -----------------------------------------------------------------------------

library(googlesheets4)
library(tidyverse)
library(parallel)
source("~/regnetR/R/utils/gemma_functions.R")

gsheets_id <- "1oXo1jfoPYcX94bq2F6Bqm1Es1glc8g9mnJvYAO37Vw4"
date <- "Apr2022"
exprs_dir <- "~/Data/Expression_files/Gemma/Resultsets/"
out_file <- paste0("~/Data/Expression_files/Gemma/TF_perturb_batch1_rslist_", date, ".RDS")
out_unfilt_file <- paste0("~/Data/Expression_files/Gemma/TF_perturb_batch1_unfiltered_rslist_", date, ".RDS")

# curated sheet to link experiments to resultset IDs + relevant columns
rs_df <- read_sheet(ss = gsheets_id, 
                    sheet = paste0("Curated_Loaded_resultset_IDs_", date), 
                    trim_ws = TRUE, 
                    col_types = "c")

# meta to match/order saved experiments
meta <- read.delim(paste0("~/Data/Metadata/Gemma/batch1_tfperturb_loaded_meta_", date, ".tsv"), stringsAsFactors = FALSE)

rs_df <- rs_df[order(match(rs_df$Experiment_ID, meta$Experiment_ID)), ]

stopifnot(identical(rs_df$Experiment_ID, meta$Experiment_ID))


# load each resultset into a list, making the following adjustments:
# 1) Only keep relevant columns and minimal FoldChange, Tstat, PValue names
# ** In a new list for filtered data: 
# 2) Remove probes that map to multiple symbols or no symbols
# 3) If there are multiple probes/entries for a symbol, keep only max abs(tstat)
# 4) Add FDR adjusted pvals
# 5) Add the percentile rank fold change: [0,1] of absolute fold change
# ------------------------------------------------------------------------------


load_rs_list <- function(rs_df, cores = 8) {
  
  l <- mclapply(1:nrow(rs_df), function(i) {
    
    rs <-
      load_result_set(GSE = rs_df$GSE[i], file = rs_df$Resultset_ID[i]) %>%
      keep_match_cols(string = rs_df$Match_Col[i]) %>%
      strip_colnames()
    
  }, mc.cores = cores)
  names(l) <- rs_df$Experiment_ID
  
  return(l)
}


prep_rs_list <- function(rs_list, cores = 8) {
  
  l <- mclapply(rs_list, function(x) {
    
    x <- x %>%
      keep_single_symbols() %>%
      filter_by_max_tstat()
    
    x$Adj_pval <- p.adjust(x$PValue, method = "fdr")
    x$PercRankFC <- get_perc_rank_fc(x)
    
    return(x)
    
  }, mc.cores = cores)
  
  return(l)
  
}


if (!all(file.exists(out_file, out_unfilt_file))) {
  
  rs_unfilt_list <- load_rs_list(rs_df)
  rs_filt_list <- prep_rs_list(rs_unfilt_list)
  
  stopifnot(all(unlist(lapply(rs_filt_list, ncol) == 9)))
  
  saveRDS(rs_unfilt_list, file = out_unfilt_file)
  saveRDS(rs_filt_list, file = out_file)
}
