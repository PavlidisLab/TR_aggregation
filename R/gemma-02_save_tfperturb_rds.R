## Save out a list of all curated/successfully loaded batch 1 TF perturbation
## resultsets (after processing) into an RDS object. 
## -----------------------------------------------------------------------------

library(googlesheets4)
library(tidyverse)
library(parallel)
source("R/setup-01_config.R")
source("R/utils/gemma_functions.R")

out_file <- paste0(expr_dir, "TF_perturb_batch1_rslist_", date, ".RDS")
out_file_unfilt <- paste0(expr_dir, "TF_perturb_batch1_unfiltered_rslist_", date, ".RDS")

# curated sheet to link experiments to resultset IDs + relevant columns
rs_df <- read_sheet(ss = gsheets_perturb, 
                    sheet = paste0("Curated_Loaded_resultset_IDs_", date), 
                    trim_ws = TRUE, 
                    col_types = "c")

# meta to match/order saved experiments
meta <- read.delim(paste0(meta_dir, "batch1_tfperturb_meta_final_", date, ".tsv"), stringsAsFactors = FALSE)

rs_df <- rs_df[order(match(rs_df$Experiment_ID, meta$Experiment_ID)), ]

stopifnot(identical(rs_df$Experiment_ID, meta$Experiment_ID))

# protein coding genes
pc_hg <- read.delim(paste0(meta_dir, "refseq_select_hg38.tsv"), stringsAsFactors = FALSE)
pc_mm <- read.delim(paste0(meta_dir, "refseq_select_mm10.tsv"), stringsAsFactors = FALSE)


# load each resultset into a list, making the following adjustments:
# 1) Only keep relevant columns and minimal FoldChange, Tstat, PValue names
# ** In a new list for filtered data:
# 2) Only keep protein coding genes
# 3) Remove probes that map to multiple symbols or no symbols
# 4) If there are multiple probes/entries for a symbol, keep only max abs(tstat)
# 5) Add FDR adjusted pvals
# 6) Add the percentile rank fold change: [0,1] of absolute fold change
# ------------------------------------------------------------------------------


load_rs_list <- function(rs_df, cores) {
  
  l <- mclapply(1:nrow(rs_df), function(i) {
    
    rs <-
      load_result_set(GSE = rs_df$GSE[i], 
                      file = rs_df$Resultset_ID[i],
                      results_dir = rs_dir) %>%
      keep_match_cols(string = rs_df$Match_Col[i]) %>%
      strip_colnames()
    
  }, mc.cores = cores)
  names(l) <- rs_df$Experiment_ID
  
  return(l)
}



prep_rs_list <- function(rs_list, meta, symbol_hg, symbol_mm, cores) {
  
  stopifnot(identical(meta$Experiment_ID, names(rs_list)))
  stopifnot(all(meta$Species %in% c("Human", "Mouse")))
  
  l <- mclapply(names(rs_list), function(x) {
    
    if (filter(meta, Experiment_ID == x)$Species == "Human") {
      symbol <- symbol_hg
    } else  {
      symbol <- symbol_mm
    }
    
    x <- rs_list[[x]] %>%
      filter(Symbol %in% symbol) %>%
      keep_single_symbols() %>%
      filter_by_max_tstat()
    
    x$Adj_pval <- p.adjust(x$PValue, method = "fdr")
    x$PercRankFC <- get_perc_rank_fc(x)
    
    return(x)
    
  }, mc.cores = cores)
  
  names(l) <- meta$Experiment_ID
  return(l)
  
}



if (!all(file.exists(out_file, out_file_unfilt))) {
  
  rs_list_unfilt <- load_rs_list(rs_df, cores = cores)
  
  rs_list_filt <- prep_rs_list(rs_list_unfilt, 
                               meta = meta,
                               symbol_hg = pc_hg$Symbol, 
                               symbol_mm = pc_mm$Symbol,
                               cores = cores)
  
  stopifnot(all(unlist(lapply(rs_list_filt, ncol) == 9)))
  
  saveRDS(rs_list_unfilt, file = out_file_unfilt)
  saveRDS(rs_list_filt, file = out_file)
}
