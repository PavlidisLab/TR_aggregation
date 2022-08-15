## This script loads a saved RDS list of perturbation resultset tables then 
## builds and saves out gene x experiment effect size matrices containing
## Tstats/FC/PRFC/Pvals/AdjPvals
## -----------------------------------------------------------------------------

library(tidyverse)
library(parallel)
source("~/regnetR/R/utils/gemma_functions.R")

date <- "Apr2022"  # latest data freeze
outdir <- "~/Data/Expression_files/Perturb_matrix/"
meta <- read.delim(file =  paste0("~/Data/Metadata/Gemma/batch1_tfperturb_meta_final_", date, ".tsv"), stringsAsFactors = FALSE)
results <- readRDS(paste0("~/Data/Expression_files/Gemma/TF_perturb_batch1_rslist_", date, ".RDS"))
pc_hg <- read.delim("~/Data/Metadata/refseq_select_hg38.tsv", stringsAsFactors = FALSE)
pc_mm <- read.delim("~/Data/Metadata/refseq_select_mm10.tsv", stringsAsFactors = FALSE)
pc_ortho <- read.delim("~/Data/Metadata/hg_mm_1to1_ortho_genes_DIOPT-v8.tsv", stringsAsFactors = FALSE)

stopifnot(identical(names(results), meta$Experiment_ID))

meta_hg <- filter(meta, Species == "Human")
meta_mm <- filter(meta, Species == "Mouse")

symbol_hg <- unique(pc_hg$Symbol)
symbol_mm <- unique(pc_mm$Symbol)


# Functions
# ------------------------------------------------------------------------------


rm_all_na <- function(mat) {
  # Remove genes/rows for which all observations are NA
  all_na <- which(apply(mat, 1, function(x) sum(is.na(x)) == ncol(mat)))
  if (length(all_na) > 0) mat <- mat[-all_na, ]
  return(mat)
}



get_mat_list <- function(meta, symbol, result_list) {
  
  # Return a list of Tstat/FC/Pval/FDR matrices for experiments contained in meta
  
  tf_mat <- matrix(nrow = length(symbol), ncol = nrow(meta))
  rownames(tf_mat) <- symbol
  colnames(tf_mat) <- meta$Experiment_ID
  
  mat_list <- list(
    Tstat_mat = tf_mat,
    FC_mat = tf_mat,
    PRFC_mat = tf_mat,
    Pval_mat = tf_mat,
    FDR_mat = tf_mat
  )
  
  # for each experiment for the TF, extract the values and fill the matrix row
  # wise with the matched symbols
  
  for (id in meta$Experiment_ID) {
    
    result <- filter(result_list[[id]], Symbol %in% symbol)
    where <- match(result$Symbol, rownames(tf_mat))
    
    mat_list$Tstat_mat[where, id] <- result$Tstat
    mat_list$FC_mat[where, id] <- result$FoldChange
    mat_list$PRFC_mat[where, id] <- result$PercRankFC
    mat_list$Pval_mat[where, id] <- result$PValue
    mat_list$FDR_mat[where, id] <- result$Adj_pval
  
  }
  return(mat_list)
}



mlist_hg <- get_mat_list(meta_hg, symbol_hg, results)
mlist_mm <- get_mat_list(meta_mm, symbol_mm, results)

# Remove genes with no measurements across all studies
mlist_mm <- lapply(mlist_mm, rm_all_na)
mlist_hg <- lapply(mlist_hg, rm_all_na)


# Create a single matrix of 1:1 orthologous genes between mouse and human
# Same procedure as above, just need to map the correct symbol
# ------------------------------------------------------------------------------


# initialize matrix

ortho_mat <- matrix(nrow = nrow(pc_ortho), ncol = nrow(meta))
rownames(ortho_mat) <- pc_ortho$ID
colnames(ortho_mat) <- meta$Experiment_ID

mlist_ortho <- list(
  Tstat_mat = ortho_mat,
  FC_mat = ortho_mat,
  PRFC_mat = ortho_mat,
  Pval_mat = ortho_mat,
  FDR_mat = ortho_mat
)


for (i in 1:nrow(meta)) {
  
  if (meta$Species[i] == "Mouse") {
    pcoding_symbol <- "Symbol_mm"
  } else if (meta$Species[i] == "Human") {
    pcoding_symbol <- "Symbol_hg"
  }
  
  id <- meta$Experiment_ID[i]
  result <- filter(results[[id]], Symbol %in% pc_ortho[, pcoding_symbol])
  result_symbol <- intersect(result$Symbol, pc_ortho[, pcoding_symbol])
  where <- match(result_symbol, pc_ortho[, pcoding_symbol])
  
  mlist_ortho$Tstat_mat[where, id] <- result$Tstat
  mlist_ortho$FC_mat[where, id] <- result$FoldChange
  mlist_ortho$PRFC_mat[where, id] <- result$PercRankFC
  mlist_ortho$Pval_mat[where, id] <- result$PValue
  mlist_ortho$FDR_mat[where, id] <- result$Adj_pval
  
}

mlist_ortho <- lapply(mlist_ortho, rm_all_na)


# Save out
saveRDS(mlist_hg, paste0(outdir, "human_list_perturb_matrix_", date, ".RDS"))
saveRDS(mlist_mm, paste0(outdir, "mouse_list_perturb_matrix_", date, ".RDS"))
saveRDS(mlist_ortho, paste0(outdir, "ortho_list_perturb_matrix_", date, ".RDS"))