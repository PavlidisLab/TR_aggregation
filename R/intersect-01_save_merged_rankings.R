## Loads and organizes the ChIP-seq and perturbation data matrices into a single
## list for export, and creates a single data frame for each TR/species of the 
# summarized/aggregated target rankings for export.
## -----------------------------------------------------------------------------

library(plyr)
library(tidyverse)
library(data.table)  # frank() for rank using two variables
source("R/setup-01_config.R")

# rank for the final ranking dfs by TR 
# dat collects all the contributing summarized data into one list
rank_outfile <- paste0(scratch_dir, date, "_ranked_target_list.RDS")
dat_outfile <- paste0(scratch_dir, date, "_all_data_list.RDS")

# Loading perturb data
perturb_meta <- read.delim(paste0(meta_dir, "Gemma/batch1_tfperturb_meta_final_", date, ".tsv"), stringsAsFactors = FALSE)
perturb_hg <- readRDS(paste0(pmat_dir, "human_list_perturb_matrix_", date, ".RDS"))
perturb_mm <- readRDS(paste0(pmat_dir, "mouse_list_perturb_matrix_", date, ".RDS"))
perturb_ortho <- readRDS(paste0(pmat_dir, "ortho_list_perturb_matrix_", date, ".RDS"))
tf_de <- readRDS(paste0(expr_dir, "TF_perturb_DE_counts_list_by_TF_FDR01_", date, ".RDS"))

# Loading ChIP-seq data
chip_type <- "QN_log"  # which chip processing scheme to use
chip_meta <- read.delim(paste0(meta_dir, "Chipseq/batch1_chip_meta_final_", date, ".tsv"), stringsAsFactors = FALSE)
chip_hg <- readRDS(paste0(cmat_dir, "Human_refseq_", date, "_processed_bindmat_list_minpeak=", min_peaks, "_ouyang_dc=5000_intensity=FALSE_binary=", binary_dist/1e3, "kb.RDS"))
chip_mm <- readRDS(paste0(cmat_dir, "Mouse_refseq_", date, "_processed_bindmat_list_minpeak=", min_peaks, "_ouyang_dc=5000_intensity=FALSE_binary=", binary_dist/1e3, "kb.RDS"))
chip_ortho <- readRDS(paste0(cmat_dir, "Ortho_refseq_", date, "_processed_bindmat_list_minpeak=", min_peaks, "_ouyang_dc=5000_intensity=FALSE_binary=", binary_dist/1e3, "kb.RDS"))
bind_summary <- readRDS(paste0(scratch_dir, date, "_refseq_bind_summary.RDS"))

# low throughput resource of curated targets
lt_curated <- read.delim(paste0(meta_dir, "Curated_targets_all_July2022.tsv"), stringsAsFactors = FALSE)

# Map of orthologous genes 
pc_ortho <- read.delim(paste0(meta_dir, "hg_mm_1to1_ortho_genes_DIOPT-v8.tsv"), stringsAsFactors = FALSE)


# Organize all effect size matrices and metadata tables into a single list
# ------------------------------------------------------------------------------


# filtering rules for each data type may lead to unequal coverage - keep common


filter_common <- function(mat_list, symbols) {
  lapply(mat_list, function(x) x[symbols, ])
}


common <- list(
  Human = intersect(rownames(chip_hg[[chip_type]]), rownames(perturb_hg$FC_mat)),
  Mouse = intersect(rownames(chip_mm[[chip_type]]), rownames(perturb_mm$FC_mat)),
  Ortho = intersect(rownames(chip_ortho[[chip_type]]), rownames(perturb_ortho$FC_mat))
)


perturb_list <- list(
  Human = filter_common(perturb_hg, common$Human),
  Mouse = filter_common(perturb_mm, common$Mouse),
  Ortho = filter_common(perturb_ortho, common$Ortho),
  Meta = perturb_meta
)


# For binding include raw binding scores, quantile norm of log10(mat+1) (used 
# for most analysis), and the binary assignment


keep_norm <- c("Raw", "QN_log", "Binary")

bind_list <- list(
  Human = filter_common(chip_hg[keep_norm], common$Human),
  Mouse = filter_common(chip_mm[keep_norm], common$Mouse),
  Ortho = filter_common(chip_ortho[keep_norm], common$Ortho),
  Meta = chip_meta
)


all_dat <- list(Binding = bind_list, Perturbation = perturb_list)


# Join summarized TR tables from ChIP-seq and perturbation experiments so that 
# each TR/species has a single dataframe of evidence per gene.
# ChIP-seq: mean binding scores, proportion binary, binding specificity model output 
# Perturbation experiments: DE counts, Average absolute FC, purity, DE prior
# Then, add column for whether the given TR-gene pair has curated evidence.
# Finally, add integer ranking of the aggregated ChIP-seq evidence (mean binding
# score), aggregated perturbation (count DE, ties broken by average abs FC), and
# an integrated ranking that takes the rank product between these two lists.
# ------------------------------------------------------------------------------


# Formats the limma summary table in preparation for joining

prepare_bind_fit <- function(bind_fit) {
  
  bind_fit <- bind_fit %>% 
    rownames_to_column(var = "Symbol") %>% 
    dplyr::rename(Bind_logFC = logFC,
                  Bind_adj_Pval = adj.P.Val) %>% 
    dplyr::select(Symbol, Bind_logFC, Bind_adj_Pval)
  return(bind_fit)
}


# Extracts the curated targets for the given TF, specifying species to handle
# case matching or the orthologous gene name structure [Human_symbol]_[Mouse_symbol]

get_targets <- function(curated_df, 
                        tf, 
                        species,  # Human|Mouse|Ortho
                        ortho_genes = pc_ortho) {
  
  if (species == "Human") {
    targets <- curated_df %>%
      mutate(TF_Symbol = str_to_upper(TF_Symbol),
             Target_Symbol = str_to_upper(Target_Symbol)) %>%
      filter(TF_Symbol == tf) %>% 
      distinct(Target_Symbol) %>% 
      pull(Target_Symbol)
    
  } else if (species == "Mouse") {
    targets <- curated_df %>%
      mutate(TF_Symbol = str_to_title(TF_Symbol),
             Target_Symbol = str_to_title(Target_Symbol)) %>%
      filter(TF_Symbol == tf) %>% 
      distinct(Target_Symbol) %>% 
      pull(Target_Symbol)
    
  } else if (species == "Ortho") {
    curated_df <- mutate(curated_df, TF_Symbol = str_to_title(TF_Symbol)) %>% 
      filter(TF_Symbol == tf)
    hg <- filter(pc_ortho, Symbol_hg %in% curated_df$Target_Symbol)$ID
    mm <- filter(pc_ortho, Symbol_mm %in% curated_df$Target_Symbol)$ID
    targets <- union(hg, mm)
  }
  
  return(targets)
}


# Append the ranks of each gene for the perturbation evidence (count DE) and
# the ChIP-seq binding signal (mean bind score). If FC_truebreak is TRUE,
# then DE count ties are broken by the average absolute FC of the gene. Then
# add the rank product of these two lines of evidence, as in the BETA algo.
# Because "top ranked" genes from the RP have arbitrary small units, take
# the rank of the rank product to yield a final integer rank where rank=1
# has the strongest evidence.
# BETA: https://pubmed.ncbi.nlm.nih.gov/24263090/
# R application: https://pubmed.ncbi.nlm.nih.gov/32894066/

add_ranks <- function(df, FC_tiebreak = TRUE) {
  
  if (FC_tiebreak) {
    
    df$Rank_perturbation <- data.table::frank(
      list(-df$Count_DE, -df$Avg_abs_FC), ties.method = "min")
  
  } else {
    
    df$Rank_perturbation <- rank(-df$Count_DE, ties.method = "min")
  
  }
  
  df$Rank_binding <- rank(-df$Mean_bind)
  rank_prod <- df$Rank_binding/nrow(df) * df$Rank_perturbation/nrow(df)
  df$Rank_integrated <- rank(rank_prod, ties.method = "min")
  
  return(arrange(df, Rank_integrated))
}



# For a single TF, merge the perturbation table, the bind fit table, and the
# bind summary table, then append presence of curated targets, and finally
# add the ranks of the lines of evidence

merge_tf <- function(de_table,
                     bind_fit = NULL, 
                     bind_summary,
                     curated_df = lt_curated,
                     species,
                     tf) {
  
  targets <- get_targets(curated_df = curated_df, species = species, tf = tf)
  
  if (!is.null(bind_fit)) {  # ortho doesn't have bind fit
    
    bind_fit <- prepare_bind_fit(bind_fit)
    merge_df <- plyr::join_all(list(de_table, bind_summary, bind_fit), by = "Symbol")
  
  } else {
    
    merge_df <- plyr::join_all(list(de_table, bind_summary), by = "Symbol")
  
  }
  
  merge_df$Curated_target <- merge_df$Symbol %in% targets
  merge_df <- add_ranks(merge_df)
  
  return(merge_df)
}



# Applies merge_tf over list of tables

merge_all <- function(de_list,
                      bind_fit_list = NULL,
                      bind_summary_list,
                      curated_df = lt_curated,
                      species) {
  
  tfs <- intersect(names(de_list), names(bind_summary_list))
  
  tf_list <- lapply(tfs, function(x) {
    
    merge_tf(de_table = de_list[[x]],
             bind_fit = bind_fit_list[[x]], 
             bind_summary = bind_summary_list[[x]],
             curated_df = lt_curated,
             species = species,
             tf = x)
  })
  names(tf_list) <- tfs
  
  return(tf_list)
}



rank_list <- list(
  
  Human = merge_all(de_list = tf_de$Human,
                    bind_fit_list = bind_summary$Human_block_fit,
                    bind_summary_list = bind_summary$Human_bind_tf,
                    species = "Human"),
  
  Mouse = merge_all(de_list = tf_de$Mouse,
                    bind_fit_list = bind_summary$Mouse_block_fit,
                    bind_summary_list = bind_summary$Mouse_bind_tf,
                    species = "Mouse"),
  
  Ortho = merge_all(de_list = tf_de$Ortho,
                    bind_summary_list = bind_summary$Ortho_bind_tf,
                    species = "Ortho")
)


# Save out

saveRDS(all_dat, file = dat_outfile)
saveRDS(rank_list, rank_outfile)
