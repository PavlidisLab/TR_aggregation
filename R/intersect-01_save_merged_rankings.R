## Loads and organizes the ChIP-seq and perturbation data matrices, and create
## a single data frame for each TR/species of the summarized/aggregated target
## rankings for export
## -----------------------------------------------------------------------------

library(plyr)
library(tidyverse)
library(data.table)  # frank() for rank using two variables

date <- "Apr2022"  # most recent data freeze
rank_outfile <- paste0("~/scratch/R_objects/", date, "_ranked_target_list.RDS")
dat_outfile <- paste0("~/scratch/R_objects/", date, "_all_data_list.RDS")

# Loading perturb data
fdr <- 0.1
perturb_meta <- read.delim(paste0("~/Data/Metadata/Gemma/batch1_tfperturb_meta_final_", date, ".tsv"), stringsAsFactors = FALSE)
perturb_hg <- readRDS(paste0("~/Data/Expression_files/Perturb_matrix/human_list_perturb_matrix_", date, ".RDS"))
perturb_mm <- readRDS(paste0("~/Data/Expression_files/Perturb_matrix/mouse_list_perturb_matrix_", date, ".RDS"))
perturb_ortho <- readRDS(paste0("~/Data/Expression_files/Perturb_matrix/ortho_list_perturb_matrix_", date, ".RDS"))
tf_de <- readRDS(paste0("~/Data/TF_perturb/TF_FDR=0.1_counts_list_", date, ".RDS"))

# Loading ChIP-seq data
intensity_flag <- FALSE  # use ouyang scores with or without macs2 intensity
binary_dist <- 25e3  # distance threshold used for binary gene scores
min_peaks <- 100  # how many peaks were required to keep an experiment
peakset <- "idr"  # idr for TF (but overlap still for mecp2) or all overlap
chip_dir <- "~/Data/Annotated_objects/Bind_matrices/Encpipe/" # bind matrix path
chip_type <- "QN_log"  # which chip processing scheme to use
chip_meta <- read.delim(paste0("~/Data/Metadata/Chipseq/batch1_chip_meta_final_", date, ".tsv"), stringsAsFactors = FALSE)
chip_hg <- readRDS(paste0(chip_dir, "Human_refseq_", date, "_processed_bindmat_list_minpeak=", min_peaks, "_ouyang_dc=5000_intensity=", intensity_flag, "_binary=", binary_dist/1e3, "kb_peakset=", peakset, ".RDS"))
chip_mm <- readRDS(paste0(chip_dir, "Mouse_refseq_", date, "_processed_bindmat_list_minpeak=", min_peaks, "_ouyang_dc=5000_intensity=", intensity_flag, "_binary=", binary_dist/1e3, "kb_peakset=", peakset, ".RDS"))
chip_ortho <- readRDS(paste0(chip_dir, "Ortho_refseq_", date, "_processed_bindmat_list_minpeak=", min_peaks, "_ouyang_dc=5000_intensity=", intensity_flag, "_binary=", binary_dist/1e3, "kb_peakset=", peakset, ".RDS"))
bind_summary <- readRDS(paste0("~/scratch/R_objects/", date, "_refseq_bind_summary_peakset=", peakset, ".RDS"))

# low throughput resource of curated targets
lt_curated <- read.delim("~/Data/Metadata/Curated_targets_all_July2022.tsv", stringsAsFactors = FALSE)

# Map of orthologous genes 
pc_ortho <- read.delim("~/Data/Metadata/hg_mm_1to1_ortho_genes_DIOPT-v8.tsv", stringsAsFactors = FALSE)


# Organize and export data matrices
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


# For binding include raw binding scores, quantile norm of log (used for all
# analysis), and the binary assignment

keep_norm <- c("Raw", "QN_log", "Binary")

bind_list <- list(
  Human = filter_common(chip_hg[keep_norm], common$Human),
  Mouse = filter_common(chip_mm[keep_norm], common$Mouse),
  Ortho = filter_common(chip_ortho[keep_norm], common$Ortho),
  Meta = chip_meta
)


all_dat <- list(Binding = bind_list, Perturbation = perturb_list)

saveRDS(all_dat, file = dat_outfile)


# Join each aggregated TR table such that each species/TR combo has a single df 
# of evidence, and add rankings
# ------------------------------------------------------------------------------


prepare_bind_fit <- function(bind_fit) {
  
  # Formats the limma summary table in preparation for joining
  
  bind_fit <- bind_fit %>% 
    rownames_to_column(var = "Symbol") %>% 
    dplyr::rename(Bind_logFC = logFC,
                  Bind_adj_Pval = adj.P.Val) %>% 
    dplyr::select(Symbol, Bind_logFC, Bind_adj_Pval)
  return(bind_fit)
}



get_targets <- function(curated_df, 
                        tf, 
                        species,  # Human|Mouse|Ortho
                        ortho_genes = pc_ortho) {
  
  # Extracts the curated targets for the given TF, specifying species to handle
  # case matching or the orthologous gene name structure [Human_symbol]_[Mouse_symbol]
  
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



add_ranks <- function(df, FC_tiebreak = TRUE) {
  
  # Append the ranks of each gene for the perturbation evidence (count DE) and
  # the ChIP-seq binding signal (mean bind score). If FC_truebreak is TRUE,
  # then DE count ties are broken by the average absolute FC of the gene. Then
  # add the rank product of these two lines of evidence, as in the BETA algo.
  # Because "top ranked" genes from the RP have arbitrary small units, take
  # the rank of the rank product to yield a final integer rank where rank=1
  # has the strongest evidence.
  # BETA: https://pubmed.ncbi.nlm.nih.gov/24263090/
  # R application: https://pubmed.ncbi.nlm.nih.gov/32894066/
  
  if (FC_tiebreak) {
    df$Rank_perturbation <- data.table::frank(
      list(-df$Count_DE, -df$Avg_abs_FC), ties.method = "min")
  } else {
    df$Rank_perturbation <- rank(-df$Count_DE, ties = "min")
  }
  
  df$Rank_binding <- rank(-df$Mean_bind)
  rank_prod <- df$Rank_binding/nrow(df) * df$Rank_perturbation/nrow(df)
  df$Rank_integrated <- rank(rank_prod, ties = "min")
  
  return(arrange(df, Rank_integrated))
}



merge_tf <- function(de_table,
                     bind_fit = NULL, 
                     bind_summary,
                     curated_df = lt_curated,
                     species,
                     tf) {
  
  # For a single TF, merge the perturbation table, the bind fit table, and the
  # bind summary table, then append presence of curated targets, and finally
  # add the ranks of the lines of evidence
  
  targets <- get_targets(curated_df = curated_df, species = species, tf = tf)
  
  if (!is.null(bind_fit)) {
    bind_fit <- prepare_bind_fit(bind_fit)
    merge_df <- plyr::join_all(list(de_table, bind_summary, bind_fit), by = "Symbol")
  } else {
    merge_df <- plyr::join_all(list(de_table, bind_summary), by = "Symbol")
  }
  
  merge_df$Curated_target <- merge_df$Symbol %in% targets
  merge_df <- add_ranks(merge_df)
  
  return(merge_df)
}



merge_all <- function(de_list,
                      bind_fit_list = NULL,
                      bind_summary_list,
                      curated_df = lt_curated,
                      species) {
  
  # Applies merge_tf over list of tables
  
  
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



all_list <- list(
  
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

saveRDS(all_list, rank_outfile)
