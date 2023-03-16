## Loads and organizes the ChIP-seq and perturbation data matrices into a single
## list for export, and creates a single data frame for each TR/species of the 
## summarized/aggregated target rankings for export.
## -----------------------------------------------------------------------------

library(plyr)
library(tidyverse)
source("R/setup-01_config.R")
source("R/utils/ranking_functions.R")

# Loading perturb data
perturb_meta <- read.delim(perturb_meta_path, stringsAsFactors = FALSE)
perturb_hg <- readRDS(pmat_path_hg)
perturb_mm <- readRDS(pmat_path_mm)
perturb_ortho <- readRDS(pmat_path_ortho)
tf_de <- readRDS(prank_path_group)

# Loading ChIP-seq data
chip_meta <- read.delim(chip_meta_path, stringsAsFactors = FALSE)
chip_hg <- readRDS(blist_hg_path)
chip_mm <- readRDS(blist_mm_path)
chip_ortho <- readRDS(blist_ortho_path)
bind_summary <- readRDS(bind_summary_path)

# low throughput resource of curated targets
lt_curated <- read.delim(curated_path_all, stringsAsFactors = FALSE)

# Map of orthologous genes 
pc_ortho <- read.delim(ortho_path, stringsAsFactors = FALSE)


# Each data type has its own global filtering steps. For simplicity, only 
# keeping genes that survived filtering across both. First describe genes 
# present in one method but not the other (and thus will be removed)
# ------------------------------------------------------------------------------


chip_not_perturb <- list(
  Human = setdiff(rownames(chip_hg$QN_log), 
                  rownames(perturb_hg$FC_mat)),
  Mouse = setdiff(rownames(chip_mm$QN_log), 
                  rownames(perturb_mm$FC_mat))
)


perturb_not_chip <- list(
  Human = setdiff(rownames(perturb_hg$FC_mat),
                  rownames(chip_hg$QN_log)),
  Mouse = setdiff(rownames(perturb_mm$FC_mat),
                  rownames(chip_mm$QN_log))
)


# Find many more genes (3200+) in mouse perturb relative to mouse chip. Olf or 
# Gm* are common.


n_chip_only <- lapply(chip_not_perturb, length)
n_perturb_only <- lapply(perturb_not_chip, length)


# Find that the genes not in the chip matrix tend to have a lot of NAs in the
# perturb matrix (median ~75% in both species). However, there are examples of
# genes that are frequently measured and have DE counts. CPS1 in human is
# measured in 71/77 experiments and DE in 18/71, but has been filtered in chip
# mat. Similarly, Slc25a37 in mouse has been measured in all 165 experiments and
# is DE in 17/165, but was filtered from chip mat due to low counts.


prop_na <- function(mat, meta, species) {
  
  apply(mat, 1, function(x) {
    sum(is.na(x)) / sum(meta$Species == species)
  }) 
}


na_list <- list(
  Human = prop_na(perturb_hg$FC_mat[perturb_not_chip$Human, ], perturb_meta, "Human"),
  Mouse = prop_na(perturb_mm$FC_mat[perturb_not_chip$Mouse, ], perturb_meta, "Mouse")
)

na_summ <- lapply(na_list, summary)
na_min <- lapply(na_list, function(x) names(x[x == min(x)]))
min_de_hg <- rowSums(perturb_hg$FDR_mat[na_min$Human, ] < 0.05, na.rm = TRUE)
min_de_mm <- rowSums(perturb_mm$FDR_mat[na_min$Mouse, ] < 0.05, na.rm = TRUE)


# Genes not in the perturb matrix have lower scores across experiments, but 
# there are examples of these genes ranking highly in individual experiments.
# KRTAP4-16 high in set of human ASCL1 experiments (note already see keratin 
# cluster genes high in ASCL1). H2ac18 high in mouse Runx1.


chip_only_mean <- list(
  Human = rowMeans(chip_hg$QN_log[chip_not_perturb$Human,]),
  Mouse = rowMeans(chip_mm$QN_log[chip_not_perturb$Mouse,])
)

mean_summ <- lapply(chip_only_mean, summary)

chip_only_rank <- list(
  Human = apply(-chip_hg$QN_log, 2, rank)[chip_not_perturb$Human, ],
  Mouse = apply(-chip_mm$QN_log, 2, rank)[chip_not_perturb$Mouse, ]
)

rank_summ <- lapply(chip_only_rank, function(x) summary(t(x)))

rank_head <- lapply(chip_only_rank, function(x) {
  apply(x, 1, function(y) head(sort(y)))
})


# head(sort(chip_hg$QN_log["KRTAP4-16", ], decreasing = TRUE), 10)
# head(sort(chip_mm$QN_log["H2ac18", ], decreasing = TRUE), 10)


# Effect matrices filtered for common genes into a list
# ------------------------------------------------------------------------------


filter_common <- function(mat_list, symbols) {
  lapply(mat_list, function(x) x[symbols, ])
}


common <- list(
  Human = intersect(rownames(chip_hg$QN_log), rownames(perturb_hg$FC_mat)),
  Mouse = intersect(rownames(chip_mm$QN_log), rownames(perturb_mm$FC_mat)),
  Ortho = intersect(rownames(chip_ortho$QN_log), rownames(perturb_ortho$FC_mat))
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


# Get a list of dataframes for each TR/species that merges the aggregated 
# evidence from ChIP-seq and perturbation experiments, and append the curation
# status for the given gene, and the aggregated + integrated rankings. 
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
    
    curated_df <- curated_df %>% 
      mutate(TF_Symbol = str_to_upper(TF_Symbol)) %>% 
      filter(TF_Symbol == tf)
    
    hg <- filter(pc_ortho, Symbol_hg %in% curated_df$Target_Symbol)$ID
    mm <- filter(pc_ortho, Symbol_mm %in% curated_df$Target_Symbol)$ID
    
    targets <- union(hg, mm)
  }
  
  return(targets)
}


# For a single TF, merge the perturbation/diff expr table, the bind fit table,
# and the bind summary table, then append presence of curated targets, and
# finally add the ranks of the lines of evidence.

# NOTE: A gene can survive the global perturbation filtering, but have all NAs
# for a given TR-specific subset. The 'na_perturb' argument determines if the 
# integrated ranking considers this information: if TRUE, genes with all perturb
# NAs will be tied with the worst integrated ranking (even if the binding 
# ranking is strong). If FALSE, the integrated ranking can still assign elevated
# importance to all perturb NA genes if the corresponding binding evidence is 
# still strong.

merge_tf <- function(de_table,
                     bind_fit = NULL, 
                     bind_summary,
                     curated_df,
                     species,
                     tf,
                     na_perturb = FALSE) {
  
  common_genes <- intersect(de_table$Symbol, bind_summary$Symbol)
  de_table <- filter(de_table, Symbol %in% common_genes)
  bind_summary <- filter(bind_summary, Symbol %in% common_genes)

  
  if (!is.null(bind_fit)) {  # ortho doesn't have bind fit
    
    bind_fit <- prepare_bind_fit(bind_fit)
    
    merge_df <- plyr::join_all(list(de_table, bind_summary, bind_fit), 
                               by = "Symbol")
  } else {
    
    merge_df <- plyr::join_all(list(de_table, bind_summary), by = "Symbol")
    
  }
  
  targets <- get_targets(curated_df = curated_df, species = species, tf = tf)
  
  merge_df <- merge_df %>% 
    mutate(Curated_target = Symbol %in% targets) %>% 
    add_ranks(na_perturb = na_perturb)
  
  return(merge_df)
}


# Applies merge_tf over list of tables each named with specific TR

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

saveRDS(all_dat, file = alldat_path)
saveRDS(rank_list, file = rank_path)
