## Functions for working with perturbation data formatted as a gene x experiment
## matrix of effect sizes
## -----------------------------------------------------------------------------

library(tidyverse)


# Given an gene x experiment FDR mat, return a dataframe summarizing DE status
# across experiments.
# Count_DE: Tally of gene being diff expr at FDR cutoff across experiments
# Proportion_DE_measured: Prop. of DE across experiments with non-NAs for gene
# Proportion_DE_all: Prop. of DE across all experiments
# Count_NA: Tally of NA measurements for the gene across experiments

count_table <- function(fdr_mat, fdr) {
  
  deg_mat <- fdr_mat < fdr
  count_de <- rowSums(deg_mat, na.rm = TRUE)
  count_na <- apply(deg_mat, 1, function(x) sum(is.na(x)))
  count_exps <- ncol(deg_mat)
  prop_measured <- round(count_de / (count_exps - count_na), 3)
  prop_all <- round(count_de / count_exps, 3)
  
  df <- data.frame(
    Symbol = names(count_de),
    Count_DE = count_de,
    Proportion_DE_measured = prop_measured,
    Proportion_DE_all = prop_all,
    Count_NA = count_na,
    stringsAsFactors = FALSE
  )
  
  return(df[order(df$Count_DE, decreasing = TRUE), ])
}



# Return a data frame that has for each gene the count of times it was down/up
# in GoF (overexpression) and LoF (KO + KD) experiments
# GoF_down: Tally when a gene had FC < 0 in overexpression experiments
# GoF_up: Tally when a gene had FC > 0 in overexpression experiments
# LoF_down: Tally when a gene had FC < 0 in KO and KD experiments
# LoF_up: Tally when a gene had FC > 0 in KO and KD experiments

tally_direction <- function(fc_mat, meta) {
  
  gof_exps <- filter(meta, Perturbation == "Overexpression")$Experiment_ID
  lof_exps <- filter(meta, Perturbation %in% c("Knockdown", "Knockout"))$Experiment_ID
  gof_mat <- fc_mat[, gof_exps, drop = FALSE]
  lof_mat <- fc_mat[, lof_exps, drop = FALSE]
  
  gof_down <- apply(gof_mat, 1, function(x) sum(x < 0, na.rm = TRUE))
  gof_up <- apply(gof_mat, 1, function(x) sum(x > 0, na.rm = TRUE))
  lof_down <- apply(lof_mat, 1, function(x) sum(x < 0, na.rm = TRUE))
  lof_up <- apply(lof_mat, 1, function(x) sum(x > 0, na.rm = TRUE))
  
  return(data.frame(
    Symbol = rownames(fc_mat),
    GoF_down = gof_down,
    GoF_up = gof_up,
    LoF_down = lof_down,
    LoF_up = lof_up
  ))
}



# Return a logical vector of the FC consistency across perturbation classes.
# A gene is "consistent across" perturbation types if its predominant change
# of direction (pos or neg FC) across experiments is opposite for gain
# of function (overexpression) and loss of function (knockdowns, knockouts)
# experiments. In other words, a gene is consistent if it predominantly has 
# positive FC in gain of functions and predominantly negative FC in loss of
# functions, or vice versa. It is not consistent if the predominant 
# change of direction (positive or negative FC) is the same across both 
# perturbation classes.

is_consistent_across <- function(count_table) {
  
  stopifnot(
    c("GoF_down", "GoF_up", "LoF_down", "LoF_up") %in% colnames(count_table))
  
  gof_class <- ifelse(
    count_table$GoF_down > count_table$GoF_up, "GoF_down", "GoF_up")
  
  lof_class <- ifelse(
    count_table$LoF_down > count_table$LoF_up, "LoF_down", "LoF_up")
  
  consistent <- ifelse(
    (gof_class == "GoF_down" & lof_class == "LoF_down") |
    (gof_class == "GoF_up" & lof_class == "LoF_up"),
    FALSE,
    TRUE)
  
  return(consistent)
}


# Appends column of FC purity to a count table.
# Treating GoF and LoF experiments as two classes, gives a measure [0,1] for
# each gene of whether most experiments had the same up/down direction within
# each class. Doesn't consider consistency (eg, purity can be 1 when all GoF
# experiments are downreg and all LoF experiments are also downreg)
# https://stats.stackexchange.com/questions/95731/how-to-calculate-purity
# https://nlp.stanford.edu/IR-book/html/htmledition/evaluation-of-clustering-1.html
# TODO: this doesn't address ties, which are just auto-assigned to down-reg...


add_purity <- function(count_table, use_total = "all") {
  
  stopifnot(use_total %in% c("all", "gene"))
  
  stopifnot(
    c("GoF_down", "GoF_up", "LoF_down", "LoF_up", "Count_NA") %in% colnames(count_table))
  
  # are genes primarily up or down regulated in each of GoF/LoF experiments
  gof_class <- ifelse(count_table$GoF_down > count_table$GoF_up, "GoF_down", "GoF_up")
  lof_class <- ifelse(count_table$LoF_down > count_table$LoF_up, "LoF_down", "LoF_up")
  
  if (use_total == "all") {  # infer total fixed # of experiments
    
    all_msrd <- which(count_table$Count_NA == 0)[1]
    total <- sum(count_table[all_msrd, c("LoF_up", "LoF_down", "GoF_up", "GoF_down")])
    
    count_table$FC_purity <- vapply(1:nrow(count_table), function(x) {
      (count_table[x, gof_class[x]] +  count_table[x, lof_class[x]]) / total
    }, FUN.VALUE = numeric(1))
    
  } else {  # each gene has its own total of non-NAs
    
    total <- rowSums(count_table[, c("LoF_up", "LoF_down", "GoF_up", "GoF_down")])
    
    count_table$FC_purity <- vapply(1:nrow(count_table), function(x) {
      (count_table[x, gof_class[x]] +  count_table[x, lof_class[x]]) / total[x]
    }, FUN.VALUE = numeric(1))
    
  }
  return(count_table)
}



# If Consistent_across is FALSE (eg, dominant FC direction is same for LoF
# and GoF), make Purity negative

add_signed_purity <- function(count_table) {
 
  stopifnot("FC_purity" %in% colnames(count_table))
  
  consistent_across <- is_consistent_across(count_table)
  
  count_table <- count_table %>% 
    mutate(FC_signed_purity = ifelse(consistent_across, FC_purity, -FC_purity))
  
  return(count_table)
}


# Get the absolute average fold change for each gene in fc_mat as a df

avg_abs_fc <- function(fc_mat) {
  
  data.frame(Avg_abs_FC = rowMeans(abs(fc_mat), na.rm = TRUE)) %>% 
    rownames_to_column(var = "Symbol")
}



# Return a dataframe that contains for each symbol:
# Count_DE: Tally of gene being diff expr at FDR cutoff across experiments
# Proportion_DE_measured: Prop. of DE across experiments with non-NAs for gene
# Proportion_DE_all: Prop. of DE across all experiments
# Count_NA: Tally of NA measurements for the gene across experiments
# GoF_down: Tally when a gene had FC < 0 in overexpression experiments
# GoF_up: Tally when a gene had FC > 0 in overexpression experiments
# LoF_down: Tally when a gene had FC < 0 in KO and KD experiments
# LoF_up: Tally when a gene had FC > 0 in KO and KD experiments
# Avg_abs_FC: A gene's average absolute fold change across experiments
# FC_purity: Measure of consistency of FC direction within GoF/LoF experiments
# FC_signed_purity: FC_purity but make negative if GoF and LoF agree in direction
# DE_prior_rank: Ranks genes by how likely they are to be DE in diverse studies

process_all <- function(fdr_mat, 
                        fdr = 0.1, 
                        fc_mat, 
                        meta, 
                        purity_arg = "gene",
                        de_prior = NULL) {
  
  
  de_count <- count_table(fdr_mat, fdr)
  
  direction <- tally_direction(fc_mat, meta)
  
  fc <- avg_abs_fc(fc_mat)
  
  df <- left_join(de_count, direction, by = "Symbol") %>% 
    left_join(fc, by = "Symbol") %>% 
    add_purity(use_total = purity_arg) %>% 
    add_signed_purity()
  
  if (!is.null(de_prior)) {
    df <- left_join(df, de_prior[, c("Symbol", "DE_prior_rank")], by = "Symbol")
  }
  
  return(df)
}



# Performs process_all() for each unique symbol's experiments in meta/FC/FDR
# and returns a list of a df for each symbol

process_tf <- function(fdr_mat, 
                       fdr = 1, 
                       fc_mat, 
                       meta, 
                       use_total = "all",
                       de_prior = NULL,
                       ortho = FALSE,
                       ncores = 8) {
  
  if (ortho) meta$Symbol <- str_to_upper(meta$Symbol)
  
  tfs <- unique(meta$Symbol)
  
  tf_list <- mclapply(tfs, function(x) {
    
    tf_meta <- filter(meta, Symbol == x)
    
    process_all(fdr_mat = fdr_mat[, tf_meta$Experiment_ID], 
                fdr = fdr, 
                fc_mat = fc_mat[, tf_meta$Experiment_ID], 
                meta = tf_meta,
                de_prior = de_prior)
  }, mc.cores = ncores)
  names(tf_list) <- tfs
  
  return(tf_list)
}



# Given a list of DE counts for ortho and each species, return the ortho list
# with two extra columns for the counts deriving from each species

merge_ortho_counts <- function(count_ortho, count_hg, count_mm, pc_ortho) {
  
  tfs <- str_to_upper(c(names(count_ortho), names(count_hg), names(count_mm)))
  tfs <- unique(tfs)
  
  tf_list <- lapply(tfs, function(x) {
    
    mm <- count_mm[[str_to_title(x)]][, c("Symbol", "Count_DE")] %>% 
      filter(Symbol %in% pc_ortho$Symbol_mm) %>% 
      dplyr::rename(Symbol_mm = Symbol, Count_DE_mm = Count_DE) %>% 
      left_join(pc_ortho, by = "Symbol_mm") %>% 
      dplyr::rename(Symbol = ID)
    
    hg <- count_hg[[str_to_upper(x)]][, c("Symbol", "Count_DE")] %>% 
      filter(Symbol %in% pc_ortho$Symbol_hg) %>% 
      dplyr::rename(Symbol_hg = Symbol, Count_DE_hg = Count_DE) %>% 
      left_join(pc_ortho, by = "Symbol_hg") %>% 
      dplyr::rename(Symbol = ID)
    
    left_join(count_ortho[[x]],
              mm[, c("Count_DE_mm", "Symbol")],
              by = "Symbol"
              ) %>%
    left_join(hg[, c("Count_DE_hg", "Symbol")],
              by = "Symbol")
    
  })
  names(tf_list) <- tfs
  
  return(tf_list)
}
