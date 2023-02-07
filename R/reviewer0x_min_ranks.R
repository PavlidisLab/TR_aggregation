## This script is for implementing an alternative ranking method, where genes
## are assigned their lowest rank across TR-specific experiments. These ranks
## are then compared to the original aggregate rankings.
## -----------------------------------------------------------------------------

library(tidyverse)
library(data.table)
source("R/setup-01_config.R")
source("R/utils/ranking_functions.R")
source("R/utils/perturbmatrix_functions.R")
source("R/utils/similarity_functions.R")
source("R/utils/plot_functions.R")

# List of all TR rankings and data matrices
rank_list <- readRDS(paste0(scratch_dir, date, "_ranked_target_list.RDS"))
dat_list <- readRDS(paste0(scratch_dir, date, "_all_data_list.RDS"))

tr_hg <- unique(filter(dat_list$Binding$Meta, Species == "Human")$Symbol)
tr_mm <- unique(filter(dat_list$Binding$Meta, Species == "Mouse")$Symbol)


# Convert matrix of effect sizes to ranks WITHIN experiments
# ------------------------------------------------------------------------------


# Perturbation can have ties due to pvals being coerced to 0. Break ties with FC

perturb_rank <- function(pval_mat, fc_mat) {
  
  rank_l <- lapply(1:ncol(fc_mat), function(x) {
    data.table::frank(list(pval_mat[, x], -fc_mat[, x]), ties.method = "min")
  })
  
  rank_mat <- do.call(cbind, rank_l)
  colnames(rank_mat) <- colnames(fc_mat)
  rownames(rank_mat) <- rownames(fc_mat)
  return(rank_mat)
}


# Binding
bind_hg <- dat_list$Binding$Human$QN_log
bind_mm <- dat_list$Binding$Mouse$QN_log
bind_rank_hg <- apply(-bind_hg, 2, rank, ties = "min")
bind_rank_mm <- apply(-bind_mm, 2, rank, ties = "min")


# Perturbation
fc_hg <- abs(dat_list$Perturbation$Human$FC_mat)
fc_mm <- abs(dat_list$Perturbation$Mouse$FC_mat)
pval_hg <- dat_list$Perturbation$Human$Pval_mat
pval_mm <- dat_list$Perturbation$Mouse$Pval_mat
perturb_rank_hg <- perturb_rank(pval_hg, fc_hg)
perturb_rank_mm <- perturb_rank(pval_mm, fc_mm)


# For each TR, get each gene's minimum rank ACROSS perturbation/binding 
# experiments. Break ties using the average FC/binding across experiments.
# ------------------------------------------------------------------------------


# Return a list of dfs, one for each TR, which has each gene's minimum rank 


get_bind_min_ranks <- function(bind_mat, rank_mat, meta, trs) {
  
  rank_l <- lapply(trs, function(x) {
    
    exps <- filter(meta, Symbol == x)$Experiment_ID
    rank_min <- apply(rank_mat[, exps, drop = FALSE], 1, min)
    which_min <- exps[apply(rank_mat[, exps, drop = FALSE], 1, which.min)]
    avg_bind <- rowMeans(bind_mat[, exps, drop = FALSE])
    rank_min_avg <- data.table::frank(list(rank_min, -avg_bind), ties.method = "min")
    
    rank_df <- data.frame(Symbol = rownames(rank_mat),
                          Rank_binding_min = rank_min,
                          Binding_min_exp = which_min,
                          Rank_binding_min_avg = rank_min_avg)
    
  })
  names(rank_l) <- trs
  return(rank_l)
}



get_perturb_min_ranks <- function(fc_mat, rank_mat, meta, trs) {
  
  rank_l <- lapply(trs, function(x) {
    
    exps <- filter(meta, Symbol == x)$Experiment_ID
    rank_min <- apply(rank_mat[, exps, drop = FALSE], 1, min, na.rm = TRUE)
    which_min <- exps[apply(rank_mat[, exps, drop = FALSE], 1, which.min)]
    avg_fc <- rowMeans(abs(fc_mat[, exps, drop = FALSE]), na.rm = TRUE)
    rank_min_avg <- data.table::frank(list(rank_min, -avg_fc), ties.method = "min")
    
    rank_df <- data.frame(Symbol = rownames(rank_mat),
                          Rank_perturbation_min = rank_min,
                          Perturb_min_exp = which_min,
                          Rank_perturbation_min_avg = rank_min_avg)
    
  })
  names(rank_l) <- trs
  return(rank_l)
  
}


# Get the list of ranks for each data type, join, then add the integrated 
# ranking (rank products) of the minimum rankings using average tie breaks

get_all_min_ranks <- function(bind_mat, 
                              bind_rank_mat, 
                              bind_meta, 
                              fc_mat, 
                              perturb_rank_mat, 
                              perturb_meta, 
                              trs) {
  
  bind_l <- get_bind_min_ranks(bind_mat, bind_rank_mat, bind_meta, trs)
  perturb_l <- get_perturb_min_ranks(fc_mat, perturb_rank_mat, perturb_meta, trs)
  
  # Join and add rank product
  
  rank_l <- lapply(trs, function(x) {
    
    df <- left_join(bind_l[[x]], perturb_l[[x]], by = "Symbol") %>% 
      mutate(
        Rank_integrated_min = rank(
          (Rank_binding_min_avg/nrow(bind_mat)) * (Rank_perturbation_min_avg/nrow(fc_mat)),
          ties.method = "min")
      )
  })
  names(rank_l) <- trs
  return(rank_l)
}




rank_hg <- get_all_min_ranks(bind_mat = bind_hg,
                             bind_rank_mat = bind_rank_hg,
                             bind_meta = dat_list$Binding$Meta,
                             fc_mat = fc_hg,
                             perturb_rank_mat = perturb_rank_hg,
                             perturb_meta = dat_list$Perturbation$Meta,
                             trs = tr_hg)


rank_mm <- get_all_min_ranks(bind_mat = bind_mm,
                             bind_rank_mat = bind_rank_mm,
                             bind_meta = dat_list$Binding$Meta,
                             fc_mat = fc_mm,
                             perturb_rank_mat = perturb_rank_mm,
                             perturb_meta = dat_list$Perturbation$Meta,
                             trs = tr_mm)


# Join minimum ranks with aggregate ranking for comparisons
# ------------------------------------------------------------------------------


join_ranks <- function(agg_list, min_rank_list, trs) {
  
  rank_l <- lapply(trs, function(x) {
    
    left_join(agg_list[[x]], min_rank_list[[x]], by = "Symbol") %>% 
      mutate(
        Diff_binding = Rank_binding - Rank_binding_min_avg,
        Diff_perturbation = Rank_perturbation - Rank_perturbation_min_avg,
        Diff_integrated = Rank_integrated - Rank_integrated_min
      ) %>% 
      arrange(Rank_integrated)
  })
  
  names(rank_l) <- trs
  return(rank_l)
}


all_rank_hg <- join_ranks(rank_list$Human, rank_hg, tr_hg)
all_rank_mm <- join_ranks(rank_list$Mouse, rank_mm, tr_mm)


# Demonstrating step-like nature of min ranks without breaking ties with
# average bind/FC

a1 <- arrange(all_rank_mm$Ascl1, Rank_binding_min)
plot(a1$Rank_binding_min[1:200], cex = 1.4, pch = 19)
plot(a1$Rank_binding_min_avg[1:200], cex = 1.4, pch = 19)


# Get the correlation of the min and aggregated ranking


get_cor <- function(rank_list) {
  
  keep_cols <- c("Rank_binding", "Rank_perturbation", "Rank_integrated")
  keep_rows <- paste0(keep_cols, "_min")
  
  cor_l <- lapply(rank_list, function(x) {
    
    cor_mat <- cor(x[, c(keep_cols, keep_rows)])
    
    data.frame(
      Integrated = cor_mat["Rank_integrated", "Rank_integrated_min"],
      Binding = cor_mat["Rank_binding", "Rank_binding_min"],
      Perturbation = cor_mat["Rank_perturbation", "Rank_perturbation_min"])
  })
  
  cor_df <- data.frame(do.call(rbind, cor_l))
  rownames(cor_df) <- names(rank_list)
  return(cor_df)
}


cor_hg <- get_cor(all_rank_hg)
cor_mm <- get_cor(all_rank_mm)


# Heatmap of cor

pal <- c('#ffffcc','#ffeda0','#fed976','#feb24c','#fd8d3c','#fc4e2a','#e31a1c','#bd0026','#800026')


plot_heatmap <- function(cor_mat, pal) {
  
  pheatmap(cor_mat,
           cluster_rows = FALSE,
           cluster_cols = FALSE,
           color = pal,
           display_numbers = TRUE,
           fontsize = 30,
           number_color = "black")
}


plot_heatmap(cor_hg, pal)
plot_heatmap(cor_mm, pal)


# Using human ASCL1 to pull examples of genes that have the biggest change in 
# integrated rankings

exps_perturb <-  filter(dat_list$Perturbation$Meta, Symbol == "ASCL1")$Experiment_ID


# PRSS35 sees biggest drop from aggregated to minimum rankings. Driven by a drop
# in the perturbation ranking... ASCL1 has a lot of DE genes. So while PRSS35
# is freq DE across ASCL1 experiments (+ good abs FC), it is "unimpressive"
# within experiments.

maxneg_change <- slice_min(all_rank_hg$ASCL1, Diff_integrated, n = 1)
perturb_rank_hg[maxneg_change$Symbol, exps_perturb]
pval_hg[maxneg_change$Symbol, exps_perturb]
mean(fc_hg[maxneg_change$Symbol, exps_perturb], na.rm = TRUE)
sum(dat_list$Perturbation$Human$FDR_mat[maxneg_change$Symbol, exps_perturb] < fdr, na.rm = TRUE)


# SLC9C2 has biggest rise from aggregated to minimum rankings. Also driven by
# the perturbation data: has 5/8 NAs, but is ranked 21st in one of the perturb
# experiments in which it is measured. However, it is not actually DE in that
# experiment, which has a wonky pval distn. Pval=0.00072, FDR=0.52714

maxpos_change <- slice_max(all_rank_hg$ASCL1, Diff_integrated, n = 1)
perturb_rank_hg[maxpos_change$Symbol, exps_perturb]
pval_hg[maxpos_change$Symbol, exps_perturb]
mean(fc_hg[maxpos_change$Symbol, exps_perturb], na.rm = TRUE)
dat_list$Perturbation$Human$FDR_mat[maxpos_change$Symbol, exps_perturb]
hist(pval_hg[, "GSE151000_ASCL1_Human_Knockdown"])


# Plot of difference in ranks

plot_df <- all_rank_hg$ASCL1 %>%
  mutate(Group = Symbol %in% c(maxpos_change$Symbol, maxneg_change$Symbol))

ggplot() +
  geom_point(aes(x = Rank_integrated, y = Diff_integrated), 
             data = filter(plot_df, !Group),
             shape = 21, alpha = 0.2) +
  geom_point(aes(x = Rank_integrated, y = Diff_integrated), 
             data = filter(plot_df, Group),
             shape = 21, size = 3, fill = "red", colour = "black") +
  ylab("Change in integrated ranking") +
  xlab("Original integrated ranking") +
  theme_classic() +
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 20))


# Compare AUPRC of aggregate vs min rank
# ------------------------------------------------------------------------------


# Repurpose existing function to retrieve AUPRC of perturb/binding/integrated

get_all_auprc <- function(rank_df) {
  
  # Order by respective rankings
  df_int <- arrange(rank_df, Rank_integrated)
  df_int_min <- arrange(rank_df, Rank_integrated_min)
  df_bind <- arrange(rank_df, Rank_binding)
  df_bind_min <- arrange(rank_df, Rank_binding_min)
  df_perturb <- arrange(rank_df, Rank_perturbation)
  df_perturb_min <- arrange(rank_df, Rank_perturbation_min)
  
  
  data.frame(
    Integrated = get_auprc(df_int),
    Integrated_min = get_auprc(df_int_min),
    Binding = get_auprc(df_bind),
    Binding_min = get_auprc(df_bind_min),
    Perturbation = get_auprc(df_perturb),
    Perturbation_min = get_auprc(df_perturb_min)
  )
}


auprc_hg <- lapply(all_rank_hg, get_all_auprc)
auprc_mm <- lapply(all_rank_mm, get_all_auprc)


# Find that the aggregate ranking typically exceeds the min rank approach.

gt_hg <- unlist(lapply(auprc_hg, function(x) x["Integrated"] > x["Integrated_min"]))
gt_mm <- unlist(lapply(auprc_mm, function(x) x["Integrated"] > x["Integrated_min"]))


# Bar plots of the AUPRC across rankings

auprc_barplot <- function(auprc_list) {
  
  tfs <- names(auprc_list)
  
  plot_l <- lapply(tfs, function(x) {
    
    auprc_vec <- round(auprc_list[[x]], 4)
    
    df <- data.frame(
      Rank = factor(names(auprc_vec), levels = unique(names(auprc_vec))),
      AUPRC = unlist(auprc_vec)
    )
    
    ggplot(df, aes(x = Rank, y = AUPRC)) +
      geom_bar(stat = "identity") +
      ggtitle(x) +
      theme_classic() +
      theme(axis.text = element_text(size = 20),
            axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
            axis.title.x = element_blank(),
            axis.title.y = element_text(size = 20))
  })
  
  names(plot_l) <- tfs
  return(plot_l)
}



bar_hg <- auprc_barplot(auprc_hg)
bar_mm <- auprc_barplot(auprc_mm)


# Using mouse Runx1 as example to inspect where min ranking yielded greater
# AUPRC than aggregate. Check Prec/recall/n recovered at different steps k.
# Find that min ranking has an extra curated target at an early step (k=10), but
# at higher ks the aggregated ranking recovers more targets.

rank_df <- all_rank_mm$Runx1
k <- c(10, 100, 500, nrow(rank_df))


auprc_k_list <- lapply(k, function(x) {
  
  pr_agg <- rank_df %>% arrange(Rank_integrated) %>% get_pr()
  pr_min <- rank_df %>% arrange(Rank_integrated_min) %>% get_pr()
  
  n_agg <- rank_df %>%
    slice_min(Rank_integrated, n = x) %>%
    pull(Curated_target) %>%
    sum()
  
  n_min <- rank_df %>% 
    slice_min(Rank_integrated_min, n = x) %>% 
    pull(Curated_target) %>% 
    sum()
  
  data.frame(
    Precison_agg = pr_agg$Precision[x],
    Recall_agg = pr_agg$Recall[x],
    nTarget_agg = n_agg,
    Precison_min = pr_min$Precision[x],
    Recall_min = pr_min$Recall[x],
    nTarget_min = n_min
  )
})


auprc_k_df <- round(do.call(rbind, auprc_k_list), 4)
rownames(auprc_k_df) <- paste0("step=", k)
