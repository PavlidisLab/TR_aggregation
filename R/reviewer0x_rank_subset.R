## This script isolates experiments of a specified cell-type for a given TR, 
## calculates similarity within and across these experiments, and generates
## aggregated rankings within these experiments.
##------------------------------------------------------------------------------

library(tidyverse)
source("R/setup-01_config.R")
source("R/utils/ranking_functions.R")
source("R/utils/perturbmatrix_functions.R")
source("R/utils/similarity_functions.R")
source("R/utils/plot_functions.R")

# List of all TR rankings and data matrices
rank_list <- readRDS(paste0(scratch_dir, date, "_ranked_target_list.RDS"))
dat_list <- readRDS(paste0(scratch_dir, date, "_all_data_list.RDS"))


# Isolate experiments for a given TR and given context
#-------------------------------------------------------------------------------


ct <- "K-562"  # "Kasumi-1"  "K-562"
tr <- "RUNX1"

exp_bind_tr <- filter(dat_list$Binding$Meta, Symbol == tr)$Experiment_ID
exp_bind_ct <- filter(dat_list$Binding$Meta, Symbol == tr & Cell_Type == ct)$Experiment_ID
exp_bind_nonct <- filter(dat_list$Binding$Meta, Symbol == tr & Cell_Type != ct)$Experiment_ID

exp_perturb_tr <- filter(dat_list$Perturbation$Meta, Symbol == tr)$Experiment_ID
exp_perturb_ct <- filter(dat_list$Perturbation$Meta, Symbol == tr & Cell_Type == ct)$Experiment_ID
exp_perturb_nonct <- filter(dat_list$Perturbation$Meta, Symbol == tr & Cell_Type != ct)$Experiment_ID

mat_bind_tr <- dat_list$Binding$Human$QN_log[, exp_bind_tr]
mat_bind_ct <- dat_list$Binding$Human$QN_log[, exp_bind_ct]
mat_bind_nonct <- dat_list$Binding$Human$QN_log[, exp_bind_nonct]

mat_fc_tr <- dat_list$Perturbation$Human$FC_mat[, exp_perturb_tr, drop = FALSE]
mat_fc_ct <- dat_list$Perturbation$Human$FC_mat[, exp_perturb_ct, drop = FALSE]
mat_fc_nonct <- dat_list$Perturbation$Human$FC_mat[, exp_perturb_nonct, drop = FALSE]

mat_fdr_tr <- dat_list$Perturbation$Human$FDR_mat[, exp_perturb_tr, drop = FALSE]
mat_fdr_ct <- dat_list$Perturbation$Human$FDR_mat[, exp_perturb_ct, drop = FALSE]
mat_fdr_nonct <- dat_list$Perturbation$Human$FDR_mat[, exp_perturb_nonct, drop = FALSE]


# Similarity
#-------------------------------------------------------------------------------


# Re-purposing sim_list->df function to make groups reflective of cell type.
# Mixed means cell type of interested paired with a different cell type.

format_and_merge2 <- function(sim_list, exp_vec) {
  
  df <- format_and_merge(sim_list)
  
  df <- df %>% 
    mutate(
      Group = case_when(
        Row %in% exp_vec & Col %in% exp_vec ~ "In_cell_type",
        !(Row %in% exp_vec) & !(Col %in% exp_vec) ~ "Out_cell_type",
        TRUE ~ "Mixed"
      ))
  
  df$Group <- factor(df$Group,
                     levels = c("In_cell_type", "Out_cell_type", "Mixed"))
  
  return(df)
  
}


bind_sim <- chip_sim_list(mat_bind_tr, topn = 500)
bind_df <- format_and_merge2(bind_sim, exp_bind_ct)


perturb_sim <- perturb_sim_list(mat_fc_tr, mat_fdr_tr, topn = 500)
perturb_df <- format_and_merge2(perturb_sim, exp_perturb_ct)


# Most similar 'Mixed' experiments: only one of experimental pair is in context
# of interest


bind_topmix <- bind_df %>%
  filter(Group == "Mixed") %>% 
  slice_max(Pcor, n = 10)

bind_topmix_ct <- dat_list$Binding$Meta %>%
  filter(Experiment_ID %in% setdiff(c(bind_topmix$Row, bind_topmix$Col), exp_bind_ct)) %>%
  pull(Cell_Type)


perturb_topmix <- perturb_df %>%
  filter(Group == "Mixed") %>% 
  slice_max(Pval_Intersect, n = 10)

perturb_topmix_ct <- dat_list$Perturbation$Meta %>%
  filter(Experiment_ID %in% setdiff(c(perturb_topmix$Row, perturb_topmix$Col), exp_perturb_ct)) %>%
  pull(Cell_Type)



# Rankings
#-------------------------------------------------------------------------------


get_rank_list <- function(mat_fdr, mat_fc, mat_bind, fdr = 0.1) {
  
  # Perturbation ranking
  count_de <- 
    left_join(count_table(mat_fdr, fdr = fdr),
              avg_abs_fc(mat_fc), 
              by = "Symbol") %>% 
    rank_perturb()
  
  # Binding ranking
  bind <- 
    data.frame(Mean_bind = rowMeans(mat_bind)) %>% 
    rownames_to_column(var = "Symbol") %>% 
    rank_binding()
  
  # Rank product
  rp <-
    left_join(count_de, bind, by = "Symbol") %>% 
    rank_product()
  
  return(rp)
}


rank_ct <- get_rank_list(mat_fdr = mat_fdr_ct,
                         mat_fc = mat_fc_ct,
                         mat_bind = mat_bind_ct)



rank_nonct <- get_rank_list(mat_fdr = mat_fdr_nonct,
                            mat_fc = mat_fc_nonct,
                            mat_bind = mat_bind_nonct)



# compare the rankings between full aggregate, context-specific, and experiments
# outside of the context


all_rank <- 
  left_join(rank_ct, rank_nonct, by = "Symbol", suffix = c("_in", "_out")) %>% 
  left_join(rank_list$Human[[tr]], by = "Symbol")
            

cor(all_rank[, c("Rank_integrated_in", "Rank_integrated_out", "Rank_integrated")])


# which genes diverge the most in their rankings

diff_rank <- all_rank %>% 
  mutate(Diff_rank = abs(Rank_integrated_in - Rank_integrated_out)) %>%
  dplyr::select(Symbol, contains("Rank")) %>% 
  arrange(desc(Diff_rank))


# GYPB most different - "major sialoglycoproteins of the human erythrocyte membrane"
# https://www.proteinatlas.org/ENSG00000250361-GYPB
# K562 noted as erythrocyte-like 
# GYPB assyaed in K562 https://pubmed.ncbi.nlm.nih.gov/7806496/

head(diff_rank, 10)



plot(diff_rank$Rank_integrated_in, diff_rank$Rank_integrated_out)


# Plots
#-------------------------------------------------------------------------------


stat_vboxplot <- function(df, 
                          x_var = "Group", 
                          y_var, 
                          y_name, 
                          title,
                          ortho = FALSE) {
  
  ggplot(df, aes(x = !!sym(x_var), y = !!sym(y_var))) +
    geom_violin(width = 0.4, fill = "lightslategrey") +
    geom_boxplot(width = 0.1, fill = "white") +
    ylab(y_name) +
    ggtitle(title) +
    ylab(y_name) +
    theme_classic() +
    theme(axis.title = element_text(size = 25),
          axis.title.x = element_blank(),
          axis.text.y = element_text(size = 25),
          axis.text.x = element_text(size = 25),
          plot.title = element_text(size = 25))
}


stat_vboxplot(bind_df, y_var = "Pcor", y_name = "Pcor", title = paste(tr, ct))

stat_vboxplot(perturb_df, y_var = "Pval_Intersect", y_name = "Pval_Intersect", title = paste(tr, ct))


pheatmap(bind_sim$Pcor[c(exp_bind_ct, exp_bind_nonct), c(exp_bind_ct, exp_bind_nonct)],
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         color = bluered_pal,
         gaps_row = length(exp_bind_ct),
         gaps_col = length(exp_bind_ct))


pheatmap(bind_sim$Pcor,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         color = bluered_pal)
