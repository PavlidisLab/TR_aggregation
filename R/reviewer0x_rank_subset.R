## This script isolates experiments of a specified cell-type for a given TR, 
## calculates similarity within and across these experiments, and generates
## aggregated rankings within these experiments.
##------------------------------------------------------------------------------

library(tidyverse)
library(ggrepel)
source("R/setup-01_config.R")
source("R/utils/ranking_functions.R")
source("R/utils/perturbmatrix_functions.R")
source("R/utils/similarity_functions.R")
source("R/utils/plot_functions.R")

# List of all TR rankings and data matrices
rank_list <- readRDS(rank_path)
dat_list <- readRDS(alldat_path)

# List of all ChIP-seq exp similarity for inspection
bind_exp_sim <- readRDS(chip_sim_path)
perturb_exp_sim <- readRDS(perturb_sim_path)


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

mat_bind_tr <- dat_list$Binding$Human$QN_log[, exp_bind_tr, drop = FALSE]
mat_bind_ct <- dat_list$Binding$Human$QN_log[, exp_bind_ct, drop = FALSE]
mat_bind_nonct <- dat_list$Binding$Human$QN_log[, exp_bind_nonct, drop = FALSE]

mat_fc_tr <- dat_list$Perturbation$Human$FC_mat[, exp_perturb_tr, drop = FALSE]
mat_fc_ct <- dat_list$Perturbation$Human$FC_mat[, exp_perturb_ct, drop = FALSE]
mat_fc_nonct <- dat_list$Perturbation$Human$FC_mat[, exp_perturb_nonct, drop = FALSE]

mat_fdr_tr <- dat_list$Perturbation$Human$FDR_mat[, exp_perturb_tr, drop = FALSE]
mat_fdr_ct <- dat_list$Perturbation$Human$FDR_mat[, exp_perturb_ct, drop = FALSE]
mat_fdr_nonct <- dat_list$Perturbation$Human$FDR_mat[, exp_perturb_nonct, drop = FALSE]


# Similarity
#-------------------------------------------------------------------------------


# Re-purposing sim_list->df function to make groups reflective of cell type.
# Mixed means cell type of interest paired with a different cell type but same TR

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


# ChIP-seq experiment similarity

bind_sim <- chip_sim_list(mat_bind_tr, topn = 500)
bind_df <- format_and_merge2(bind_sim, exp_bind_ct)
# view(tf_summary(bind_df))

bind_top_mix <- bind_df %>%
  filter(Group == "Mixed") %>% 
  slice_max(Pcor, n = 10)

bind_top_nontr <- bind_exp_sim$df_list$Human %>% 
  filter((Row %in% exp_bind_ct | Col %in% exp_bind_ct) & Group == "Out") %>% 
  slice_max(Pcor, n = 10)

# Perturbation experiment similarity

perturb_sim <- perturb_sim_list(mat_fc_tr, mat_fdr_tr, topn = 500)
perturb_df <- format_and_merge2(perturb_sim, exp_perturb_ct)
# view(tf_summary(perturb_df))

perturb_topmix <- perturb_df %>%
  filter(Group == "Mixed") %>% 
  # slice_max(Pval_Intersect, n = 10)
  slice_max(Pcor_abs, n = 10)

perturb_top_nontr <- perturb_exp_sim$df_list$Human %>% 
  filter((Row %in% exp_perturb_ct | Col %in% exp_perturb_ct) & Group == "Out") %>% 
  slice_max(Pval_Intersect, n = 10)
  # slice_max(Pcor_abs, n = 10)


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
  left_join(rank_list$Human[[tr]], by = "Symbol") %>% 
  dplyr::rename(
    Rank_binding_all = Rank_binding,
    Rank_perturbation_all = Rank_perturbation,
    Rank_integrated_all = Rank_integrated
  ) %>% 
  arrange(Rank_integrated_in)
            


# Cor of ranks (removed all for downstream plotting)

cols <- colnames(all_rank)
cor_cols <- str_detect(cols, "^Rank.*") & !str_detect(cols, "_all$")
rank_cor <- cor(all_rank[, cor_cols])


# Genes most divergent in their rankings


diff_rank <- all_rank %>% 
  mutate(Diff_rank = Rank_integrated_in - Rank_integrated_out) %>%
  dplyr::select(Symbol, contains("Rank")) %>% 
  arrange(desc(Diff_rank))


# Looking at the count of genes mutually in the top k of ranks for in vs out
# cell type, as well as for the full rankings of the rest of the TRs

k <- 1000

n_in_out <- mclapply(1:k, function(x) {
  rank_in <- all_rank$Symbol[1:x]
  rank_out <- arrange(all_rank, Rank_integrated_out)$Symbol[1:x]
  data.frame(K = x, Common = length(intersect(rank_in, rank_out)))
}, mc.cores = cores)

n_in_out <- do.call(rbind, n_in_out)


# n_in_trs <- lapply(rank_list$Human, function(x) {
#   length(intersect(all_rank$Symbol[1:k], arrange(x, Rank_integrated)$Symbol[1:k]))
# })


n_in_trs <- mclapply(rank_list$Human, function(df) {
  
  n <- mclapply(1:k, function(x) {
    rank_in <- all_rank$Symbol[1:x]
    rank_out <- arrange(df, Rank_integrated)$Symbol[1:x]
    data.frame(K = x, Common = length(intersect(rank_in, rank_out)))
  })
  
  do.call(rbind, n)

}, mc.cores = cores)



# Performance of rankings


pr_list <- list(
  In_cell_type = arrange(all_rank, Rank_integrated_in) %>% get_pr(),
  Out_cell_type = arrange(all_rank, Rank_integrated_out) %>% get_pr(),
  All = arrange(all_rank, Rank_integrated_all) %>%  get_pr()
)


auprc_list <- list(
  In_cell_type = arrange(all_rank, Rank_integrated_in) %>% get_auprc(),
  Out_cell_type = arrange(all_rank, Rank_integrated_out) %>% get_auprc(),
  All = arrange(all_rank, Rank_integrated_all) %>%  get_auprc()
)


# Isolating specific gene effect sizes


get_gene_stats <- function(gene,
                           exp_bind_tr, 
                           exp_bind_ct, 
                           mat_bind_tr, 
                           exp_perturb_tr,
                           exp_perturb_ct, 
                           mat_fc_tr, 
                           mat_fdr, 
                           fdr = 0.1) {
  list(
    Binding = data.frame(
      Experiment = exp_bind_tr,
      Binding = mat_bind_tr[gene, ],
      CT = exp_bind_tr %in% exp_bind_ct
    ),
    Perturbation = data.frame(
      Experiment = exp_perturb_tr,
      FC = mat_fc_tr[gene, ],
      Abs_FC = abs(mat_fc_tr[gene, ]),
      DE = mat_fdr_tr[gene, ] < fdr,
      CT = exp_perturb_tr %in% exp_perturb_ct
    )
  )
}


# LGALS12  NTRK1 TGFBR2 GBX2  TBX19  EXOC3L2  ALX4  CR1

gene <- "TGFBR2"  

gene_stats <- get_gene_stats(
  gene,
  exp_bind_tr,
  exp_bind_ct,
  mat_bind_tr,
  exp_perturb_tr,
  exp_perturb_ct,
  mat_fc_tr,
  mat_fdr
)


# Plots
#-------------------------------------------------------------------------------


# Violin boxplot of similarity statistics between experiments in/out/mixed cell types

stat_vboxplot <- function(df, 
                          x_var = "Group", 
                          y_var, 
                          y_name, 
                          title,
                          ortho = FALSE) {
  
  df$Group <- str_replace_all(df$Group, "_", " ")
  
  ggplot(df, aes(x = !!sym(x_var), y = !!sym(y_var))) +
    geom_violin(width = 0.4, fill = "lightslategrey") +
    geom_boxplot(width = 0.1, fill = "white") +
    ylab(y_name) +
    ggtitle(title) +
    ylab(y_name) +
    theme_classic() +
    theme(axis.title = element_text(size = 30),
          axis.title.x = element_blank(),
          axis.text.y = element_text(size = 25),
          axis.text.x = element_text(size = 25),
          plot.title = element_text(size = 25))
}


p1a <- stat_vboxplot(bind_df, y_var = "Pcor", y_name = "Pearson's correlation", title = paste(tr, ct, "binding"))
p1b <- stat_vboxplot(bind_df, y_var = "Intersect", y_name = "Top 500 overlap", title = paste(tr, ct, "binding"))
p1c <- stat_vboxplot(perturb_df, y_var = "Pcor_abs", y_name = "Pearson's correlation", title = paste(tr, ct, "perturbation"))
p1d <- stat_vboxplot(perturb_df, y_var = "Pval_Intersect", y_name = "Top 500 overlap", title = paste(tr, ct, "perturbation"))
p1 <- plot_grid(p1a, p1b, p1c, p1d, nrow = 2, ncol = 2)

ggsave(p1, height = 12, width = 18, dpi = 300, device = "png",
       filename = file.path(iplot_dir, "Experiment_similarity", "runx1_k562_similarity.png"))


# Boxplot of effect gene effect sizes in/out cell type


gene_boxplot <- function(df, yvar, yname, title, de = FALSE) {
  
  df <- filter(df, !is.na(!!sym(yvar)))
  
  p <- 
    ggplot(df, aes(x = CT, y = !!sym(yvar))) +
    geom_boxplot(width = 0.2, fill = "slategrey", outlier.shape = NA)
  
  if (de) {
    p <- p + 
      geom_jitter(aes(colour = DE), size = 3, shape = 21, width = 0.2) +
      scale_colour_manual(values = c("black", "red"))
    
  } else {
      
    p <- p + 
      geom_jitter(colour = "black",
                  size = 3, shape = 21,  width = 0.3, alpha = 0.5)
  }
  
  p <- p +
    # geom_label_repel(data = filter(df, Group),
                    # aes(x = CT, y = !!sym(yvar), label = Cell_Type)) +
    ylab(yname) +
    ggtitle(title) +
    theme_classic() +
    theme(axis.title = element_text(size = 30),
          axis.title.x = element_blank(),
          axis.text.y = element_text(size = 25),
          axis.text.x = element_text(size = 20),
          plot.title = element_text(size = 30, face = "italic"),
          legend.position = "top",
          legend.text = element_text(size = 20),
          legend.title = element_text(size = 20),
          legend.background = element_blank(),
          legend.box.background = element_rect(colour = "black"))
  
  return(p)
}


p2a <- gene_boxplot(gene_stats$Binding, yvar = "Binding", yname = "Binding score", title = gene)
p2b <- gene_boxplot(gene_stats$Perturbation, yvar = "Abs_FC", yname = "Absolute FC", title = gene, de = TRUE) + theme(plot.title = element_blank())
# p2b <- gene_boxplot(gene_stats$Perturbation, yvar = "FC", yname = "FC", title = gene, de = TRUE)
p2 <- plot_grid(p2a, p2b, nrow = 1)
p2 <- ggdraw(add_sub(p2, ct, vpadding = grid::unit(0, "lines"), y = 6, x = 0.5, vjust = 5.5, size = 30))


ggsave(p2, height = 6, width = 10, dpi = 300, device = "png", bg = "white",
       filename = file.path(iplot_dir, "Gene_rankings", paste(tr, ct, gene, "genestats.png", sep = "_")))


# Precision vs recall of in/out/all rankings

auprc_labels <- c(
  In_cell_type = paste0("In cell type (AUC=", signif(auprc_list$In_cell_type, 3), ")"),
  Out_cell_type = paste0("Out cell type (AUC=", signif(auprc_list$Out_cell_type, 3), ")"),
  All = paste0("All (AUC=", signif(auprc_list$All, 3), ")"))


colours <- c(
  "In_cell_type" = "#e41a1c",
  "Out_cell_type" = "#4daf4a",
  "All" = "#377eb8"
)


p3 <- ggplot() +
  geom_path(data = pr_list$In_cell_type, aes(x = Recall, y = Precision, col = "In_cell_type"), size = 1) +
  geom_path(data = pr_list$Out_cell_type, aes(x = Recall, y = Precision, col = "Out_cell_type"), size = 1) +
  geom_path(data = pr_list$All, aes(x = Recall, y = Precision, col = "All"), size = 1) +
  ggtitle(paste(tr, ct)) +
  scale_color_manual(labels = auprc_labels, values = colours) +
  theme_classic() +
  theme(axis.text = element_text(size = 25),
        axis.title = element_text(size = 30),
        plot.title = element_text(size = 30),
        legend.title = element_blank(),
        legend.text = element_text(size = 25),
        legend.position = c(0.55, 0.85))


ggsave(p3, height = 8, width = 8, dpi = 300, device = "png",
       filename = file.path(iplot_dir, "Gene_rankings", "runx1_k562_prcurve.png"))



# Line plot showing common genes for given k cut-off between rankings


plot_df <- data.frame(
  K = 1:k,
  TR = c(
    rep("RUNX1", 2 * k),
    rep("ASCL1", k),
    rep("HES1", k),
    rep("MECP2", k),
    rep("MEF2C", k),
    rep("NEUROD1", k),
    rep("PAX6", k),
    rep("TCF4", k)
  ),
  Non_K562 = c(rep(FALSE, k), rep(TRUE, k), rep(FALSE, k * 7)),
  Common = c(
    n_in_trs$RUNX1$Common,
    n_in_out$Common,
    n_in_trs$ASCL1$Common,
    n_in_trs$HES1$Common,
    n_in_trs$MECP2$Common,
    n_in_trs$MEF2C$Common,
    n_in_trs$NEUROD1$Common,
    n_in_trs$PAX6$Common,
    n_in_trs$TCF4$Common
  )
)



p4 <- ggplot(plot_df, aes(x = K, y = Common, colour = TR, linetype = Non_K562)) +
  geom_line(linewidth = 1.4) +
  ggtitle(paste(tr, ct)) +
  xlab("List cutoff") +
  ylab("Common genes") +
  scale_color_manual(values = tf_pal_hg) +
  theme_classic() +
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 25),
        plot.title = element_text(size = 25),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20),
        plot.margin = margin(10, 15, 10, 10))



ggsave(p4, height = 8, width = 10, dpi = 300, device = "png",
       filename = file.path(iplot_dir, "Experiment_similarity", "runx1_k562_topn_common.png"))



# Heatmap of ranking correlation

rank_cor_subset <- round(rank_cor[str_detect(rownames(rank_cor), "_out$"), str_detect(colnames(rank_cor), "_in$")], 3)
colnames(rank_cor_subset) <- str_to_title(str_replace(colnames(rank_cor_subset), "^Rank_", ""))
rownames(rank_cor_subset) <- str_to_title(str_replace(rownames(rank_cor_subset), "^Rank_", ""))

pheatmap(rank_cor_subset, 
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         color = bluered_pal,
         display_numbers = TRUE,
         number_color = "black",
         fontsize = 20,
         cellwidth = 40,
         cellheight = 40,
         filename = file.path(iplot_dir, "Gene_rankings", paste(tr, ct, "rank_cor.png", sep = "_")))

# Heatmaps of experiment similarity

pheatmap(bind_sim$Intersect[c(exp_bind_ct, exp_bind_nonct), c(exp_bind_ct, exp_bind_nonct)],
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         color = bluered_pal,
         gaps_row = length(exp_bind_ct),
         gaps_col = length(exp_bind_ct))

pheatmap(bind_sim$Intersect,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         color = bluered_pal)

pheatmap(perturb_sim$Pcor_abs,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         color = bluered_pal)
