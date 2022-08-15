## This script looks at the overlap of the top scoring genes between ChIP-seq
## and perturbation experiments
## -----------------------------------------------------------------------------

library(tidyverse)
library(GeneOverlap)
library(pheatmap)
library(RColorBrewer)
library(viridis)
library(cowplot)
source("~/regnetR/R/utils/similarity_functions.R")
source("~/regnetR/R/utils/plot_functions.R")

topn <- 500  # number of top genes to consider when doing overlap
common_arg <- TRUE  # should comparison only be done for mutually measured genes?
date <- "Apr2022"  # most recent data freeze
outfile <- paste0("~/scratch/R_objects/", date, "_intersect_similarity.RDS")
plot_dir <- "~/Plots/Intersect/"

# Loading all data
dat <- readRDS(paste0("~/scratch/R_objects/", date, "_all_data_list.RDS"))


# General workflow is to generate exp x exp matrices where elements represent
# different metrics of similarity/overlap between the experiments. These are then
# turned into dfs of the unique pairs for plotting/summaries
# FOR EACH OF SORTING BY PVAL, UPREGULATED, AND DOWNREGULATED FOR PERTURBATION:
# 1) Length of the intersect of the topn overlap between experiments
# 2) Jaccard of the topn overlap
# 3) FET pval of the topn overlap
# 4) Log odds ratio of the topn overlap
#-------------------------------------------------------------------------------


# NOTE: Perturbation experiments with low coverage (lots of NAs) results in
# inflated overlap. Filter experiments with low coverage

min_msrd <- 10e3
count_msrd <- apply(dat$Perturbation$Ortho$Pval_mat, 2, function(x) sum(!is.na(x)))
keep_exp <- names(which(count_msrd > min_msrd))


fc_list <- lapply(dat$Perturbation, function(x) {
  mat <- x$FC_mat
  mat[, intersect(colnames(mat), keep_exp)]
})


pval_list <- lapply(dat$Perturbation, function(x) {
  mat <- x$Pval_mat
  mat[, intersect(colnames(mat), keep_exp)]
})



if (!file.exists(outfile)) {
  
  sim_list <- list(
    Human = intersect_sim_list(
      chip_mat = dat$Binding$Human$QN_log,
      fc_mat = fc_list$Human,
      pval_mat = pval_list$Human,
      topn = topn,
      common_arg = common_arg
    ),
    Mouse = intersect_sim_list(
      chip_mat = dat$Binding$Mouse$QN_log,
      fc_mat = fc_list$Mouse,
      pval_mat = pval_list$Mouse,
      topn = topn,
      common_arg = common_arg
    ),
    Ortho = intersect_sim_list(
      chip_mat = dat$Binding$Ortho$QN_log,
      fc_mat = fc_list$Ortho,
      pval_mat = pval_list$Ortho,
      topn = topn,
      common_arg = common_arg
    )
  )
  
  df_list <- lapply(sim_list, format_and_merge, symmetric_arg = FALSE)
  
  saveRDS(list(sim_list = sim_list, df_list = df_list), file = outfile)
  
} else {
  
  temp_load <- readRDS(outfile)
  df_list <- temp_load$df_list
  sim_list <- temp_load$sim_list
  rm(temp_load)
  
}



# Summarize/inspect group similarities (focus on what is in paper)
# GEO1/TF1/Row == ChIP-seq & GEO2/TF2/Col == Perturb
#-------------------------------------------------------------------------------


focus_cols <- c("Group", 
                "Pval_Intersect_Mean", "Pval_Intersect_Max",
                "Upreg_Intersect_Mean", "Upreg_Intersect_Max",
                "Downreg_Intersect_Mean", "Downreg_Intersect_Max")

all_summ <- lapply(df_list, get_summary)
tf_summ <- lapply(df_list, tf_summary)

# Group summaries by all experiments or TF-specific
lapply(all_summ, `[`, focus_cols)
lapply(tf_summ, function(tf) lapply(tf, `[`, focus_cols))

# top pairs by Pcor/Intersect
top_inter <- lapply(df_list, function(x) arrange(x, desc(Pval_Intersect)) %>% head(10))

# split by df TF
tf_df <- lapply(df_list, function(x) {
  tfs <- unique(c(as.character(x$TF1), as.character(x$TF2)))
  tf_list <- lapply(tfs, function(y) {
    filter(x, as.character(TF1) == y & as.character(TF2) == y)
  })
  names(tf_list) <- tfs
  return(tf_list)
})


# Function that returns tables tallying the count of times a gene
# was found to overlapping for the requested TF and summary/statlist
# ------------------------------------------------------------------------------


overlap_counts <- function(gene_list, pair_df, gene_mat) {
  
  tfs <- unique(c(as.character(pair_df$TF1), as.character(pair_df$Human$TF2)))
  gene_vec <- setNames(rep(0, nrow(gene_mat)), rownames(gene_mat))
  
  counts <- lapply(tfs, function(x) {
    tf_df <- filter(pair_df, TF1 == x & TF2 == x)
    pairs <- paste(tf_df$Row, tf_df$Col, sep = ":")
    gene_list <- gene_list[pairs]
    gene_count <- table(unlist(gene_list))
    gene_vec[names(gene_count)] <- gene_count
    return(gene_vec)
  })
  
  mat <- do.call(cbind, counts)
  colnames(mat) <- tfs
  return(mat)
}


count_hg <- overlap_counts(gene_list = sim_list$Human$Pval_Genes, pair_df = df_list$Human, gene_mat = dat$Binding$Human$Raw)
count_mm <- overlap_counts(gene_list = sim_list$Mouse$Pval_Genes, pair_df = df_list$Mouse, gene_mat = dat$Binding$Mouse$Raw)
count_ortho <- overlap_counts(gene_list = sim_list$Ortho$Pval_Genes, pair_df = df_list$Ortho, gene_mat = dat$Binding$Ortho$Raw)



# Plot
# ------------------------------------------------------------------------------


# Violin+boxplot of group differences for all experiments


p1a <- stat_vboxplot(df_list$Human, y_var = "Pval_Intersect", y_name = paste0("Top-", topn, " overlap"), title = "Human: P-values")
p1b <- stat_vboxplot(df_list$Human, y_var = "Upreg_Intersect", y_name = paste0("Top-", topn, " overlap"), title = "Human: Up-regulated")
p1c <- stat_vboxplot(df_list$Human, y_var = "Downreg_Intersect", y_name = paste0("Top-", topn, " overlap"), title = "Human: Down-regulated")

p2a <- stat_vboxplot(df_list$Mouse, y_var = "Pval_Intersect", y_name = paste0("Top-", topn, " overlap"), title = "Mouse: P-values")
p2b <- stat_vboxplot(df_list$Mouse, y_var = "Upreg_Intersect", y_name = paste0("Top-", topn, " overlap"), title = "Mouse: Up-regulated")
p2c <- stat_vboxplot(df_list$Mouse, y_var = "Downreg_Intersect", y_name = paste0("Top-", topn, " overlap"), title = "Mouse: Down-regulated")

p3a <- stat_vboxplot(df_list$Ortho, y_var = "Pval_Intersect", y_name = paste0("Top-", topn, " overlap"), title = "Ortho: P-values")
p3b <- stat_vboxplot(df_list$Ortho, y_var = "Upreg_Intersect", y_name = paste0("Top-", topn, " overlap"), title = "Ortho: Up-regulated")
p3c <- stat_vboxplot(df_list$Ortho, y_var = "Downreg_Intersect", y_name = paste0("Top-", topn, " overlap"), title = "Ortho: Down-regulated")


ggsave(plot_grid(p1a, p1b, p1c, nrow = 1),
       dpi = 300, device = "png", height = 8, width = 16,
       filename = paste0(plot_dir, "Vbplot_all_human_", date, ".png"))

ggsave(plot_grid(p2a, p2b, p2c, nrow = 1),
       dpi = 300, device = "png", height = 8, width = 16,
       filename = paste0(plot_dir, "Vbplot_all_mouse_", date, ".png"))

ggsave(plot_grid(p3a, p3b, p3c, nrow = 1),
       dpi = 300, device = "png", height = 8, width = 20,
       filename = paste0(plot_dir, "Vbplot_all_ortho_", date, ".png"))



# Violin+boxplot of group differences by TR


p4a <- tf_vboxplot(df_list$Human, y_var = "Pval_Intersect", y_name = paste0("Top-", topn, " overlap"), title_prefix = "Human")
p4b <- tf_vboxplot(df_list$Mouse, y_var = "Pval_Intersect", y_name = paste0("Top-", topn, " overlap"), title_prefix = "Mouse")
p4c <- tf_vboxplot(df_list$Ortho, y_var = "Pval_Intersect", y_name = paste0("Top-", topn, " overlap"), title_prefix = "Ortho")


ggsave(plot_grid(plotlist = p4a, nrow = 2),
       dpi = 300, device = "png", height = 12, width = 20,
       filename = paste0(plot_dir, "Vbplot_tf_human_intersect_", date, ".png"))

ggsave(plot_grid(plotlist = p4b, nrow = 2),
       dpi = 300, device = "png", height = 12, width = 20,
       filename = paste0(plot_dir, "Vbplot_tf_mouse_intersect_", date, ".png"))

ggsave(plot_grid(plotlist = p4c, nrow = 2),
       dpi = 300, device = "png", height = 12, width = 24,
       filename = paste0(plot_dir, "Vbplot_tf_ortho_intersect_", date, ".png"))



p5a <- dplot(df_list$Human, stat = "Pval_Intersect", stat_name = "Top-500 overlap", species = "Human")
p5b <- dplot(df_list$Mouse, stat = "Pval_Intersect", stat_name = "Top-500 overlap", species = "Mouse")
p5c <- dplot(df_list$Ortho, stat = "Pval_Intersect", stat_name = "Top-500 overlap", species = "Ortho")


ggsave(p5a, dpi = 300, device = "png", height = 6, width = 9,
       filename = paste0(plot_dir, "Densplot_pval-intersect_all_human_", date, ".png"))

ggsave(p5b, dpi = 300, device = "png", height = 6, width = 9,
       filename = paste0(plot_dir, "Densplot_pval-intersect_all_mouse_", date, ".png"))

ggsave(p5c, dpi = 300, device = "png", height = 6, width = 9,
       filename = paste0(plot_dir, "Densplot_pval-intersect_all_ortho_", date, ".png"))


# Heatmap


anno_colour = list(Symbol = tf_pal)
# Species = c(Human = "royalblue", Mouse = "goldenrod")  # pretty cramped using TF and species colours

pal_length <- 11
heatmap_pal <- viridis::magma(pal_length)


plot_heatmap <- function(plot_mat, 
                         chip_meta, 
                         perturb_meta,
                         pal_length = 11,
                         file) {
  
  color_min <- min(plot_mat, na.rm = TRUE)
  color_max <- max(plot_mat, na.rm = TRUE)
  color_breaks <- seq(color_min, color_max, length.out = pal_length)
  
  # annotations for heatmap
  row_meta <- chip_meta %>% 
    distinct(Experiment_ID, .keep_all = TRUE) %>% 
    filter(Experiment_ID %in% rownames(plot_mat)) %>% 
    mutate(Symbol = str_to_title(Symbol)) %>% 
    # select(Symbol, Species)
    select(Symbol)
  rownames(row_meta) <- rownames(plot_mat)
  
  col_meta <- perturb_meta %>% 
    filter(Experiment_ID %in% colnames(plot_mat)) %>% 
    mutate(Symbol = str_to_title(Symbol)) %>% 
    select(Symbol)
  # select(Symbol, Species)
  rownames(col_meta) <- colnames(plot_mat)  
  
  # tf row/col breaks
  row_breaks <- sapply(unique(row_meta$Symbol), function(x) {
    tail(which(row_meta$Symbol == x), n = 1)
  })
  
  col_breaks <- sapply(unique(col_meta$Symbol), function(x) {
    tail(which(col_meta$Symbol == x), n = 1)
  })
  
  
  pheatmap(
    plot_mat,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    show_rownames = FALSE,
    show_colnames = FALSE,
    border_color = NA,
    color = heatmap_pal,
    breaks = color_breaks,
    annotation_row = row_meta,
    annotation_col = col_meta,
    annotation_names_row = FALSE,
    annotation_names_col = FALSE,
    annotation_colors = anno_colour,
    gaps_row = row_breaks,
    gaps_col = col_breaks,
    height = 4,
    width = 6,
    filename = file)
  
}


plot_heatmap(plot_mat = sim_list$Human$Pval_Intersect, chip_meta = dat$Binding$Meta, perturb_meta = dat$Perturbation$Meta,
             file = paste0(plot_dir, "Human_Intersect_heatmap_", date, "_chip=QNlog_perturb=pval_top=", topn, "common_only=", common_arg, ".png"))

plot_heatmap(plot_mat = sim_list$Mouse$Pval_Intersect, chip_meta = dat$Binding$Meta, perturb_meta = dat$Perturbation$Meta,
             file = paste0(plot_dir, "Mouse_Intersect_heatmap_", date, "_chip=QNlog_perturb=pval_top=", topn, "common_only=", common_arg, ".png"))

plot_heatmap(plot_mat = sim_list$Ortho$Pval_Intersect, chip_meta = dat$Binding$Meta, perturb_meta = dat$Perturbation$Meta,
             file = paste0(plot_dir, "Ortho_Intersect_heatmap_", date, "_chip=QNlog_perturb=pval_top=", topn, "common_only=", common_arg, ".png"))