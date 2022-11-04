## This script looks at the similarity of perturbation experiments, using cor
## across all pairwise FCs (+/- abs), as well as metrics of overlap for the 
## topn scoring (by up/downreg as well as pval) genes
## -----------------------------------------------------------------------------

library(tidyverse)
library(parallel)
library(pheatmap)
library(cowplot)
library(RColorBrewer)
library(WGCNA)
source("R/setup-01_config.R")
source("R/utils/similarity_functions.R")
source("R/utils/plot_functions.R")

topn <- 500 # how many of the top genes to keep
common_arg <- TRUE  # should comparison only be done for mutually measured genes?
date <- "Apr2022"  # most recent data freeze
outfile <- paste0(scratch_dir, date, "_perturb_similarity.RDS")
plot_dir <- paste0(pplot_dir, "/Experiment_similarity/")

# Load meta and lists of perturb effect size matrices
meta <- read.delim(file =  paste0(meta_dir, "batch1_tfperturb_meta_final_", date, ".tsv"), stringsAsFactors = FALSE)
mlist_hg <- readRDS(paste0(pmat_dir, "human_list_perturb_matrix_", date, ".RDS"))
mlist_mm <- readRDS(paste0(pmat_dir, "mouse_list_perturb_matrix_", date, ".RDS"))
mlist_ortho <- readRDS(paste0(pmat_dir, "ortho_list_perturb_matrix_", date, ".RDS"))


# General workflow is to generate exp x exp matrices where elements represent
# different metrics of similarity/overlap between the experiments. These are then
# turned into dfs of the unique pairs for plotting/summaries
# 1) Pearson cor of fold changes
# 2) Pearson cor of absolute fold changes
# FOR EACH OF SORTING BY PVAL, UPREGULATED, AND DOWNREGULATED:
# 3) Length of the intersect of the topn overlap between experiments
# 4) Jaccard of the topn overlap
# 5) FET pval of the topn overlap
# 6) Log odds ratio of the topn overlap
#-------------------------------------------------------------------------------


# NOTE: slow! get_overlap_matrices needs to be sped up


if (!file.exists(outfile)) {
  
  sim_list <- list(
    Human = perturb_sim_list(fc_mat = mlist_hg$FC_mat, pval_mat = mlist_hg$Pval_mat, topn = topn),
    Mouse = perturb_sim_list(fc_mat = mlist_mm$FC_mat, pval_mat = mlist_mm$Pval_mat, topn = topn),
    Ortho = perturb_sim_list(fc_mat = mlist_ortho$FC_mat, pval_mat = mlist_ortho$Pval_mat, topn = topn))
  
  df_list <- list(
    Human = format_and_merge(sim_list$Human, gene_mat = mlist_hg$FC_mat, add_common = TRUE),
    Mouse = format_and_merge(sim_list$Mouse, gene_mat = mlist_mm$FC_mat, add_common = TRUE),
    Ortho = format_and_merge(sim_list$Ortho, gene_mat = mlist_ortho$FC_mat, add_common = TRUE))
  
  saveRDS(list(sim_list = sim_list, df_list = df_list), file = outfile)
  
} else {
  
  dat <- readRDS(outfile)
  df_list <- dat$df_list
  sim_list <- dat$sim_list
  rm(dat)
  
}


# Summarize/inspect group similarities (focus on what is in paper)
#-------------------------------------------------------------------------------


focus_cols <- c("Group", 
                "Pcor_Mean", "Pcor_Max",
                "Pcor_abs_Mean", "Pcor_abs_Max",
                "Pval_Intersect_Mean", "Pval_Intersect_Max",
                "Upreg_Intersect_Mean", "Upreg_Intersect_Max",
                "Downreg_Intersect_Mean", "Downreg_Intersect_Max")


all_summ <- lapply(df_list, get_summary)
tf_summ <- lapply(df_list, tf_summary)

# Group summaries by all experiments or TF-specific
lapply(all_summ, `[`, focus_cols)
lapply(tf_summ, function(tf) lapply(tf, `[`, focus_cols))


# split by df TF
tf_df <- lapply(df_list, function(x) {
  tfs <- unique(c(as.character(x$TF1), as.character(x$TF2)))
  tf_list <- lapply(tfs, function(y) {
    filter(x, as.character(TF1) == y & as.character(TF2) == y)
  })
  names(tf_list) <- tfs
  return(tf_list)
})


# top pairs by Pcor/Intersect and anti-correlated pairs (non-abs)
top_pcor <- lapply(df_list, function(x) arrange(x, desc(Pcor_abs)) %>% head(10))
top_apcor <- lapply(df_list, function(x) arrange(x, Pcor) %>% head(10))
top_inter <- lapply(df_list, function(x) arrange(x, desc(Pval_Intersect)) %>% head(10))

# Top experiments from different GEOs
diff_geo_hg <- filter(df_list$Human, as.character(GEO1) != as.character(GEO2))
diff_geo_mm <- filter(df_list$Mouse, as.character(GEO1) != as.character(GEO2))

# Example of negative correlation from same TR - typically Mecp2 OE and LoF
mecp2_1 <- "GSE126640_MECP2_Human_Overexpression"
mecp2_2 <- "GSE126640_MECP2_Human_Knockout"

# Inspect top examples from inter-species & intra-TR
top_ortho <- filter(df_list$Ortho, Group == "In_TR_out_species")


diff_df1 <- data.frame(
  exp1 = mlist_hg$FC_mat[, mecp2_1],
  exp2 = mlist_hg$FC_mat[, mecp2_2]
) %>%
  mutate(Diff = exp2 - exp1)


# Example of negative correlation from different TRs: 2 Neurod1 OE and Mecp2 KO
neurod1_1 <- "GSE104435_Neurod1_Mouse_Overexpression"
neurod1_2 <- "GSE135981_Neurod1_Mouse_Overexpression"
mecp2_3 <- "GSE124521_Mecp2_Mouse_Knockout"

diff_df2 <- data.frame(
  ND1_1 = mlist_mm$FC_mat[, neurod1_1],
  ND1_2 = mlist_mm$FC_mat[, neurod1_2],
  M = mlist_mm$FC_mat[, mecp2_3]
) %>%
  mutate(Diff1 = ND1_1 - M,
         Diff2 = ND1_2 - M)

plot(mlist_mm$FC_mat[, neurod1_1], mlist_mm$FC_mat[, mecp2_3])
plot(mlist_mm$FC_mat[, neurod1_2], mlist_mm$FC_mat[, mecp2_3])

intersect(
  get_topn_genes(mlist_mm$Pval_mat[, neurod1_1], topn, decrease = FALSE),
  get_topn_genes(mlist_mm$Pval_mat[, neurod1_2], topn, decrease = FALSE)
)


# Plots
# ------------------------------------------------------------------------------


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
p4c <- tf_vboxplot(df_list$Ortho, y_var = "Pval_Intersect", y_name = paste0("Top-", topn, " overlap"), title_prefix = "Ortho", ortho = TRUE)


ggsave(plot_grid(plotlist = p4a, nrow = 2),
       dpi = 300, device = "png", height = 12, width = 20,
       filename = paste0(plot_dir, "Vbplot_tf_human_intersect_", date, ".png"))

ggsave(plot_grid(plotlist = p4b, nrow = 2),
       dpi = 300, device = "png", height = 12, width = 20,
       filename = paste0(plot_dir, "Vbplot_tf_mouse_intersect_", date, ".png"))

ggsave(plot_grid(plotlist = p4c, nrow = 2),
       dpi = 300, device = "png", height = 12, width = 24,
       filename = paste0(plot_dir, "Vbplot_tf_ortho_intersect_", date, ".png"))



p5a <- dplot(df_list$Human, stat = "Pval_Intersect", stat_name = "Top-500 Overlap by P-value", species = "Human")
p5b <- dplot(df_list$Mouse, stat = "Pval_Intersect", stat_name = "Top-500 Overlap by P-value", species = "Mouse")
p5c <- dplot(df_list$Ortho, stat = "Pval_Intersect", stat_name = "Top-500 Overlap by P-value", species = "Ortho")

p5d <- dplot(df_list$Human, stat = "Pcor", stat_name = "Pearson correlation", species = "Human")
p5e <- dplot(df_list$Mouse, stat = "Pcor", stat_name = "Pearson correlation", species = "Mouse")
p5f <- dplot(df_list$Ortho, stat = "Pcor", stat_name = "Pearson correlation", species = "Ortho")


ggsave(p5a, dpi = 300, device = "png", height = 6, width = 8,
       filename = paste0(plot_dir, "Densplot_pval-intersect_all_human_", date, ".png"))

ggsave(p5b, dpi = 300, device = "png", height = 6, width = 8,
       filename = paste0(plot_dir, "Densplot_pval-intersect_all_mouse_", date, ".png"))

ggsave(p5c, dpi = 300, device = "png", height = 6, width = 8,
       filename = paste0(plot_dir, "Densplot_pval-intersect_all_ortho_", date, ".png"))


ggsave(p5d, dpi = 300, device = "png", height = 6, width = 8,
       filename = paste0(plot_dir, "Densplot_Pcor_human_", date, ".png"))

ggsave(p5e, dpi = 300, device = "png", height = 6, width = 8,
       filename = paste0(plot_dir, "Densplot_Pcor_mouse_", date, ".png"))

ggsave(p5f, dpi = 300, device = "png", height = 6, width = 8,
       filename = paste0(plot_dir, "Densplot_Pcor_ortho_", date, ".png"))


# format matrices for heatmaps: round, triangular, and remove diag

cor_list <- lapply(sim_list, `[[`, "Pcor")

cor_list <- lapply(cor_list, function(mat) {
  mat <- round(mat, 4)
  diag(mat) <- NA
  mat[upper.tri(mat)] <- NA
  return(mat)
})


# colours

tf_pal = c(
  Ascl1 = "#00188f",
  Hes1 = "#fff100",
  Mecp2 = "#009e49",
  Mef2c = "#00bcf2",
  Neurod1 = "#68217a",
  Pax6 = "#ffd8b1",
  Runx1 = "#e81123",
  Tcf4 = "#ff8c00"
)

tf_colour = list(Symbol = tf_pal)
# Species = c(Human = "royalblue", Mouse = "goldenrod")


pal_length <- 100
heatmap_pal <- colorRampPalette(c("#0571b0", "white", "#ca0020"))(pal_length)


plot_heatmap <- function(mat, 
                         meta, 
                         file, 
                         col_min = -1, 
                         col_max = 1,
                         ortho = FALSE) {
  
  # Anno for colours
  anno <- meta %>% 
    mutate(Symbol = str_to_title(Symbol)) %>% 
    # dplyr::select(Symbol, Species)
    dplyr::select(Symbol)
  
  # Because otherwise anno colour for single species
  # if(!ortho) anno$Species <- NULL
  
  rownames(anno) <- rownames(mat)
  
  color_breaks <- seq(col_min, col_max, length.out = pal_length)
  
  # Cut locations for gaps
  cuts <- sapply(unique(anno$Symbol), function(x) {
    tail(which(anno$Symbol == x), n = 1)
  }) 
  
  pheatmap(
    mat,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    show_rownames = FALSE,
    show_colnames = FALSE,
    na_col = "white",
    border_color = "white",
    color = heatmap_pal,
    breaks = color_breaks,
    annotation_row = anno,
    annotation_col = anno,
    annotation_names_row = FALSE,
    annotation_names_col = FALSE,
    annotation_colors = tf_colour,
    gaps_row = cuts,
    gaps_col = cuts,
    legend = TRUE,
    height = 4,
    width = 6,
    filename = file)
  
}


plot_heatmap(mat = cor_list$Mouse, 
             meta = filter(meta, Species == "Mouse"),
             file = paste0(plot_dir, "cor_heatmap_mouse_", date, ".png"))

plot_heatmap(mat = cor_list$Human, 
             meta = filter(meta, Species == "Human"),
             file = paste0(plot_dir, "cor_heatmap_human_", date, ".png"))

plot_heatmap(mat = cor_list$Ortho, 
             meta = meta, 
             ortho = TRUE,
             file = paste0(plot_dir, "cor_heatmap_ortho_", date, ".png"))
