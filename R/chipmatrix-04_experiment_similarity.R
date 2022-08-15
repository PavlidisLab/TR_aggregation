## This script looks at different measures of similarity between the batch1
## ChIP-seq experiments
## -----------------------------------------------------------------------------

library(tidyverse)
library(cowplot)
library(pheatmap)
library(RColorBrewer)
library(WGCNA)
source("~/regnetR/R/utils/similarity_functions.R")
source("~/regnetR/R/utils/plot_functions.R")

topn <- 500 # how many of the top genes to keep
date <- "Apr2022"  # most recent data freeze
plot_dir <- "~/Plots/Chipseq/Binding_similarity/"  # plot export dir
outfile <- paste0("~/scratch/R_objects/", date, "_chip_similarity_refseq.RDS")

# Loading ChIP-seq data
intensity_flag <- FALSE  # use ouyang scores with or without macs2 intensity
binary_dist <- 25e3  # distance threshold used for binary gene scores
min_peaks <- 100  # how many peaks were required to keep an experiment
peakset <- "idr"  # idr for TF (but overlap still for mecp2) or all overlap
chip_dir <- "~/Data/Annotated_objects/Bind_matrices/Encpipe/" # bind matrix path
chip_type <- "QN_log"  # which chip processing scheme to use
meta <- read.delim(paste0("~/Data/Metadata/Chipseq/batch1_chip_meta_final_", date, ".tsv"), stringsAsFactors = FALSE)
chip_hg <- readRDS(paste0(chip_dir, "Human_refseq_", date, "_processed_bindmat_list_minpeak=", min_peaks, "_ouyang_dc=5000_intensity=", intensity_flag, "_binary=", binary_dist/1e3, "kb_peakset=", peakset, ".RDS"))
chip_mm <- readRDS(paste0(chip_dir, "Mouse_refseq_", date, "_processed_bindmat_list_minpeak=", min_peaks, "_ouyang_dc=5000_intensity=", intensity_flag, "_binary=", binary_dist/1e3, "kb_peakset=", peakset, ".RDS"))
chip_ortho <- readRDS(paste0(chip_dir, "Ortho_refseq_", date, "_processed_bindmat_list_minpeak=", min_peaks, "_ouyang_dc=5000_intensity=", intensity_flag, "_binary=", binary_dist/1e3, "kb_peakset=", peakset, ".RDS"))


# General workflow is to generate exp x exp matrices where elements represent
# different metrics of similarity/overlap between the experiments. These are then
# turned into dfs of the unique pairs for plotting/summaries
# 1) Pearson cor of continuous bind scores
# 2) Length of the intersect of the topn overlap between experiments
# 3) Jaccard of the topn overlap
# 4) FET pval of the topn overlap
# 5) Log odds ratio of the topn overlap
#-------------------------------------------------------------------------------


# NOTE: slow! get_overlap_matrices needs to be sped up


if (!file.exists(outfile)) {
  
  sim_list <- list(
    Human = chip_sim_list(chip_hg[[chip_type]], topn),
    Mouse = chip_sim_list(chip_mm[[chip_type]], topn),
    Ortho = chip_sim_list(chip_ortho[[chip_type]], topn))
  
  df_list <- lapply(sim_list, format_and_merge)
  
  saveRDS(list(sim_list = sim_list, df_list = df_list), file = outfile)
  
} else {
  
  dat <- readRDS(outfile)
  df_list <- dat$df_list
  sim_list <- dat$sim_list
  rm(dat)
  
}


# Summarize/inspect group similarities (focus on what is in paper)
#-------------------------------------------------------------------------------


focus_cols <- c("Group", "Pcor_Mean", "Pcor_Max", "Intersect_Mean", "Intersect_Max")

all_summ <- lapply(df_list, get_summary)
tf_summ <- lapply(df_list, tf_summary)

# Group summaries by all experiments or TF-specific
lapply(all_summ, `[`, focus_cols)
lapply(tf_summ, function(tf) lapply(tf, `[`, focus_cols))

# top pairs by Pcor/Intersect
top_pcor <- lapply(df_list, function(x) arrange(x, desc(Pcor)) %>% head(10))
top_inter <- lapply(df_list, function(x) arrange(x, desc(Intersect)) %>% head(10))

# split by df TF
tf_df <- lapply(df_list, function(x) {
  tfs <- unique(c(as.character(x$TF1), as.character(x$TF2)))
  tf_list <- lapply(tfs, function(y) {
    filter(x, as.character(TF1) == y & as.character(TF2) == y)
  })
  names(tf_list) <- tfs
  return(tf_list)
})


# Top experiments from different GEOs
diff_geo_hg <- filter(df_list$Human, as.character(GEO1) != as.character(GEO2))
diff_geo_mm <- filter(df_list$Mouse, as.character(GEO1) != as.character(GEO2))


# Point out intra-vs-inter HES1 comparisons: higher K562 and MCF-7 similarity
# between inter-HES1 pairs than intra-HES1
hes1_intra_hg <- df_list$Human %>%
  filter((as.character(TF1) == "Hes1" | as.character(TF2) == "Hes1") & Group == "In_TR_in_species")

hes1_inter_hg <- df_list$Human %>%
  filter((as.character(TF1) == "Hes1" | as.character(TF2) == "Hes1") & Group == "Out") %>% 
  arrange(desc(Pcor))

k562 <- filter(meta, Cell_Type == "K562")$Experiment_ID
mcf7 <- filter(meta, Cell_Type == "MCF-7")$Experiment_ID

hes1_k562 <- filter(hes1_inter_hg, Row %in% k562 & Col %in% k562) 
hes1_mcf7 <- filter(hes1_inter_hg, Row %in% mcf7 & Col %in% mcf7) 

# Runx1 dominates cross species intra-TF comparisons

top_ortho <- df_list$Ortho %>% 
  filter(Group == "In_TR_out_species") %>% 
  arrange(desc(Intersect))

mean_pcor <- mean(df_list$Ortho[df_list$Ortho$Group == "In_TR_out_species", "Pcor"])
mean_inter <- mean(df_list$Ortho[df_list$Ortho$Group == "In_TR_out_species", "Intersect"])

pcor_gt <- df_list$Ortho %>% 
  filter(Group == "In_TR_out_species" & Pcor > mean_pcor) %>% 
  dplyr::count(TF1)

inter_gt <- df_list$Ortho %>% 
  filter(Group == "In_TR_out_species" & Intersect > mean_inter) %>% 
  dplyr::count(TF1)


# Mef2c example of high overlap with low correlation. Experiments have low
# peak counts and thus sparsely bound. Means that experiments with no binding
# are included in the top 500, resulting in overlaps of zero bound genes.

df <- data.frame(
  Symbol = rownames(chip_ortho$QN_log),
  # Exp1 = chip_ortho$QN_log[, "GSE112315_Mef2c_Mouse_Day0"], 
  Exp1 = chip_ortho$QN_log[, "GSE139509_Mecp2_Mouse_Ab1-Mecp2-KO_HISTONE"],
  Exp2 = chip_ortho$QN_log[, "GSE32644_MEF2C_Human_Control"])

cor(df$Exp1, df$Exp2)
plot(df$Exp1, df$Exp2)
df <- df[df$Exp1 > 0.2 & df$Exp2 > 0.2, ]

# Mecp2 example of decent raw/qn_log cor/overlap, but non-existant binary cor. 
# Near max peaks, resulting in binary binding for almost every gene
df <- data.frame(Exp1 = chip_ortho$Binary[, "GSE122364_MECP2_Human_GSK343_HISTONE"], 
                 Exp2 = chip_ortho$Binary[, "GSE125585_Mecp2_Mouse_Control_HISTONE"])
df <- arrange(df, desc(Exp1))
cor(df$Exp1, df$Exp2)
plot(df$Exp1, df$Exp2)
sum(df$Exp1 == 1)
sum(df$Exp2 == 1)


# Plot
#-------------------------------------------------------------------------------


# Violin+boxplot of group differences for all experiments


p1a <- stat_vboxplot(df_list$Human, y_var = "Pcor", y_name = "Pearson correlation", title = "Human ChIP-seq")
p1b <- stat_vboxplot(df_list$Human, y_var = "Intersect", y_name = paste0("Top-", topn, " overlap"), title = "Human ChIP-seq")

p1c <- stat_vboxplot(df_list$Mouse, y_var = "Pcor", y_name = "Pearson correlation", title = "Mouse ChIP-seq")
p1d <- stat_vboxplot(df_list$Mouse, y_var = "Intersect", y_name = paste0("Top-", topn, " overlap"), title = "Mouse ChIP-seq")

p1e <- stat_vboxplot(df_list$Ortho, y_var = "Pcor", y_name = "Pearson correlation", title = "Ortho ChIP-seq", ortho = TRUE)
p1f <- stat_vboxplot(df_list$Ortho, y_var = "Intersect", y_name = paste0("Top-", topn, " overlap"), title = "Ortho ChIP-seq", ortho = TRUE)


ggsave(plot_grid(p1a, p1b, nrow = 1),
       dpi = 300, device = "png", height = 8, width = 12,
       filename = paste0(plot_dir, "Vbplot_all_human_", date, ".png"))

ggsave(plot_grid(p1c, p1d, nrow = 1),
       dpi = 300, device = "png", height = 8, width = 12,
       filename = paste0(plot_dir, "Vbplot_all_mouse_", date, ".png"))

ggsave(plot_grid(p1e, p1f, nrow = 1),
       dpi = 300, device = "png", height = 8, width = 12,
       filename = paste0(plot_dir, "Vbplot_all_ortho_", date, ".png"))



# Violin+boxplot of group differences by TR


p2a <- tf_vboxplot(df_list$Human, y_var = "Intersect", y_name = paste0("Top-", topn, " overlap"), title_prefix = "Human")
p2b <- tf_vboxplot(df_list$Mouse, y_var = "Intersect", y_name = paste0("Top-", topn, " overlap"), title_prefix = "Mouse")
p2c <- tf_vboxplot(df_list$Ortho, y_var = "Intersect", y_name = paste0("Top-", topn, " overlap"), title_prefix = "Ortho", ortho = TRUE)


ggsave(plot_grid(plotlist = p2a, nrow = 2),
       dpi = 300, device = "png", height = 12, width = 20,
       filename = paste0(plot_dir, "Vbplot_tf_human_intersect_", date, ".png"))

ggsave(plot_grid(plotlist = p2b, nrow = 2),
       dpi = 300, device = "png", height = 12, width = 20,
       filename = paste0(plot_dir, "Vbplot_tf_mouse_intersect_", date, ".png"))

ggsave(plot_grid(plotlist = p2c, nrow = 2),
       dpi = 300, device = "png", height = 12, width = 24,
       filename = paste0(plot_dir, "Vbplot_tf_ortho_intersect_", date, ".png"))


# density plots of pcor. keep consistent xscale

xmin <- -0.4
xmax <- 0.9

p3a <- dplot(df_list$Human, stat = "Pcor", stat_name = "Pearson correlation", species = "Human") + xlim(xmin, xmax)
p3b <- dplot(df_list$Mouse, stat = "Pcor", stat_name = "Pearson correlation", species = "Mouse") + xlim(xmin, xmax)
p3c <- dplot(df_list$Ortho, stat = "Pcor", stat_name = "Pearson correlation", species = "Ortho") + xlim(xmin, xmax)


ggsave(p3a, dpi = 300, device = "png", height = 6, width = 8,
       filename = paste0(plot_dir, "Densplot_all_human_", date, ".png"))

ggsave(p3b, dpi = 300, device = "png", height = 6, width = 8,
       filename = paste0(plot_dir, "Densplot_all_mouse_", date, ".png"))

ggsave(p3c, dpi = 300, device = "png", height = 6, width = 8,
       filename = paste0(plot_dir, "Densplot_all_ortho_", date, ".png"))


# density plot of ortho with and w/o RUNX1

p4a <- dplot(df_list$Ortho, stat = "Intersect", stat_name = "Count of top 500 overlap", species = "Ortho")
p4a <- p4a + ggtitle("Ortho including RUNX1") + theme(legend.position = "none")

p4b <- dplot(filter(df_list$Ortho, TF1 != "Runx1" & TF2 != "Runx1"),
             stat = "Intersect", stat_name = "Count of top 500 overlap", species = "Ortho")
p4b <- p4b + ggtitle("Ortho excluding RUNX1")


p4c <- dplot(df_list$Ortho, stat = "Pcor", stat_name = "Pearson correlation", species = "Ortho")
p4c <- p4c + ggtitle("Ortho including RUNX1") + theme(legend.position = "none")

p4d <- dplot(filter(df_list$Ortho, TF1 != "Runx1" & TF2 != "Runx1"),
             stat = "Pcor", stat_name = "Pearson correlation", species = "Ortho")
p4d <- p4d + ggtitle("Ortho excluding RUNX1")


ggsave(plot_grid(p4a, p4b), dpi = 300, device = "png", height = 8, width = 20,
       filename = paste0(plot_dir, "Densplot_inter_runx1_ortho_", date, ".png"))

ggsave(plot_grid(p4c, p4d), dpi = 300, device = "png", height = 8, width = 20,
       filename = paste0(plot_dir, "Densplot_cor_runx1_ortho_", date, ".png"))