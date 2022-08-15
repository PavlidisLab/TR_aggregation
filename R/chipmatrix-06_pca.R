## This script performs PCA on the ChIP-seq bind matrices
## TODO: plot functions for cor heatmap and technical var scatter
## TODO: fix ortho PC scatter showing TR by species
## -----------------------------------------------------------------------------

library(tidyverse)
library(factoextra)
library(grid)
library(gridExtra)
source("~/regnetR/R/utils/plot_functions.R")

date <- "Apr2022"  # most recent data freeze
intensity_flag <- FALSE  # use ouyang scores with or without macs2 intensity
binary_dist <- 25e3  # distance threshold used for binary gene scores
min_peaks <- 100  # how many peaks were required to keep an experiment
peakset <- "idr"  # idr for TF (but overlap still for mecp2) or all overlap
plot_dir <- "~/Plots/Chipseq/PCA/"  # plot export dir
mat_dir <- "~/Data/Annotated_objects/Bind_matrices/Encpipe/" # bind matrix path

# batch 1 ChIP-seq meta
meta <- read.delim(paste0("~/Data/Metadata/Chipseq/batch1_chip_meta_final_", date, ".tsv"), stringsAsFactors = FALSE)

# batch 1 binding matrices
bind_hg <- readRDS(paste0(mat_dir, "Human_refseq_", date, "_processed_bindmat_list_minpeak=", min_peaks, "_ouyang_dc=5000_intensity=", intensity_flag, "_binary=", binary_dist/1e3, "kb_peakset=", peakset, ".RDS"))
bind_mm <- readRDS(paste0(mat_dir, "Mouse_refseq_", date, "_processed_bindmat_list_minpeak=", min_peaks, "_ouyang_dc=5000_intensity=", intensity_flag, "_binary=", binary_dist/1e3, "kb_peakset=", peakset, ".RDS"))
bind_ortho <- readRDS(paste0(mat_dir, "Ortho_refseq_", date, "_processed_bindmat_list_minpeak=", min_peaks, "_ouyang_dc=5000_intensity=", intensity_flag, "_binary=", binary_dist/1e3, "kb_peakset=", peakset, ".RDS"))

tfs_hg <- unique(str_to_title(meta$Symbol))
tfs_mm <- unique(filter(meta, Species == "Mouse")$Symbol)


# PCA
#-------------------------------------------------------------------------------


remove_nonvar <- function(mat) {
  # Remove rows of matrix that show no variation
  no_sd <- which(apply(mat, 1, sd) == 0)
  if (length(no_sd) > 1) {
    mat <- mat[-no_sd, ]
  }
  return(mat)
}


pca_and_var <- function(mat) {
  
  # Performs PCA with prcomp and returns list of the resulting 
  # object as well as the variance explained
  
  # prcomp expects samples as rows, features (genes) as columns so transpose
  pcmat <- prcomp(t(mat), scale = TRUE)
  
  # variance explained by the PCs
  prc_var <- pcmat$sdev ^ 2
  var_explained <- round(prc_var / sum(prc_var) * 100, 2)
  cumvar_explained <- cumsum(var_explained)/sum(var_explained)
  return(list(PC = pcmat, 
              Var_explained = var_explained, 
              Cumvar_explained = cumvar_explained))
}


pc_hg <- lapply(bind_hg, function(x) pca_and_var(remove_nonvar(x)))
pc_mm <- lapply(bind_mm, function(x) pca_and_var(remove_nonvar(x)))
pc_ortho <- lapply(bind_ortho, function(x) pca_and_var(remove_nonvar(x)))


# Variance explained by the PCs

lapply(pc_hg, `[[`, "Var_explained")
lapply(pc_hg, `[[`, "Cumvar_explained")

lapply(pc_mm, `[[`, "Var_explained")
lapply(pc_mm, `[[`, "Cumvar_explained")

lapply(pc_ortho, `[[`, "Var_explained")
lapply(pc_ortho, `[[`, "Cumvar_explained")


# Inspect the PC scores (samples in PC space) and loading scores (gene contributions)
#-------------------------------------------------------------------------------


# top genes

sort(abs(pc_hg$QN_log$PC$rotation[, 1]), decreasing = TRUE)[1:50]
sort(abs(pc_hg$QN_log$PC$rotation[, 2]), decreasing = TRUE)[1:50]

sort(abs(pc_mm$QN_log$PC$rotation[, 1]), decreasing = TRUE)[1:50]
sort(abs(pc_mm$QN_log$PC$rotation[, 2]), decreasing = TRUE)[1:50]

sort(abs(pc_ortho$QN_log$PC$rotation[, 1]), decreasing = TRUE)[1:50]
sort(abs(pc_ortho$QN_log$PC$rotation[, 2]), decreasing = TRUE)[1:50]

# top samples

sort(abs(pc_hg$QN_log$PC$x[, 1]), decreasing = TRUE)[1:50]
sort(abs(pc_hg$QN_log$PC$x[, 2]), decreasing = TRUE)[1:50]

sort(abs(pc_mm$QN_log$PC$x[, 1]), decreasing = TRUE)[1:50]
sort(abs(pc_mm$QN_log$PC$x[, 2]), decreasing = TRUE)[1:50]

sort(abs(pc_ortho$QN_log$PC$x[, 1]), decreasing = TRUE)[1:50]
sort(abs(pc_ortho$QN_log$PC$x[, 2]), decreasing = TRUE)[1:50]


# Correlate PCs with technical variables
# NOTE: Mecp2 Histone pipeline did not produce all QC metrics... so when doing
# pairiwise cor, there will be lots of drop outs for Mecp2 samples. 
#-------------------------------------------------------------------------------


num_pcs <- 10


pc_df_hg <- lapply(pc_hg, function(x) {
  # merge with meta
  df <- as.data.frame(x$PC$x[, 1:num_pcs]) %>%
    rownames_to_column("Experiment_ID") %>%
    merge(meta, sort = FALSE) %>% 
    mutate(Symbol = str_to_title(Symbol))
})


pc_df_mm <- lapply(pc_mm, function(x) {
  # merge with run QC df
  df <- as.data.frame(x$PC$x[, 1:num_pcs]) %>%
    rownames_to_column("Experiment_ID") %>%
    merge(meta, sort = FALSE)
  
})


pc_df_ortho <- lapply(pc_ortho, function(x) {
  # merge with run QC df
  df <- as.data.frame(x$PC$x[, 1:num_pcs]) %>%
    rownames_to_column("Experiment_ID") %>%
    merge(meta, sort = FALSE) %>% 
    mutate(Symbol = str_to_title(Symbol))
  
})

# only want numeric variables
nums <- unlist(lapply(pc_df_hg$Binary, is.numeric))  


# correlate

pc_cor_hg <- lapply(pc_df_hg, function(x) {
  round(cor(x[, nums], use = "pairwise.complete.obs"), 2)
})

pc_cor_mm <- lapply(pc_df_mm, function(x) {
  round(cor(x[, nums], use = "pairwise.complete.obs"), 2)
})

pc_cor_ortho <- lapply(pc_df_ortho, function(x) {
  round(cor(x[, nums], use = "pairwise.complete.obs"), 2)
})


# Plots
#-------------------------------------------------------------------------------


# scree plot - focus on quantile norm log (used for analysis)

p1a <- factoextra::fviz_screeplot(pc_hg$QN_log$PC, addlabels = TRUE, ylim = c(0, 50))
p1b <- factoextra::fviz_screeplot(pc_mm$QN_log$PC, addlabels = TRUE, ylim = c(0, 50))
p1c <- factoextra::fviz_screeplot(pc_ortho$QN_log$PC, addlabels = TRUE, ylim = c(0, 50))


# PC scatter plots - looking at pairwise PCs, colored by TFs


pc_scatter <- function(df,
                       pc_list,
                       pc_x, 
                       pc_y,
                       title) {
  
  ggplot(df, aes(x = !!sym(paste0("PC", pc_x)),
                 y = !!sym(paste0("PC", pc_y)),
                 fill = Symbol)) +
    geom_point(colour = "black", shape = 21, size = 3) +
    xlab(paste0("PC", pc_x, "(", pc_list$Var_explained[pc_x], "%)")) +
    ylab(paste0("PC", pc_y, "(", pc_list$Var_explained[pc_y], "%)")) +
    ggtitle(title) +
    theme_classic() +
    theme(axis.title = element_text(size = 15)) +
    scale_fill_manual(values = tf_pal)
  
}



# Human

for (i in names(pc_df_hg)) {
  
  plot_df <- pc_df_hg[[i]]
  pc_list <- pc_hg[[i]]
  
  plist <- lapply(1:(num_pcs-1), function(x) {
    pc_scatter(plot_df, pc_list, pc_x = x, pc_y = x + 1, title = "Human")
  })
  
  ggsave(
    filename = paste0(plot_dir, "PCplots_by_TF_", date, "_Human_", i, ".pdf"), 
    plot = marrangeGrob(plist, nrow = 1, ncol = 1), 
    width = 7, height = 7
  )
}


# Mouse

for (i in names(pc_df_mm)) {
  
  plot_df <- pc_df_mm[[i]]
  pc_list <- pc_mm[[i]]
  
  plist <- lapply(1:(num_pcs-1), function(x) {
    pc_scatter(plot_df, pc_list, pc_x = x, pc_y = x + 1, title = "Mouse")
  })
  
  ggsave(
    filename = paste0(plot_dir, "PCplots_by_TF_", date, "_Mouse_", i, ".pdf"), 
    plot = marrangeGrob(plist, nrow = 1, ncol = 1), 
    width = 7, height = 7
  )
}


# Ortho

for (i in names(pc_df_ortho)) {
  
  plot_df <- pc_df_ortho[[i]]
  pc_list <- pc_ortho[[i]]
  
  plist <- lapply(1:(num_pcs-1), function(x) {
    pc_scatter(plot_df, pc_list, pc_x = x, pc_y = x + 1, title = "Ortho")
  })
  
  ggsave(
    filename = paste0(plot_dir, "PCplots_by_TF_", date, "_Ortho", i, ".pdf"), 
    plot = marrangeGrob(plist, nrow = 1, ncol = 1), 
    width = 7, height = 7
  )
}


# heatmap of pc correlations


# set colour range
color_min <- -1
color_max <- 1
pal_length <- 100
heatmap_pal <- colorRampPalette(c("#0571b0", "white", "#ca0020"))(pal_length)
color_breaks <- seq(color_min, color_max, length.out = pal_length)



# Human

for (i in names(pc_cor_hg)) {
  
  plot_mat <- 
    as.matrix(pc_cor_hg[[i]][(num_pcs+1):nrow(pc_cor_hg[[i]]), 1:num_pcs])
  
  pheatmap(
    plot_mat,
    color = heatmap_pal,
    breaks = color_breaks,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    na_col = "black",
    border_color = "black",
    display_numbers = TRUE,
    number_color = "black",
    fontsize_number = 14,
    fontsize_row = 14,
    fontsize_col = 14,
    main = paste0("Human ChIP-seq experiments - ", i),
    height = 10,
    width = 12,
    filename = paste0(plot_dir, "Human_", date, "_", i, "_pcmat_cor_technical_vars.png")
  )
}

# Mouse

for (i in names(pc_cor_mm)) {
  
  plot_mat <- 
    as.matrix(pc_cor_mm[[i]][(num_pcs+1):nrow(pc_cor_mm[[i]]), 1:num_pcs])
  
  pheatmap(
    plot_mat,
    color = heatmap_pal,
    breaks = color_breaks,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    na_col = "black",
    border_color = "black",
    display_numbers = TRUE,
    number_color = "black",
    fontsize_number = 12,
    fontsize_row = 14,
    fontsize_col = 14,
    main = paste0("Mouse ChIP-seq experiments - ", i),
    height = 10,
    width = 12,
    filename = paste0(plot_dir, "Mouse_", date, "_", i, "_pcmat_cor_technical_vars.png")
  )
}


# Ortho

for (i in names(pc_cor_ortho)) {
  
  plot_mat <- 
    as.matrix(pc_cor_ortho[[i]][(num_pcs+1):nrow(pc_cor_ortho[[i]]), 1:num_pcs])
  
  pheatmap(
    plot_mat,
    color = heatmap_pal,
    breaks = color_breaks,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    na_col = "black",
    border_color = "black",
    display_numbers = TRUE,
    number_color = "black",
    fontsize_number = 12,
    fontsize_row = 14,
    fontsize_col = 14,
    main = paste0("All ChIP-seq experiments - ", i),
    height = 10,
    width = 12,
    filename = paste0(plot_dir, "Ortho_", date, "_", i, "_pcmat_cor_technical_vars.png")
  )
}


# Explore PC 1 vs 2 colored by technical factors


# Human


for (i in names(pc_df_hg)) {
  
  plot_df <- pc_df_hg[[i]]
  var_explained <- pc_hg[[i]]$Var_explained
  
  # presence of input
  p1 <- 
    mutate(plot_df, Has_input = Count_input > 0) %>% 
    ggplot(., aes(x = PC1, y = PC2, fill = Has_input)) +
    geom_point(colour = "black", shape = 21, size = 3, alpha = 0.5) +
    coord_fixed() +
    xlab(paste0("PC1 (", var_explained[1], "%)")) +
    ylab(paste0("PC2 (", var_explained[2], "%)")) +
    ggtitle(paste0("Human ChIP-seq experiments - ", i)) +
    theme_classic() +
    theme(axis.title = element_text(size = 15),
          legend.text = element_text(size = 15),
          legend.title = element_text(size = 15),
          legend.justification = c(0,0.5)) +
    scale_fill_manual(values = c("TRUE" = "black",
                                 "FALSE" = "red"))
  
  # count of reproducible peaks
  p2 <- 
    ggplot(plot_df, aes(x = PC1, y = PC2, fill = log10(N_peaks))) +
    geom_point(colour = "black", shape = 21, size = 3) +
    coord_fixed() +
    xlab(paste0("PC1 (", var_explained[1], "%)")) +
    ylab(paste0("PC2 (", var_explained[2], "%)")) +
    ggtitle(paste0("Human ChIP-seq experiments - ", i)) +
    theme_classic() +
    theme(axis.title = element_text(size = 15),
          legend.text = element_text(size = 15),
          legend.title = element_text(size = 15),
          legend.justification = c(0,0.5)) +
    scale_fill_viridis_c(option = "magma")
  
  
  # count of average experiment reads
  p3 <- 
    ggplot(plot_df, aes(x = PC1, y = PC2, fill = log10(Avg_exp_mapped_reads_nodup))) +
    geom_point(colour = "black", shape = 21, size = 3) +
    coord_fixed() +
    xlab(paste0("PC1 (", var_explained[1], "%)")) +
    ylab(paste0("PC2 (", var_explained[2], "%)")) +
    ggtitle(paste0("Human ChIP-seq experiments - ", i)) +
    theme_classic() +
    theme(axis.title = element_text(size = 15),
          legend.text = element_text(size = 15),
          legend.title = element_text(size = 15),
          legend.justification = c(0,0.5))
  
  # average NSC
  p4 <- 
    ggplot(plot_df, aes(x = PC1, y = PC2, fill = Avg_NSC)) +
    geom_point(colour = "black", shape = 21, size = 3) +
    coord_fixed() +
    xlab(paste0("PC1 (", var_explained[1], "%)")) +
    ylab(paste0("PC2 (", var_explained[2], "%)")) +
    ggtitle(paste0("Human ChIP-seq experiments - ", i)) +
    theme_classic() +
    theme(axis.title = element_text(size = 15),
          legend.text = element_text(size = 15),
          legend.title = element_text(size = 15),
          legend.justification = c(0,0.5)) +
    scale_fill_viridis_c()

  ggsave(
    filename = paste0(plot_dir, "PCplots_by_techvar_", date, "_Human_", i, ".pdf"),
    plot = marrangeGrob(list(p1, p2, p3, p4), nrow = 1, ncol = 1), 
    width = 7, height = 7
  )
}


# Mouse


for (i in names(pc_df_mm)) {
  
  plot_df <- pc_df_mm[[i]]
  var_explained <- pc_mm[[i]]$Var_explained
  
  # presence of input
  p1 <- 
    mutate(plot_df, Has_input = Count_input > 0) %>% 
    ggplot(., aes(x = PC1, y = PC2, fill = Has_input)) +
    geom_point(colour = "black", shape = 21, size = 3, alpha = 0.5) +
    coord_fixed() +
    xlab(paste0("PC1 (", var_explained[1], "%)")) +
    ylab(paste0("PC2 (", var_explained[2], "%)")) +
    ggtitle(paste0("Mouse ChIP-seq experiments - ", i)) +
    theme_classic() +
    theme(axis.title = element_text(size = 15),
          legend.text = element_text(size = 15),
          legend.title = element_text(size = 15),
          legend.justification = c(0,0.5)) +
    scale_fill_manual(values = c("TRUE" = "black",
                                 "FALSE" = "red"))
  
  # count of reproducible peaks
  p2 <- 
    ggplot(plot_df, aes(x = PC1, y = PC2, fill = log10(N_peaks))) +
    geom_point(colour = "black", shape = 21, size = 3) +
    coord_fixed() +
    xlab(paste0("PC1 (", var_explained[1], "%)")) +
    ylab(paste0("PC2 (", var_explained[2], "%)")) +
    ggtitle(paste0("Mouse ChIP-seq experiments - ", i)) +
    theme_classic() +
    theme(axis.title = element_text(size = 15),
          legend.text = element_text(size = 15),
          legend.title = element_text(size = 15),
          legend.justification = c(0,0.5)) +
    scale_fill_viridis_c(option = "magma")
  
  
  # count of average experiment reads
  p3 <- 
    ggplot(plot_df, aes(x = PC1, y = PC2, fill = log10(Avg_exp_mapped_reads_nodup))) +
    geom_point(colour = "black", shape = 21, size = 3) +
    coord_fixed() +
    xlab(paste0("PC1 (", var_explained[1], "%)")) +
    ylab(paste0("PC2 (", var_explained[2], "%)")) +
    ggtitle(paste0("Mouse ChIP-seq experiments - ", i)) +
    theme_classic() +
    theme(axis.title = element_text(size = 15),
          legend.text = element_text(size = 15),
          legend.title = element_text(size = 15),
          legend.justification = c(0,0.5))
  
  # average NSC
  p4 <- 
    ggplot(plot_df, aes(x = PC1, y = PC2, fill = Avg_NSC)) +
    geom_point(colour = "black", shape = 21, size = 3) +
    coord_fixed() +
    xlab(paste0("PC1 (", var_explained[1], "%)")) +
    ylab(paste0("PC2 (", var_explained[2], "%)")) +
    ggtitle(paste0("Mouse ChIP-seq experiments - ", i)) +
    theme_classic() +
    theme(axis.title = element_text(size = 15),
          legend.text = element_text(size = 15),
          legend.title = element_text(size = 15),
          legend.justification = c(0,0.5)) +
    scale_fill_viridis_c()
  
  ggsave(
    filename = paste0(plot_dir, "PCplots_by_techvar_", date, "_Mouse_", i, ".pdf"),
    plot = marrangeGrob(list(p1, p2, p3, p4), nrow = 1, ncol = 1), 
    width = 7, height = 7
  )
}


# Ortho


for (i in names(pc_df_ortho)) {
  
  plot_df <- pc_df_ortho[[i]]
  var_explained <- pc_ortho[[i]]$Var_explained
  
  # presence of input
  p1 <- 
    mutate(plot_df, Has_input = Count_input > 0) %>% 
    ggplot(., aes(x = PC1, y = PC2, fill = Has_input)) +
    geom_point(colour = "black", shape = 21, size = 3, alpha = 0.5) +
    coord_fixed() +
    xlab(paste0("PC1 (", var_explained[1], "%)")) +
    ylab(paste0("PC2 (", var_explained[2], "%)")) +
    ggtitle(paste0("All ChIP-seq experiments - ", i)) +
    theme_classic() +
    theme(axis.title = element_text(size = 15),
          legend.text = element_text(size = 15),
          legend.title = element_text(size = 15),
          legend.justification = c(0,0.5)) +
    scale_fill_manual(values = c("TRUE" = "black",
                                 "FALSE" = "red"))
  
  # count of reproducible peaks
  p2 <- 
    ggplot(plot_df, aes(x = PC1, y = PC2, fill = log10(N_peaks))) +
    geom_point(colour = "black", shape = 21, size = 3) +
    coord_fixed() +
    xlab(paste0("PC1 (", var_explained[1], "%)")) +
    ylab(paste0("PC2 (", var_explained[2], "%)")) +
    ggtitle(paste0("All ChIP-seq experiments - ", i)) +
    theme_classic() +
    theme(axis.title = element_text(size = 15),
          legend.text = element_text(size = 15),
          legend.title = element_text(size = 15),
          legend.justification = c(0,0.5)) +
    scale_fill_viridis_c(option = "magma")
  
  
  # count of average experiment reads
  p3 <- 
    ggplot(plot_df, aes(x = PC1, y = PC2, fill = log10(Avg_exp_mapped_reads_nodup))) +
    geom_point(colour = "black", shape = 21, size = 3) +
    coord_fixed() +
    xlab(paste0("PC1 (", var_explained[1], "%)")) +
    ylab(paste0("PC2 (", var_explained[2], "%)")) +
    ggtitle(paste0("All ChIP-seq experiments - ", i)) +
    theme_classic() +
    theme(axis.title = element_text(size = 15),
          legend.text = element_text(size = 15),
          legend.title = element_text(size = 15),
          legend.justification = c(0,0.5))
  
  # average NSC
  p4 <- 
    ggplot(plot_df, aes(x = PC1, y = PC2, fill = Avg_NSC)) +
    geom_point(colour = "black", shape = 21, size = 3) +
    coord_fixed() +
    xlab(paste0("PC1 (", var_explained[1], "%)")) +
    ylab(paste0("PC2 (", var_explained[2], "%)")) +
    ggtitle(paste0("All ChIP-seq experiments - ", i)) +
    theme_classic() +
    theme(axis.title = element_text(size = 15),
          legend.text = element_text(size = 15),
          legend.title = element_text(size = 15),
          legend.justification = c(0,0.5)) +
    scale_fill_viridis_c()
  
  ggsave(
    filename = paste0(plot_dir, "PCplots_by_techvar_", date, "_Ortho_", i, ".pdf"),
    plot = marrangeGrob(list(p1, p2, p3, p4), nrow = 1, ncol = 1), 
    width = 7, height = 7
  )
}


# Ortho PC plot by species+TF
# TODO: can't get 2 group legend to work https://aosmith.rbind.io/2020/07/09/ggplot2-override-aes/

# for (bind_type in names(pc_df_ortho)) {
#   
#   plot_df <- pc_df_ortho[[bind_type]]
#   
#   plist <- lapply(1:(num_pcs-1), function(i) {
#     
#     p <-
#       ggplot(plot_df, aes(x = !!sym(paste0("PC", i)),
#                           y = !!sym(paste0("PC", i+1)),
#                           fill = Symbol,
#                           shape = Species)) +
#       geom_point(size = 3, colour = "black") +
#       # coord_fixed() +
#       xlab(paste0("PC", i, "(", pc_ortho[[bind_type]]$Var_explained[i], "%)")) +
#       ylab(paste0("PC", i+1, "(", pc_ortho[[bind_type]]$Var_explained[i+1], "%)")) +
#       ggtitle(paste0("All - ", bind_type)) +
#       theme_classic() +
#       theme(axis.title = element_text(size = 15)) +
#       scale_fill_manual(values = tf_pal) +
#       scale_shape_manual(values = c(21, 25)) +
#       guides(fill = guide_legend(override.aes = list(fill = tf_pal)))
#     
#   })
#   
#   ggsave(
#     filename = paste0(plot_dir, "PCplots_by_TFandSpecies", date, "_Ortho_", bind_type, ".pdf"), 
#     plot = marrangeGrob(plist, nrow = 1, ncol = 1), 
#     width = 7, height = 7
#   )
# }
