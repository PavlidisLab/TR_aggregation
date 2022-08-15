## This script summarizes bind scores across ChIP-seq experiments and fits
## a limma model for TF group differences
## -----------------------------------------------------------------------------

library(tidyverse)
library(cowplot)
library(pheatmap)
library(viridis)
library(grid)
library(gridExtra)
library(limma)
library(edgeR)
library(VennDiagram)
library(egg)
library(GeneOverlap)

date <- "Apr2022"  # most recent data freeze
intensity_flag <- FALSE  # use ouyang scores with or without macs2 intensity
binary_dist <- 25e3  # distance threshold used for binary gene scores
min_peaks <- 100  # how many peaks were required to keep an experiment
peakset <- "idr"  # idr for TF (but overlap still for mecp2) or all overlap
plot_dir <- "~/Plots/Chipseq/Binding_summary/"  # plot export dir
chip_dir <- "~/Data/Annotated_objects/Bind_matrices/Encpipe/" # bind matrix path
outfile <- paste0("~/scratch/R_objects/", date, "_refseq_bind_summary_peakset=", peakset, ".RDS")
chip_type <- "QN_log"  # which chip processing scheme to use

# batch 1 ChIP-seq meta and GEO groups (for blocking in the model)
meta <- read.delim(paste0("~/Data/Metadata/Chipseq/batch1_chip_meta_final_", date, ".tsv"), stringsAsFactors = FALSE)
geo_groups <- read.delim(paste0("~/Data/Metadata/Chipseq/batch1_chip_geo_groups_", date, ".tsv"), stringsAsFactors = FALSE)

# ortho genes
pc_ortho <- read.delim("~/Data/Metadata/hg_mm_1to1_ortho_genes_DIOPT-v8.tsv")

# batch 1 binding matrices
chip_hg <- readRDS(paste0(chip_dir, "Human_refseq_", date, "_processed_bindmat_list_minpeak=", min_peaks, "_ouyang_dc=5000_intensity=", intensity_flag, "_binary=", binary_dist/1e3, "kb_peakset=", peakset, ".RDS"))
chip_mm <- readRDS(paste0(chip_dir, "Mouse_refseq_", date, "_processed_bindmat_list_minpeak=", min_peaks, "_ouyang_dc=5000_intensity=", intensity_flag, "_binary=", binary_dist/1e3, "kb_peakset=", peakset, ".RDS"))
chip_ortho <- readRDS(paste0(chip_dir, "Ortho_refseq_", date, "_processed_bindmat_list_minpeak=", min_peaks, "_ouyang_dc=5000_intensity=", intensity_flag, "_binary=", binary_dist/1e3, "kb_peakset=", peakset, ".RDS"))

stopifnot(all(colnames(chip_ortho$Raw) %in% meta$Experiment_ID))


# Formatting meta for linear model + join GEO groups
# GEO_group_bin lumps experiments from same group but different TR, while
# GEO_group splits same group different TR
#-------------------------------------------------------------------------------


meta_hg <- meta %>%
  left_join(geo_groups[, c("Experiment_ID", "GEO_Group", "GEO_Group_Bin")],
            by = "Experiment_ID") %>% 
  filter(Species == "Human" & Experiment_ID %in% colnames(chip_hg$Raw)) %>%
  arrange(Symbol) %>%
  mutate(
    Has_input = as.logical(Count_input > 0),
    Has_replicates = as.logical(Count_samples > 1),
    Symbol = factor(Symbol, levels = unique(Symbol)),
    N_peaks = log10(N_peaks),
    GEO_group = factor(GEO_Group, levels = unique(GEO_Group)),
    GEO_group_bin = factor(GEO_Group_Bin, levels = unique(GEO_Group_Bin)))
  


meta_mm <- meta %>%
  left_join(geo_groups[, c("Experiment_ID", "GEO_Group", "GEO_Group_Bin")],
            by = "Experiment_ID") %>% 
  filter(Species == "Mouse" & Experiment_ID %in% colnames(chip_mm$Raw)) %>%
  arrange(Symbol) %>%
  mutate(
    Has_input = as.logical(Count_input > 0),
    Has_replicates = as.logical(Count_samples > 1),
    Symbol = factor(Symbol, levels = unique(Symbol)),
    N_peaks = log10(N_peaks),
    GEO_group = factor(GEO_Group, levels = unique(GEO_Group)),
    GEO_group_bin = factor(GEO_Group_Bin, levels = unique(GEO_Group_Bin)))


meta_ortho <- meta %>% 
  filter(Experiment_ID %in% colnames(chip_ortho[[chip_type]])) %>% 
  distinct(Experiment_ID, .keep_all = TRUE) %>% 
  mutate(Symbol = str_to_title(Symbol))


# Gene wise summaries  (mean bind score across experiments (+/- grouped by TR
# status) and the proportion of binary binding assignments)
#-------------------------------------------------------------------------------


bind_df <- function(score_mat, binary_mat, common = NULL) {
  
  # Returns a data frame of the gene-wise mean continuous bind score and 
  # fraction/decimal of experiments bound in binary
  
  if (is.null(common)) {
    common <- intersect(rownames(score_mat), rownames(binary_mat))
  }
  
  df <- data.frame(
    Symbol = common,
    Mean_bind = rowMeans(score_mat[common, , drop = FALSE]),
    Proportion_binary = rowSums(binary_mat[common, , drop = FALSE]) / ncol(binary_mat))
}


tf_bind_list <- function(score_mat, binary_mat, meta) {
  
  # Returns a list of bind dfs for each unique TF in meta
  
  tfs <- unique(meta$Symbol)
  common <- intersect(rownames(score_mat), rownames(binary_mat))
  
  tf_list <- lapply(tfs, function(x) {
    ids <- unique(filter(meta, Symbol == x)$Experiment_ID)
    ids <- intersect(ids, colnames(score_mat))
    bind_df(score_mat = score_mat[common, ids, drop = FALSE], 
            binary_mat = binary_mat[common, ids, drop = FALSE], 
            common = common)
  })
  
  names(tf_list) <- tfs
  return(tf_list)
  
}


top_bound <- function(all_df, tf_list) {
  
  # Return the bind df of each TF for the top bound across all experiments
  
  top_symbol <- filter(all_df, Proportion_binary == max(Proportion_binary))$Symbol
  top_tf <- lapply(tf_list,  function(x) filter(x, Symbol %in% top_symbol))
  return(do.call(rbind, top_tf))
}



# Across all experiments

bind_all <- list(
  Human = bind_df(score_mat = chip_hg[[chip_type]], binary_mat = chip_hg$Binary),
  Mouse = bind_df(score_mat = chip_mm[[chip_type]], binary_mat = chip_mm$Binary),
  Ortho = bind_df(score_mat = chip_ortho[[chip_type]], binary_mat = chip_ortho$Binary)
)


# By TF

bind_tf <- list(
  Human = tf_bind_list(
    score_mat = chip_hg[[chip_type]],
    binary_mat = chip_hg$Binary,
    meta = meta_hg
  ),
  Mouse = tf_bind_list(
    score_mat = chip_mm[[chip_type]],
    binary_mat = chip_mm$Binary,
    meta = meta_mm
  ),
  Ortho = tf_bind_list(
    score_mat = chip_ortho[[chip_type]],
    binary_mat = chip_ortho$Binary,
    meta = meta_ortho
  )
)


# looking at the top bound genes across experiments
# Human - GPAA1 and BCL3 (obvious RUNX1 influence but still intermediate-high
# in others). Mouse - Socs3. Ortho - DMWD

top_bound_list <- list(
  Human = top_bound(all_df = bind_all$Human, tf_list = bind_tf$Human),
  Mouse = top_bound(all_df = bind_all$Mouse, tf_list = bind_tf$Mouse),
  Ortho = top_bound(all_df = bind_all$Ortho, tf_list = bind_tf$Ortho)
)


# Prepare raw bind score matrices for limma voom modeling (which perform quantile
# norm and log).
# NOTE: Currently just examining hist to filter genes by raw counts across
# experiments. Mouse more zero-bound genes compared to human (likely due to
# bad gene annotations).
#-------------------------------------------------------------------------------


# prepare data matrices - min count filter on raw/no-norm score matrix (will be 
# QN+log in voom). for now currently just from examining distn of row/gene sums

sum_hg <- rowSums(chip_hg$Raw)
hist(sum_hg, breaks = 100)
min_hg <- 15
abline(v = min_hg, col = "red")

# Subset matrix to only minimum raw counts across experiments
keep_hg <- sum_hg > min_hg
mat_hg <- chip_hg$Raw[keep_hg, meta_hg$Experiment_ID]
# hist(rowSums(mat_hg), breaks = 100)

# NOTE: for mouse exclude Tcf4 as only has 1 sample

meta_mm <- filter(meta_mm, Symbol != "Tcf4") %>% droplevels()

sum_mm <- rowSums(chip_mm$Raw)
hist(sum_mm, breaks = 100)
min_mm <- 18
abline(v = min_mm, col = "red")

# Subset matrix to only minimum raw counts across experiments
keep_mm <- sum_mm > min_mm
mat_mm <- chip_mm$Raw[keep_mm, meta_mm$Experiment_ID]
hist(rowSums(mat_mm), breaks = 100)


# Design matrices, voom, limma model
# means model (~0) so each TF has a coef corresponding to group mean with no intercept
#-------------------------------------------------------------------------------


design_hg <- model.matrix(
  ~ 0 + Symbol + Has_replicates + Has_input + N_peaks, data = meta_hg)

rownames(design_hg) <- meta_hg$Experiment_ID
colnames(design_hg) <- str_replace(colnames(design_hg), "Symbol", "")


design_mm <- model.matrix(
  ~ 0 + Symbol + Has_replicates + Has_input + N_peaks, data = meta_mm)

rownames(design_mm) <- meta_mm$Experiment_ID
colnames(design_mm) <- str_replace(colnames(design_mm), "Symbol", "")


# voom + fit +/- block on GEO group (requires voom to be run twice)
# https://support.bioconductor.org/p/59700/

# Human
voom_hg <- voom(mat_hg, 
                design = design_hg, 
                normalize.method = "quantile")

# account for increased cor from experiments from the same lab
cor_hg <- duplicateCorrelation(voom_hg, 
                               design = design_hg, 
                               block = meta_hg$GEO_group)

voom_b_hg <- voom(mat_hg, 
                  design = design_hg,
                  normalize.method = "quantile",
                  block = meta_hg$GEO_group,
                  correlation = cor_hg$consensus.correlation)

fit_hg <- lmFit(voom_hg,
                design = design_hg)

fit_b_hg <- lmFit(voom_b_hg, 
                  design = design_hg,
                  block = meta_hg$GEO_group,
                  correlation = cor_hg$consensus.correlation)

# Mouse
voom_mm <- voom(mat_mm, 
                design = design_mm, 
                normalize.method = "quantile")

# account for increased cor from experiments from the same lab
cor_mm <- duplicateCorrelation(voom_mm, 
                               design = design_mm, 
                               block = meta_mm$GEO_group)

voom_b_mm <- voom(mat_mm, 
                  design = design_mm,
                  normalize.method = "quantile",
                  block = meta_mm$GEO_group,
                  correlation = cor_mm$consensus.correlation)

fit_mm <- lmFit(voom_mm,
                design = design_mm)

fit_b_mm <- lmFit(voom_b_mm, 
                  design = design_mm,
                  block = meta_mm$GEO_group,
                  correlation = cor_mm$consensus.correlation)



# interested in each TR vs the rest, so must iteratively construct the 
# appropriate contrast vector. (Law et al., 2020) used as reference
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7873980/ 
#-------------------------------------------------------------------------------


contr_list <- function(meta, contr_vec) {
  # Make a contrast vector for each TF vs the rest: 1 for TF of interest, 
  # -(1/(#TF - 1)) for the rest, leaving nuisance variables as 0
  
  tfs <- unique(meta$Symbol)
  
  clist <- lapply(tfs, function(x) {
    contr_vec[x] <- 1
    contr_vec[names(contr_vec) %in% setdiff(tfs, x)] <- -(1/(length(tfs)-1))
    return(contr_vec)
  })
  names(clist) <- tfs
  return(clist)
}


contr_hg <- rep(0, length(colnames(coef(fit_hg))))
names(contr_hg) <- colnames(coef(fit_hg))

contr_mm <- rep(0, length(colnames(coef(fit_mm))))
names(contr_mm) <- colnames(coef(fit_mm))


clist_hg <- contr_list(meta_hg, contr_hg)
clist_mm <- contr_list(meta_mm, contr_mm)


# For each symbol, get the top results for the symbol vs all contrast
#-------------------------------------------------------------------------------


top_fit <- function(contr_list, fit) {
  # Apply eBayes + contrast fit to the model and extract toptable for each contr
  lapply(contr_list, function(x) {
    contr_fit <- eBayes(contrasts.fit(fit, x))
    topTable(contr_fit, n = Inf)
  })
}


top_hg <- top_fit(clist_hg, fit_hg)  # -block by GEO group
top_b_hg <- top_fit(clist_hg, fit_b_hg) # +block by GEO group

top_mm <- top_fit(clist_mm, fit_mm)  # -block by GEO group
top_b_mm <- top_fit(clist_mm, fit_b_mm) # +block by GEO group


# count of pos FC sig genes at FDR05
n_hg <- sapply(top_hg, function(x) nrow(filter(x, logFC > 0 & adj.P.Val < 0.05)))
n_mm <- sapply(top_mm, function(x) nrow(filter(x, logFC > 0 & adj.P.Val < 0.05)))
n_b_hg <- sapply(top_b_hg, function(x) nrow(filter(x, logFC > 0 & adj.P.Val < 0.05)))
n_b_mm <- sapply(top_b_mm, function(x) nrow(filter(x, logFC > 0 & adj.P.Val < 0.05)))
summary(c(n_b_hg, n_b_mm))


# Significant shared ortho genes
#-------------------------------------------------------------------------------


tfs_ortho <- unique(meta_mm$Symbol)


top_ortho <- lapply(tfs_ortho, function(x) {
  
  hg <- filter(top_b_hg[[str_to_upper(x)]], logFC > 0 & adj.P.Val < 0.05) %>% 
    rownames_to_column(var = "Symbol_hg") %>% 
    filter(Symbol_hg %in% pc_ortho$Symbol_hg) %>% 
    left_join(pc_ortho[, c("Symbol_hg", "ID")], by = "Symbol_hg")
  
  mm <- filter(top_b_mm[[x]], logFC > 0 & adj.P.Val < 0.05) %>% 
    rownames_to_column(var = "Symbol_mm") %>% 
    filter(Symbol_mm %in% pc_ortho$Symbol_mm) %>% 
    left_join(pc_ortho[, c("Symbol_mm", "ID")], by = "Symbol_mm")
  
  ortho <- hg %>% 
    filter(ID %in% mm$ID) %>% 
    left_join(mm, by = "ID", suffix = c("_hg", "_mm"))
  
  list(Human = hg, Mouse = mm, Ortho = ortho)
  
})
names(top_ortho) <- tfs_ortho


n_ortho <- unlist(lapply(top_ortho, function(x) nrow(x$Ortho)))

# FET for overlapping significant genes

ortho_overlap <- lapply(top_ortho, function(x) {
  a <- x$Human$ID
  b <- x$Mouse$ID
  c <- testGeneOverlap(newGeneOverlap(a, b, genome.size = nrow(pc_ortho)))
})


ortho_pvals <- lapply(ortho_overlap, function(x) x@pval)


# Look at +/- blocking. Specifically, rank cor of adjusted pvals
#-------------------------------------------------------------------------------


# Human
blist_hg <- lapply(unique(meta_hg$Symbol), function(x) {
  
  noblock <- top_hg[[x]] %>% 
   dplyr::select(logFC, t, adj.P.Val) %>% 
    setNames(paste0('Noblock_', names(.))) %>% 
    rownames_to_column(var = "Symbol")

  block <- top_b_hg[[x]] %>% 
    dplyr::select(logFC, t, adj.P.Val) %>% 
    setNames(paste0('Block_', names(.))) %>% 
    rownames_to_column(var = "Symbol")
  
  left_join(noblock, block, by = "Symbol")

})
names(blist_hg) <- unique(meta_hg$Symbol)


# Mouse
blist_mm <- lapply(unique(meta_mm$Symbol), function(x) {
  
  noblock <- top_mm[[x]] %>% 
    dplyr::select(logFC, t, adj.P.Val) %>% 
    setNames(paste0('Noblock_', names(.))) %>% 
    rownames_to_column(var = "Symbol")
  
  block <- top_b_mm[[x]] %>% 
    dplyr::select(logFC, t, adj.P.Val) %>% 
    setNames(paste0('Block_', names(.))) %>% 
    rownames_to_column(var = "Symbol")
  
  left_join(noblock, block, by = "Symbol")
  
})
names(blist_mm) <- unique(meta_mm$Symbol)


scor_hg <- sapply(blist_hg, function(x) {
  cor(x[, "Block_adj.P.Val",], x[, "Noblock_adj.P.Val"], method = "spearman")
})

scor_mm <- sapply(blist_mm, function(x) {
  cor(x[, "Block_adj.P.Val",], x[, "Noblock_adj.P.Val"], method = "spearman")
})


# Explore genes associated with nuisance variables
#-------------------------------------------------------------------------------

# Human
top_peaks_hg <- topTable(eBayes(fit_hg), n = Inf, coef = "N_peaks")
top_input_hg <- topTable(eBayes(fit_hg), n = Inf, coef = "Has_inputTRUE")
top_reps_hg <- topTable(eBayes(fit_hg), n = Inf, coef = "Has_replicatesTRUE")

plot(x = meta_hg$N_peak, y = log2(mat_hg["PRAMEF18",]+1))
boxplot(log2(mat_hg["BTRC",]+1) ~ meta_hg$Has_input)
boxplot(log2(mat_hg["TNS3",]+1) ~ meta_hg$Has_replicates)

# Mouse
top_peaks_mm <- topTable(eBayes(fit_mm), n = Inf, coef = "N_peaks")
top_input_mm <- topTable(eBayes(fit_mm), n = Inf, coef = "Has_inputTRUE")
top_reps_mm <- topTable(eBayes(fit_mm), n = Inf, coef = "Has_replicatesTRUE")

plot(x = meta_mm$N_peaks, y = log2(mat_mm["Kif23",]+1))
boxplot(log2(mat_mm["Sptan1",]+1) ~ meta_mm$Has_input)
boxplot(log2(mat_mm["Hipk2",]+1) ~ meta_mm$Has_replicates)


# Top diff bound genes dataframe
# For ASCL1 human, find KRTAP[4/9]-* gene cluster as top results and MECP2 has
# PCDHG*. Anno influenced by distance, so this behavior is not unexpected. Just 
# take most sig for plotting
#-------------------------------------------------------------------------------


top_n <- 8


top_cluster <- function(df) {
  
  # Return a df where only the most significant clustered genes are kept
  
  cluster_genes <- c("^KRTAP.*$", "^PCDHG.*$")
  
  top_cluster <- lapply(cluster_genes, function(x) {
    df %>% 
    rownames_to_column(var = "Symbol") %>%
    filter(str_detect(Symbol, x)) %>%
    arrange(adj.P.Val) %>%
    slice_head(n = 1)
  })
  
  top_cluster <- do.call(rbind, top_cluster)
  
  out <- df %>% 
    rownames_to_column(var = "Symbol") %>%
    filter(!str_detect(Symbol, paste(cluster_genes, collapse = "|"))) %>%
    rbind(top_cluster) %>% 
    arrange(adj.P.Val)

  rownames(out) <- out$Symbol
  dplyr::select(out, -Symbol)
}


# Human
genes_hg <- unlist(lapply(top_b_hg, function(x) {
  x %>%
    top_cluster() %>% 
    filter(logFC > 0) %>% 
    arrange(adj.P.Val) %>% 
    slice_head(n = top_n) %>% 
    rownames()
}), use.names = FALSE)

fc_hg <- lapply(top_b_hg, function(x) x[genes_hg, "logFC", drop = FALSE])
topdf_hg <- do.call(cbind, fc_hg)
colnames(topdf_hg) <- unique(meta_hg$Symbol)


# Mouse
genes_mm <- unlist(lapply(top_b_mm, function(x) {
  x %>% 
    filter(logFC > 0) %>% 
    arrange(adj.P.Val) %>% 
    slice_head(n = top_n) %>% 
    rownames()
}), use.names = FALSE)
fc_mm <- lapply(top_b_mm, function(x) x[genes_mm, "logFC", drop = FALSE])
topdf_mm <- do.call(cbind, fc_mm)
colnames(topdf_mm) <- unique(meta_mm$Symbol)



# Save out objects

saveRDS(
  list(
    Human_fit = top_hg,
    Human_block_fit = top_b_hg,
    Human_bind_all = bind_all$Human,
    Human_bind_tf = bind_tf$Human,
    Mouse_fit = top_mm,
    Mouse_block_fit = top_b_mm,
    Mouse_bind_all = bind_all$Mouse,
    Mouse_bind_tf = bind_tf$Mouse,
    Ortho_bind_all = bind_all$Ortho,
    Ortho_bind_tf = bind_tf$Ortho
  ),
  file = outfile
)


# Plots
#-------------------------------------------------------------------------------


# heatmap of top binding scores per TF

pal_length <- 11
heatmap_pal <- viridis::inferno(pal_length)
color_breaks_hg <- seq(min(topdf_hg), max(topdf_hg), length.out = pal_length)
color_breaks_mm <- seq(min(topdf_mm), max(topdf_mm), length.out = pal_length)


pheatmap(
  topdf_hg,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  cellwidth = 40,
  border_color = NA,
  color = heatmap_pal,
  gaps_row = seq(0, top_n*n_distinct(meta_hg$Symbol), top_n),
  gaps_col = 1:n_distinct(meta_hg$Symbol),
  height = 12,
  width = 10,
  angle_col = 90,
  fontsize_row = 13,
  fontsize_col = 15,
  filename = paste0(plot_dir, "Heatmap_Human_", date, "_minpeak=", min_peaks, "_ouyang_dc=5000_intensity=", intensity_flag, "_peakset=", peakset, ".png")
)


pheatmap(
  topdf_mm,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  cellwidth = 40,
  border_color = NA,
  color = heatmap_pal,
  gaps_row = seq(0, top_n*n_distinct(meta_mm$Symbol), top_n),
  gaps_col = 1:n_distinct(meta_mm$Symbol),
  height = 12,
  width = 10,
  angle_col = 90,
  fontsize_row = 13,
  fontsize_col = 15,
  filename = paste0(plot_dir, "Heatmap_Mouse_", date, "_minpeak=", min_peaks, "_ouyang_dc=5000_intensity=", intensity_flag, "_peakset=", peakset, ".png")
)


# Violin+Box plot of bind scores across TFs


vbplot <- function(mat, plot_genes, ncol = 1, species) {
  
  # Long df for plotting
  all_df <- data.frame(mat[plot_genes, ]) %>% 
    rownames_to_column(var = "Symbol") %>% 
    reshape2::melt(id = "Symbol", var = "Experiment_ID")
  
  if (species == "Human") {
    all_df$TF <- str_to_upper(str_split(all_df$Experiment_ID, "_", simplify = TRUE)[, 2])
  } else if (species == "Mouse") {
    all_df$TF <- str_split(all_df$Experiment_ID, "_", simplify = TRUE)[, 2]
  }

  # List of plots (one gene per element) to be combined after
  plot_list <- lapply(plot_genes, function(x) {
    
    gene_df <- filter(all_df, Symbol == x)
    
    # violin plot won't work for n < 3; use points instead
    which_lt3 <- names(which(table(gene_df$TF) < 3))
    
    # Build plot
    p <-
      ggplot(gene_df, aes(x = TF, y = value)) +
      geom_violin(data = gene_df[!(gene_df$TF %in% which_lt3),], 
                  width = 0.6,
                  fill = "lightslategrey") +
      geom_boxplot(data = gene_df[!(gene_df$TF %in% which_lt3),], 
                   width = 0.1,
                   fill = "white") +
      geom_point(data = gene_df[gene_df$TF %in% which_lt3,],
                 shape = 21,
                 colour = "black",
                 fill = "lightslategrey") +
      theme_classic() +
      ggtitle(x) +
      ylab("Binding score") +
      theme(
        axis.text = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.title.x = element_blank(),
        plot.title = element_text(size = 20, hjust = 0.5, colour = "darkred")
      )
    
    # Only want TF x axis labels on bottom plot
    if (x == plot_genes[length(plot_genes)]) {
      p <- p + theme(axis.text.x = element_text(angle = 90, vjust = 0.45))
    } else {
      p <- p + theme(axis.text.x = element_blank())
    }
    return(p)
  })
  
  # egg:: preserves consistent plot panel sizing, unlike cowplot
  return(egg::ggarrange(plots = plot_list, ncol = ncol))
  
}



# Plot genes as most sig from blocked model
pgenes_hg <- unlist(lapply(top_b_hg, function(x) {
  x %>%
    top_cluster() %>% 
    filter(logFC > 0) %>% 
    arrange(adj.P.Val) %>% 
    slice_head(n = 1) %>% 
    rownames()
}), use.names = FALSE)


p1 <- vbplot(mat = chip_hg$QN_log, plot_genes = pgenes_hg, species = "Human")
p1 <- vbplot(mat = mat_hg, plot_genes = pgenes_hg, species = "Human")

ggsave(p1, dpi = 300, height = 20, width = 8,
       filename = paste0(plot_dir, "Bplot_bindscore_Human_", date, "_minpeak=", min_peaks, "_QNL_ouyang_dc=5000_intensity=", intensity_flag, "_peakset=", peakset, ".png"))


# Plot genes as most sig from blocked model
pgenes_mm <- unlist(lapply(top_b_mm, function(x) {
  x %>%
    top_cluster() %>% 
    filter(logFC > 0) %>% 
    arrange(adj.P.Val) %>% 
    slice_head(n = 1) %>% 
    rownames()
}), use.names = FALSE)


p2 <- vbplot(mat = chip_mm$QN_log, plot_genes = pgenes_mm, species = "Mouse")

ggsave(p2, dpi = 300, height = 20, width = 8,
       filename = paste0(plot_dir, "Bplot_bindscore_Mouse_", date, "_minpeak=", min_peaks, "_QNL_ouyang_dc=5000_intensity=", intensity_flag, "_peakset=", peakset, ".png"))


# Looking at the top bound genes across all experiments


top_bind_hg <- c("GPAA1", "GSE1")

plot_df_hg <- bind_all$Human %>% 
  mutate(
    Group = Symbol %in% top_bind_hg,
    Top = as.character(Symbol) %in% c(
      as.character(slice_max(bind_all$Human, Mean_bind, n = 50)$Symbol),
      as.character(slice_max(bind_all$Human, Proportion_binary, n = 50)$Symbol)
    )) %>% 
  arrange(desc(Proportion_binary))


p3a <- 
  ggplot() +
  geom_point(data = plot_df_hg, 
             aes(y = Proportion_binary, x = Mean_bind),
             shape = 21, size = 2.3, alpha = 0.4, colour = "black") +
  geom_text(data = filter(plot_df_hg, Top),
            check_overlap = TRUE,
            aes(y = Proportion_binary, x = Mean_bind, label = Symbol)) +
  ylab("Proportion of bound experiments (+/- 25kb)") +
  xlab("Mean binding score") +
  ggtitle(paste0("Human n=", nrow(meta_hg), " experiments")) +
  theme_classic() +
  theme(
    axis.text = element_text(size = 20),
    axis.title = element_text(size = 20),
    plot.title = element_text(size = 20, hjust = 0.5),
    legend.position = "none"
  )

p3b <- vbplot(mat = chip_hg$QN_log, plot_genes = top_bind_hg, species = "Human")

ggsave(plot_grid(p3a, p3b, nrow = 1),
       dpi = 300, device = "png", width = 16, height = 8,
       filename = paste0(plot_dir, "Human_top_bound.png"))


top_bind_mm <- c("Socs3", "Gse1")


plot_df_mm <- bind_all$Mouse %>% 
  mutate(
    Group = Symbol %in% top_bind_mm,
    Top = as.character(Symbol) %in% c(
      as.character(slice_max(bind_all$Mouse, Mean_bind, n = 50)$Symbol),
      as.character(slice_max(bind_all$Mouse, Proportion_binary, n = 50)$Symbol)
    )) %>% 
  arrange(desc(Mean_bind))

p4a <- 
  ggplot() +
  geom_point(data = plot_df_mm, 
             aes(y = Proportion_binary, x = Mean_bind),
             shape = 21, size = 2.3, alpha = 0.4, colour = "black") +
  geom_text(data = filter(plot_df_mm, Top),
            check_overlap = TRUE,
            aes(y = Proportion_binary, x = Mean_bind, label = Symbol)) +
  ylab("Proportion of bound experiments (+/- 25kb)") +
  xlab("Mean binding score") +
  ggtitle(paste0("Mouse n=", nrow(meta_mm), " experiments")) +
  theme_classic() +
  theme(
    axis.text = element_text(size = 20),
    axis.title = element_text(size = 20),
    plot.title = element_text(size = 20, hjust = 0.5),
    legend.position = "none"
  )

p4b <- vbplot(mat = chip_mm$QN_log, plot_genes = top_bind_mm, species = "Mouse")

ggsave(plot_grid(p4a, p4b, nrow = 1),
       dpi = 300, device = "png", width = 16, height = 8,
       filename = paste0(plot_dir, "Mouse_top_bound.png"))


# Venn diagram of ortho overlap


# set_cols <- c("royalblue", "goldenrod")
# 
# for (tf in names(top_ortho)) {
#   
#   hg <- ortho_overlap[[tf]]@listA
#   mm <- ortho_overlap[[tf]]@listB
#   
#   venn.diagram(x = list(hg, mm),
#                category.names = c("Human", "Mouse"),
#                filename = paste0(plot_dir, date, "_venn_ortho_overlap_", tf, ".png"),
#                imagetype= "png",
#                height = 1024, 
#                width = 1024, 
#                resolution = 300,
#                lwd = 2,
#                lty = 'blank',
#                fill = set_cols,
#                main = str_to_upper(tf),
#                main.pos = c(0.5, 0.9),
#                cat.default.pos = "outer",
#                cat.pos = c(-27, 27),
#                ext.text = FALSE)
# }
