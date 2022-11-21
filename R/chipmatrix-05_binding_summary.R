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
library(GeneOverlap)
source("R/setup-01_config.R")
source("R/utils/plot_functions.R")

plot_dir <- paste0(cplot_dir, "Binding_summary/")
outfile <- paste0(scratch_dir, date, "_refseq_bind_summary.RDS")

# Loading ChIP-seq data
chip_type <- "QN_log"  # which chip processing scheme to use
meta <- read.delim(paste0(meta_dir, "Chipseq/batch1_chip_meta_final_", date, ".tsv"), stringsAsFactors = FALSE)
chip_hg <- readRDS(paste0(cmat_dir, "Human_refseq_", date, "_processed_bindmat_list_minpeak=", min_peaks, "_ouyang_dc=5000_intensity=FALSE_binary=", binary_dist/1e3, "kb.RDS"))
chip_mm <- readRDS(paste0(cmat_dir, "Mouse_refseq_", date, "_processed_bindmat_list_minpeak=", min_peaks, "_ouyang_dc=5000_intensity=FALSE_binary=", binary_dist/1e3, "kb.RDS"))
chip_ortho <- readRDS(paste0(cmat_dir, "Ortho_refseq_", date, "_processed_bindmat_list_minpeak=", min_peaks, "_ouyang_dc=5000_intensity=FALSE_binary=", binary_dist/1e3, "kb.RDS"))
stopifnot(all(colnames(chip_ortho$Raw) %in% meta$Experiment_ID))

# Curated GEO groups for blocking variable in the model
geo_groups <- read.delim(paste0(meta_dir, "Chipseq/batch1_chip_geo_groups_", date, ".tsv"), stringsAsFactors = FALSE)

# Ortho genes
pc_ortho <- read.delim(paste0(meta_dir, "hg_mm_1to1_ortho_genes_DIOPT-v8.tsv"), stringsAsFactors = FALSE)


# Formatting meta for linear model + join GEO groups
# GEO_group_bin lumps experiments from same group but different TR
# GEO_group splits same group testing different TRs
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


# Gene wise summaries: mean bind score across experiments (+/- grouped by TR
# status) and the proportion of binary binding assignments
#-------------------------------------------------------------------------------


# Dataframe of binding summaries across all experiment in provided matrices

bind_df <- function(score_mat, binary_mat) {
  
  common <- intersect(rownames(score_mat), rownames(binary_mat))
  
  df <- data.frame(
    Symbol = common,
    Mean_bind = rowMeans(score_mat[common, , drop = FALSE]),
    Proportion_binary = rowSums(binary_mat[common, , drop = FALSE]) / ncol(binary_mat))
  
  return(df)
}


# Returns a list of bind dfs for each unique TF in meta


tf_bind_list <- function(score_mat, binary_mat, meta) {
  
  tfs <- unique(meta$Symbol)
  common <- intersect(rownames(score_mat), rownames(binary_mat))
  
  tf_list <- lapply(tfs, function(x) {
    
    ids <- unique(filter(meta, Symbol == x)$Experiment_ID)
    ids <- intersect(ids, colnames(score_mat))
    bind_df(score_mat = score_mat[common, ids, drop = FALSE], 
            binary_mat = binary_mat[common, ids, drop = FALSE])
  })
  
  names(tf_list) <- tfs
  
  return(tf_list)
}


# Return the bind df of each TF for the top bound gene across all experiments

top_bound <- function(all_df, tf_list) {
  
  top_symbol <- filter(all_df, Proportion_binary == max(Proportion_binary))$Symbol
  top_tf <- lapply(tf_list,  function(x) filter(x, Symbol %in% top_symbol))
  
  return(do.call(rbind, top_tf))
}



# Binding summary across all experiments - genes that are commonly bound

bind_all <- list(
  Human = bind_df(score_mat = chip_hg[[chip_type]], binary_mat = chip_hg$Binary),
  Mouse = bind_df(score_mat = chip_mm[[chip_type]], binary_mat = chip_mm$Binary),
  Ortho = bind_df(score_mat = chip_ortho[[chip_type]], binary_mat = chip_ortho$Binary)
)


# Binding summary by TF

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


# Prepare raw bind score matrices for limma voom modeling, since it performs 
# quantile norm and log. Examining hist to filter genes by raw counts across
# experiments. Mouse more zero-bound genes compared to human (likely due to
# bad gene annotations).
#-------------------------------------------------------------------------------


sum_hg <- rowSums(chip_hg$Raw)
hist(sum_hg, breaks = 100)
min_hg <- 15
abline(v = min_hg, col = "red")

# Subset matrix to only minimum raw counts across experiments
keep_hg <- sum_hg > min_hg
mat_hg <- chip_hg$Raw[keep_hg, meta_hg$Experiment_ID]
hist(rowSums(mat_hg), breaks = 100)

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
# ~ 0 + Symbol + Has_replicates + Has_input + N_peaks
# Means model (~0) so each TF has a coef corresponding to group mean with no intercept
#-------------------------------------------------------------------------------


design_mat <- function(meta) {
  
  mat <- model.matrix(
    ~ 0 + Symbol + Has_replicates + Has_input + N_peaks, data = meta)
  
  rownames(mat) <- meta$Experiment_ID
  colnames(mat) <- str_replace(colnames(mat), "Symbol", "")
  
  return(mat)
}


design_hg <- design_mat(meta_hg)
design_mm <- design_mat(meta_mm)


# voom data and fit model +/- blocking on GEO group

# voom is a transformation for count matrices with a mean-variance relationship,
# as seen in the binding matrices. While the bind matrices are not integer
# counts like in RNA-seq, G. Smyth says it is appropriate to use on numeric
# matrices https://support.bioconductor.org/p/45695/. For blocking it is also
# recommended to run voom twice https://support.bioconductor.org/p/59700/

# +block: considers the enhanced correlation between samples from the same lab 
# (inferred from GEO group) by fitting a mixed effect linear model via REML.
#-------------------------------------------------------------------------------


noblock_fit <- function(dat_mat, design_mat) {
  
  v <- voom(dat_mat, design = design_mat, normalize.method = "quantile")
  fit <- lmFit(v, design = design_mat)
  
  return(fit)
}


block_fit <- function(dat_mat, design_mat, meta) {
  
  v <- voom(dat_mat, design = design_mat, normalize.method = "quantile")
  
  dc <- duplicateCorrelation(v, design = design_mat, block = meta$GEO_group)
  
  v_b <- voom(dat_mat, 
              design = design_mat,
              normalize.method = "quantile",
              block = meta$GEO_group,
              correlation = dc$consensus.correlation)
  
  fit <- lmFit(v_b, 
               design = design_mat,
               block = meta$GEO_group,
               correlation = dc$consensus.correlation)
  
  return(fit)
}



fit_l <- list(
  Human_block = block_fit(mat_hg, design_hg, meta_hg),
  Human = noblock_fit(mat_hg, design_hg),
  Mouse_block = block_fit(mat_mm, design_mm, meta_mm),
  Mouse = noblock_fit(mat_mm, design_mm)
)


# Interested in each TR vs the rest, so must construct the appropriate contrast
# vectors. This is setting 1 for the TR of interest, and -1/(nTR-1) for the 
# other nTR-1 TRs, averaging their contributions as the "out" group. Nuisance
# variables are kept as 0 (the output still corrects for them)
# Law et al., 2020 used as reference 
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7873980/ 
#-------------------------------------------------------------------------------


contr_list <- function(meta, fit) {
  
  tfs <- unique(meta$Symbol)
  
  contr_vec <- rep(0, length(colnames(coef(fit))))
  names(contr_vec) <- colnames(coef(fit))
  
  clist <- lapply(tfs, function(x) {
    contr_vec[x] <- 1
    contr_vec[names(contr_vec) %in% setdiff(tfs, x)] <- -(1/(length(tfs) - 1))
    return(contr_vec)
  })
  names(clist) <- tfs
  
  return(clist)
}


clist_hg <- contr_list(meta_hg, fit_l$Human_block)
clist_mm <- contr_list(meta_mm, fit_l$Mouse_block)


# For each symbol, get the model estimates for the symbol vs all contrast
#-------------------------------------------------------------------------------

# Apply eBayes + contrast fit to the model and extract toptable for each contr

top_fit <- function(contr_list, fit) {
  lapply(contr_list, function(x) {
    contr_fit <- eBayes(contrasts.fit(fit, x))
    topTable(contr_fit, adjust.method = "BH", number = Inf)
  })
}


count_sig <- function(top_l, fdr = 0.05, lfc = 0) {
  vapply(top_l, function(x) {
    nrow(filter(x, logFC > lfc & adj.P.Val < fdr))
  }, FUN.VALUE = numeric(1))
}


top_l <- list(
  Human_block = top_fit(clist_hg, fit_l$Human_block),
  Human = top_fit(clist_hg, fit_l$Human),
  Mouse_block = top_fit(clist_mm, fit_l$Mouse_block),
  Mouse = top_fit(clist_mm, fit_l$Mouse)
)


# Count of positive FC sig genes at FDR05
n_l <- lapply(top_l, count_sig)


# Paper reports blocked model counts combining both species
summary(c(n_l$Human_block, n_l$Mouse_block))


# Significant shared ortho genes
#-------------------------------------------------------------------------------


shared_sig <- function(tfs, pc_ortho, top_hg, top_mm) {
  
  top_ortho <- lapply(tfs, function(x) {
    
    hg <- top_hg[[str_to_upper(x)]] %>% 
      filter(logFC > 0 & adj.P.Val < 0.05) %>% 
      rownames_to_column(var = "Symbol_hg") %>% 
      filter(Symbol_hg %in% pc_ortho$Symbol_hg) %>% 
      left_join(pc_ortho[, c("Symbol_hg", "ID")], by = "Symbol_hg")
    
    mm <- top_mm[[str_to_title(x)]] %>% 
      filter(logFC > 0 & adj.P.Val < 0.05) %>% 
      rownames_to_column(var = "Symbol_mm") %>% 
      filter(Symbol_mm %in% pc_ortho$Symbol_mm) %>% 
      left_join(pc_ortho[, c("Symbol_mm", "ID")], by = "Symbol_mm")
    
    ortho <- hg %>% 
      filter(ID %in% mm$ID) %>% 
      left_join(mm, by = "ID", suffix = c("_hg", "_mm"))
    
    list(Human = hg, Mouse = mm, Ortho = ortho)
  })
  names(top_ortho) <- tfs
 
  return(top_ortho) 
}
  

# mouse b/c excluding TCF4
tfs_ortho <- unique(meta_mm$Symbol)
top_ortho <- shared_sig(tfs_ortho, pc_ortho, top_l$Human_block, top_l$Mouse_block)

# Count shared sig
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


join_top <- function(top_b, top_nob, meta) {
  
  tfs <- unique(meta$Symbol)
  
  l <- lapply(unique(tfs), function(x) {
    
    block <- top_b[[x]] %>% 
      dplyr::select(logFC, t, adj.P.Val) %>% 
      setNames(paste0('Noblock_', names(.))) %>% 
      rownames_to_column(var = "Symbol")
    
    noblock <- top_nob[[x]] %>% 
      dplyr::select(logFC, t, adj.P.Val) %>% 
      setNames(paste0('Block_', names(.))) %>% 
      rownames_to_column(var = "Symbol")
    
    left_join(noblock, block, by = "Symbol")
    
  })
  names(l) <- unique(tfs)
  
  return(l)
}


block_scor <- function(join_l) {
  vapply(join_l, function(x) {
    cor(x[, "Block_adj.P.Val",], x[, "Noblock_adj.P.Val"], method = "spearman")
  }, FUN.VALUE = numeric(1))
}



# Human: See that ASCL1 has lowest cor (~0.93), makes sense as this TF had the 
# most data sets from same group
join_l_hg <- join_top(top_l$Human_block, top_l$Human, meta_hg)
scor_hg <- block_scor(join_l_hg)

# Mouse: Hes lowest cor (~0.22), makes sense as this only has 2 data sets from
# the same group (+/- LPS)
join_l_mm <- join_top(top_l$Mouse_block, top_l$Mouse, meta_mm)
scor_mm <- block_scor(join_l_mm)


# Explore genes associated with nuisance variables
#-------------------------------------------------------------------------------

# Human
top_peaks_hg <- topTable(eBayes(fit_l$Human), number = Inf, coef = "N_peaks")
top_input_hg <- topTable(eBayes(fit_l$Human), number = Inf, coef = "Has_inputTRUE")
top_reps_hg <- topTable(eBayes(fit_l$Human), number = Inf, coef = "Has_replicatesTRUE")

plot(x = meta_hg$N_peak, y = chip_hg$QN_log[rownames(top_peaks_hg)[1], ])
boxplot(chip_hg$QN_log[rownames(top_input_hg)[1], ] ~ meta_hg$Has_input)
boxplot(chip_hg$QN_log[rownames(top_reps_hg)[1], ] ~ meta_hg$Has_replicates)

# Mouse
top_peaks_mm <- topTable(eBayes(fit_l$Mouse), number = Inf, coef = "N_peaks")
top_input_mm <- topTable(eBayes(fit_l$Mouse), number = Inf, coef = "Has_inputTRUE")
top_reps_mm <- topTable(eBayes(fit_l$Mouse), number = Inf, coef = "Has_replicatesTRUE")

plot(x = meta_mm$N_peaks, y = chip_mm$QN_log[rownames(top_peaks_mm)[1], meta_mm$Experiment_ID])
boxplot(chip_mm$QN_log[rownames(top_input_mm)[1], meta_mm$Experiment_ID] ~ meta_mm$Has_input)
boxplot(chip_mm$QN_log[rownames(top_reps_mm)[1], meta_mm$Experiment_ID] ~ meta_mm$Has_replicates)


# Top diff bound genes dataframe
# Get a dataframe for plotting that has the FC of each TRs n most significant
# positive genes. Get these FCs for each TR to illustrate within-TR enrichment
# For ASCL1 human, find KRTAP[4/9]-* gene cluster as top results and MECP2 has
# PCDHG*. Anno influenced by distance, so this behavior is not unexpected. Just 
# take most sig gene of these clusters for plotting
#-------------------------------------------------------------------------------


# Return a df where only the most significant clustered genes are kept

top_cluster <- function(df) {
  
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


# Return a df of the FC of the most significant positive genes for each TR,
# as organized in the list of top tables

top_n <- 8

top_genes <- function(top_l, top_n) {
  
  genes <- unlist(lapply(top_l, function(x) {
    x %>%
      top_cluster() %>% 
      filter(logFC > 0) %>% 
      arrange(adj.P.Val) %>% 
      slice_head(n = top_n) %>% 
      rownames()
  }), use.names = FALSE)
  
  fc <- lapply(top_l, function(x) x[genes, "logFC", drop = FALSE])
  top_df <- do.call(cbind, fc)
  colnames(top_df) <- names(top_l)
  
  return(top_df)
}


genes_hg <- top_genes(top_l$Human_block, top_n)
genes_mm <- top_genes(top_l$Mouse_block, top_n)


# Save out objects


saveRDS(
  list(
    Human_fit = top_l$Human,
    Human_block_fit = top_l$Human_block,
    Human_bind_all = bind_all$Human,
    Human_bind_tf = bind_tf$Human,
    Mouse_fit = top_l$Mouse,
    Mouse_block_fit = top_l$Mouse_block,
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


color_breaks_hg <- seq(min(genes_hg), max(genes_hg), length.out = pal_length)
color_breaks_mm <- seq(min(genes_mm), max(genes_mm), length.out = pal_length)


pheatmap(
  genes_hg,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  cellwidth = 40,
  border_color = NA,
  color = heatmap_pal,
  gaps_row = seq(0, top_n * n_distinct(meta_hg$Symbol), top_n),
  gaps_col = 1:n_distinct(meta_hg$Symbol),
  height = 12,
  width = 10,
  angle_col = 90,
  fontsize_row = 13,
  fontsize_col = 15,
  filename = paste0(plot_dir, "Bind_specificity_heatmap_human_", date, ".png")
)



pheatmap(
  genes_mm,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  cellwidth = 40,
  border_color = NA,
  color = heatmap_pal,
  gaps_row = seq(0, top_n * n_distinct(meta_mm$Symbol), top_n),
  gaps_col = 1:n_distinct(meta_mm$Symbol),
  height = 12,
  width = 10,
  angle_col = 90,
  fontsize_row = 13,
  fontsize_col = 15,
  filename = paste0(plot_dir, "Bind_specificity_heatmap_mouse_", date, ".png")
)



# Plot genes as most sig from blocked model. Note that this is using the
# QNL matrix (used for rest of analysis and mean summaries), not the voom
# transformed matrix. Find that these are equivalent and want to keep same
# scale as QNL as used for rest of analysis


top_gene_hg <- rownames(top_genes(top_l$Human_block, 1))

p1 <- bind_vboxplot(mat = chip_hg$QN_log, plot_genes = top_gene_hg, species = "Human")

ggsave(p1, dpi = 300, height = 20, width = 8,
      filename = paste0(plot_dir, "VBplot_bindscore_human_", date, ".png"))

# Hacky: because Hes6 top for Ascl1 and Neurod1, must remove suffix that was
# coerced due to duplicate rownames generated by top_genes()
top_gene_mm <- rownames(top_genes(top_l$Mouse_block, 1))
top_gene_mm <- str_replace(top_gene_mm, "\\.1", "")

p2 <- bind_vboxplot(mat = chip_mm$QN_log, plot_genes = top_gene_mm, species = "Mouse")

ggsave(p2, dpi = 300, height = 20, width = 8,
       filename = paste0(plot_dir, "VBplot_bindscore_mouse_", date, ".png"))


# Looking at the top bound genes across all experiments. 


top_bind_hg <- c("GPAA1", "GSE1")


p3 <- topbound_plot(bind_df = bind_all$Human, 
                    topbound_symbols = top_bind_hg, 
                    bind_mat = chip_hg$QN_log,
                    title = paste0("Human n=", nrow(meta_hg), " experiments"),
                    species = "Human")


ggsave(p3, dpi = 300, device = "png", width = 16, height = 8,
       filename = paste0(plot_dir, "Human_top_bound.png"))


top_bind_mm <- c("Socs3", "Gse1")


p4 <- topbound_plot(bind_df = bind_all$Mouse, 
                    topbound_symbols = top_bind_mm, 
                    bind_mat = chip_mm$QN_log,
                    title = paste0("Mouse n=", nrow(meta_mm), " experiments"),
                    species = "Mouse")


ggsave(p4, dpi = 300, device = "png", width = 16, height = 8,
       filename = paste0(plot_dir, "Human_top_bound.png"))
