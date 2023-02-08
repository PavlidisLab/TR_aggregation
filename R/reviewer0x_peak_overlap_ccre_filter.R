## Comparing overlap of ChIP-seq experiments, with and without filtering peaks
## for cCRE overlap.
## -----------------------------------------------------------------------------

library(tidyverse)
library(GenomicRanges)
library(jaccard)
library(dendextend)
source("R/setup-01_config.R")
source("R/utils/range_table_functions.R")
source("R/utils/similarity_functions.R")
source("R/utils/plot_functions.R")

pipeline_dir <- paste0(pipeout_dir, "chip/")
plot_dir <- paste0(cplot_dir, "GRanges/")

# GRanges objects
gr_hg <- readRDS(paste0(gr_dir, "human_batch1_grlist_", date, ".RDS"))
gr_mm <- readRDS(paste0(gr_dir, "mouse_batch1_grlist_", date, ".RDS"))

# cCRE tables -> GR objects
ccre_hg <- read.delim(ccre_path_hg, stringsAsFactors = FALSE)
ccre_hg <- makeGRangesFromDataFrame(ccre_hg, keep.extra.columns = TRUE)
ccre_hg$Group <- str_replace_all(ccre_hg$Group, ",|-", "_")

ccre_mm <- read.delim(ccre_path_mm, stringsAsFactors = FALSE)
ccre_mm <- makeGRangesFromDataFrame(ccre_mm, keep.extra.columns = TRUE)
ccre_mm$Group <- str_replace_all(ccre_mm$Group, ",|-", "_")

# batch 1 ChIP-seq meta
meta <- read.delim(paste0(meta_dir, "Chipseq/batch1_chip_meta_final_", date, ".tsv"), stringsAsFactors = FALSE)

stopifnot(all(meta$Experiment_ID %in% c(names(gr_hg), names(gr_mm))))

meta_hg <- meta %>% 
  filter(Experiment_ID %in% names(gr_hg)) %>% 
  distinct(Experiment_ID, .keep_all = TRUE) %>% 
  arrange(match(Experiment_ID, names(gr_hg)))

meta_mm <- meta %>% 
  filter(Experiment_ID %in% names(gr_mm)) %>%
  distinct(Experiment_ID, .keep_all = TRUE) %>% 
  arrange(match(Experiment_ID, names(gr_mm)))

stopifnot(identical(meta_hg$Experiment_ID, names(gr_hg)))
stopifnot(identical(meta_mm$Experiment_ID, names(gr_mm)))


# Re-size peaks to fixed size and create subset of cCRE filtered.
# ------------------------------------------------------------------------------


window_size <- 150

# Fix range to the summit position and pad with a fixed window

summit_window <- function(gr, window_size) {
  start(gr) <- end(gr) <- start(gr) + gr$Summit_from_start
  gr <- gr + window_size
  return(gr)
}

gr_rs_hg <- GRangesList(lapply(gr_hg, summit_window, window_size))
gr_rs_mm <- GRangesList(lapply(gr_mm, summit_window, window_size))


# Subset GRs to those overlapping cCREs

gr_sub_hg <- GRangesList(lapply(gr_rs_hg, function(x) {
  x[unique(findOverlaps(x, ccre_hg)@from)]
}))

gr_sub_mm <- GRangesList(lapply(gr_rs_mm, function(x) {
  x[unique(findOverlaps(x, ccre_mm)@from)]
}))


# Proportion of peak overlap between pairs of experiments
# ------------------------------------------------------------------------------


get_prop_mat <- function(exp_gr, cores) {
  
  l <- mclapply(names(exp_gr), function(i) {
    
    unlist(lapply(names(exp_gr), function(j) {
      if (i == j) return(NA)
      ol <- findOverlaps(exp_gr[[i]], exp_gr[[j]])
      n_distinct(ol@from) / length(exp_gr[[i]])
    }))
    
  }, mc.cores = cores)
 
  prop_mat <- do.call(cbind, l)
  colnames(prop_mat) <- rownames(prop_mat) <- names(exp_gr)
  
  return(prop_mat)
}


# Slow and ended up using Jaccard

# prop_list <- list(
#   Human = get_prop_mat(gr_rs_hg, cores),
#   Human_sub = get_prop_mat(gr_sub_hg, cores),
#   Mouse = get_prop_mat(gr_rs_mm, cores),
#   Mouse_sub = get_prop_mat(gr_sub_mm, cores)
# )


# Jaccard coefficient between pairs of experiments. First build global matrix
# of all regions bound at least once
# ------------------------------------------------------------------------------


# region by experiment matrix of overlap with these regions


ol_matrix <- function(all_gr, gr_list, ncores = cores) {
  
  ol_vec <- rep(0, length(all_gr))
  
  ol_l <- mclapply(gr_list, function(x) {
    ol <- findOverlaps(x, all_gr)
    ol_vec[unique(ol@to)] <- 1
    return(ol_vec)
  }, mc.cores = cores)
  
  ol_mat <- do.call(cbind, ol_l)
  colnames(ol_mat) <- names(gr_list)
  rownames(ol_mat) <- paste(seqnames(all_gr), start(all_gr), end(all_gr), sep = ":")
  return(ol_mat)
}


# get_jaccard

get_jac_mat <- function(ol_mat) {
  
  jac_mat <- matrix(NA, ncol = ncol(ol_mat), nrow = ncol(ol_mat))
  colnames(jac_mat) <- rownames(jac_mat) <- colnames(ol_mat)
  
  for (i in 1:ncol(ol_mat)) {
    for (j in 1:ncol(ol_mat)) {
      if (i == j) next
      jac_mat[i, j] <- jaccard(ol_mat[, i], ol_mat[, j])
    }
  }
  
  return(jac_mat)
}


# List of global bound regions +/- cCRE filter

all_list <- list(
  Human = reduce(unlist(gr_rs_hg)),
  Human_sub = reduce(unlist(gr_sub_hg)),
  Mouse = reduce(unlist(gr_rs_mm)),
  Mouse_sub = reduce(unlist(gr_sub_mm))
)


# Binary matrices of region overlap. Here using the full region as reference for
# both +/- cCRE - when using Jaccard, this gives the same result as using the 
# cCRE region as reference for the +cCRE set. 

olmat_list <- list(
  Human = ol_matrix(all_list$Human, gr_rs_hg),
  Human_sub = ol_matrix(all_list$Human, gr_sub_hg),
  Mouse = ol_matrix(all_list$Mouse, gr_rs_mm),
  Mouse_sub = ol_matrix(all_list$Mouse, gr_sub_mm)
)


# Jaccard between experiments (slow!)

outfile <- file.path(scratch_dir, "gr_bind_experiment_jaccard_Jan2023.RDS")

if (!file.exists(outfile)) {
  jac_list <- lapply(olmat_list, get_jac_mat)
  saveRDS(jac_list, outfile)
} else {
  jac_list <- readRDS(outfile)
}


# Summary of intra- vs inter-TR Jaccard coefs of experiment peak overlap 

jac_df <- lapply(jac_list, format_pair_df, symmetric = TRUE)
jac_summary <- lapply(jac_df, get_summary)

# Violin boxplots of Jaccard +/- cCRE filter

p1 <- stat_vboxplot(jac_df$Human, y_var = "Value", y_name = "Jaccard coef of Overlap", title = "All peaks")
p2 <- stat_vboxplot(jac_df$Human_sub, y_var = "Value", y_name = "Jaccard coef Overlap", title = "cCRE peaks")
plot_grid(p1, p2)

# Scatter of Jaccard for each experiment pair +/- cCRE filter

plot_df <-
  left_join(
    x = mutate(jac_df$Human, ID = paste(Row, Col, sep = ":")),
    y = mutate(jac_df$Human_sub, ID = paste(Row, Col, sep = ":")),
    by = "ID"
  ) %>%
  dplyr::rename(Jaccard_nocCRE = Value.x,
                Jaccard_cCRE = Value.y)

cor(plot_df$Jaccard_nocCRE, plot_df$Jaccard_cCRE)


p3 <- 
  ggplot(plot_df, aes(x = Jaccard_nocCRE, y = Jaccard_cCRE)) +
  geom_point(alpha = 0.4, shape = 21, size = 2.4) +
  xlab("Jaccard -cCRE filter") +
  ylab("Jaccard +cCRE filter") +
  ggtitle("Similarity between human ChIP-seq experiments +/- cCRE filter") +
  theme_classic() +
  theme(axis.text = element_text(size = 25),
        axis.title = element_text(size = 25),
        plot.title = element_text(size = 25))


# Hclust -> dendogram object for plotting experiment similarity +/- ccRE filter
# ------------------------------------------------------------------------------


get_hc <- function(mat) {
  
  # Reduce names for plotting
  names <- str_split(colnames(mat), "_", simplify = TRUE)[, 1:2]
  names <- make.unique(paste(names[, 1], names[, 2], sep = "_"))
  colnames(mat) <- rownames(mat) <- names
  
  hc <- hclust(as.dist(1 - mat), method = "ward.D")
  return(hc)
}


jac_hc_list <- lapply(jac_list, function(x) as.dendrogram(get_hc(x)))



# Plot dendogram, with rectangles over k cuts and the labels coloured by TR

plot_hc <- function(dend, tf_pal, k = 8) {
  
  # Isolate TF names for label colours
  tf_cols <- str_split(labels(dend), "_", simplify = TRUE)[, 2]
  tf_cols <- str_to_upper(str_replace(tf_cols, "\\..*", ""))
  tf_cols <- tf_pal[tf_cols]
  
  labels_colors(dend) <- tf_cols
  dend %>% set("labels_cex", 0.8) %>% plot()
  rect.dendrogram(dend, k = k)
}


# -cCRE filter
plot_hc(jac_hc_list$Human, tf_pal = tf_pal_hg)

# +cCRE filter
plot_hc(jac_hc_list$Human_sub, tf_pal = tf_pal_hg)


# Compare the group membership at the same cut +/- cCRE filter: get Jaccard of
# experiments between clusters generated +/- cCRE filter at the same k cutoff
# ------------------------------------------------------------------------------


# Return a matrix of Jaccard coefs for experiment overlap between clusters +/-
# filtered for cCRE overlap

get_cluster_jac <- function(ct_list, species = "Human", k = 8) {
  
  mat <- matrix(NA, ncol = k, nrow = k)
  rownames(mat) <- paste0("No_cCRE_", 1:k)
  colnames(mat) <- paste0("cCRE_", 1:k)
  
  for (i in 1:k) {
    for (j in 1:k) {
      noccre <- names(ct_list[[species]][jac_ct_list[[species]] == i])
      ccre <- names(ct_list[[paste0(species, "_sub")]][ct_list[[paste0(species, "_sub")]] == j])
      mat[i, j] <- length(intersect(noccre, ccre)) / length(union(noccre, ccre))
    }
  }
  return(mat)
}



# Cut dendogram at k=8 (there are 8 TFs/TRs)
jac_ct_list <- lapply(jac_hc_list, cutree, k = 8)

# Matrix of Jaccard coef of experiment overlap between clusters at given cut
cluster_jac_hg <- get_cluster_jac(jac_ct_list, species = "Human", k = 8)



# Heatmap that shows the overlap (Jaccard) of experiments between clusters

cluster_jac_heatmap <- function(mat) {
  
  pheatmap(mat, 
           color = bluered_pal,
           display_numbers = TRUE,
           number_color = "black",
           fontsize_number = 20,
           border_color = "black")
  
}


cluster_jac_heatmap(cluster_jac_hg)
