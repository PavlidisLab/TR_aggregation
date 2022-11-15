## Plot different processing strategies for ChIP-seq bind matrix, and save
## a final processed object for downstream analysis
## -----------------------------------------------------------------------------

library(tidyverse)
library(preprocessCore)
library(scales)
library(cowplot)
source("R/setup-01_config.R")

binary_dist <- 25e3  # window size in base pairs to consider a peak->TSS assignment
plot_dir <- paste0(cplot_dir, "Preprocess/")

# ChIP-seq meta and QC info
meta <- read.delim(paste0(meta_dir, "Chipseq/batch1_chip_meta_final_", date, ".tsv"), stringsAsFactors = FALSE)

# binmat = binary matrix based on distance thresholds.
# bsmat = binding score matrix
binmat_hg <- readRDS(paste0(cmat_dir, "binary_refseq_human_batch1_", date, "_", "distance=", binary_dist/1e3, "kb.RDS"))
binmat_mm <- readRDS(paste0(cmat_dir, "binary_refseq_mouse_batch1_", date, "_", "distance=", binary_dist/1e3, "kb.RDS"))
binmat_ortho <- readRDS(paste0(cmat_dir, "binary_refseq_ortho_batch1_", date, "_", "distance=", binary_dist/1e3, "kb.RDS"))
bsmat_hg <- readRDS(paste0(cmat_dir, "ouyang_refseq_human_batch1_", date, "_", "dc=5000_intensity=FALSE.RDS"))
bsmat_mm <- readRDS(paste0(cmat_dir, "ouyang_refseq_mouse_batch1_", date, "_", "dc=5000_intensity=FALSE.RDS"))
bsmat_ortho <- readRDS(paste0(cmat_dir, "ouyang_refseq_ortho_batch1_", date, "_", "dc=5000_intensity=FALSE.RDS"))

bs_list <- list("Human" = bsmat_hg,
                "Mouse" = bsmat_mm,
                "Ortho" = bsmat_ortho)

bin_list <- list("Human" = binmat_hg,
                 "Mouse" = binmat_mm,
                 "Ortho" = binmat_ortho)


# Remove genes with no/trivial binding across all experiments. Note that some 
# genes have row sums == 1e-86, effectively 0. Mouse has more 0 bound genes.
# Here, using an arbitrary rounding threshold and rowsum filter based on 
# inspecting the distribution of rowsums. As implemented, rounding threshold is 
# redundant as rowsum filter removes genes anyways. Also note that a more 
# stringent rowsum filter is applied prior to the limma analysis in 05_
#-------------------------------------------------------------------------------


# remove rows that have no binding 

filter_nobind <- function(mat, rowsum_arg, round_arg = 10) {
  
  mat <- round(mat, round_arg)
  no_bind <- which(rowSums(mat) < rowsum_arg)
  if (length(no_bind) > 0) mat <- mat[-no_bind,]
  
  return(mat)
}


# Remove no bind
bs_filt <- list(
  Human = filter_nobind(bs_list$Human, rowsum_arg = 1),
  Mouse = filter_nobind(bs_list$Mouse, rowsum_arg = 3),
  Ortho = filter_nobind(bs_list$Ortho, rowsum_arg = 1)
)

bin_filt <- lapply(bin_list, function(x) x[rowSums(x) != 0, ])


# Sums before and after filter
bs_sums <- lapply(bs_list, rowSums)
bs_filt_sums <- lapply(bs_filt, rowSums)
bin_sums <- lapply(bin_list, rowSums)
bin_filt_sums <- lapply(bin_filt, rowSums)


# Human: 398 genes removed; Mouse: 3436; Ortho: 13

bs_low <- lapply(names(bs_list), function(x) {
  bs_sums[[x]][setdiff(names(bs_sums[[x]]), names(bs_filt_sums[[x]]))]
})

n_rm <- lapply(bs_low, length)


# Inspect distributions before/after


hist(bs_sums$Human, breaks = 100)
hist(bs_filt_sums$Human, breaks = 100)
plot(density(log10(bs_sums$Human + 1)))
plot(density(log10(bs_filt_sums$Human + 1)))

hist(bs_sums$Mouse, breaks = 100)
hist(bs_filt_sums$Mouse, breaks = 100)
plot(density(log10(bs_sums$Mouse + 1)))
plot(density(log10(bs_filt_sums$Mouse + 1)))

hist(bs_sums$Ortho, breaks = 100)
hist(bs_filt_sums$Ortho, breaks = 100)
plot(density(log10(bs_sums$Ortho + 1)))
plot(density(log10(bs_filt_sums$Ortho + 1)))


# Process continuous binding scores
# 1) Raw scores 
# 2) Log10+1 
# 3) Quantile norm log10+1 **works on samples, not features**
# 4) Standardize log10+1 (features mean=0 sd=1)
# 5) Min max norm log10+1 (features bound [0, 1])
# 6) Scale features after Quantile norm log10+1
#-------------------------------------------------------------------------------


minmax_norm <- function(mat) {
  # normalize rows to be bound [0,1]
  t(apply(mat, 1, function(x) (x - min(x)) / (max(x) - min(x))))
}


# take in a matrix of genes x samples and return a list of different
# processing steps

process_list <- function(mat) {

  list(
    Raw = mat,
    Log = log10(mat + 1),
    QN_log = normalize.quantiles(log10(mat + 1), keep.names = TRUE),
    Scale_log = t(scale(t(log10(mat + 1)))),
    Minmax = minmax_norm(log10(mat + 1)),
    Scale_QN_log = t(scale(t(normalize.quantiles(log10(mat + 1), keep.names = TRUE))))
  )
}


bs_norm_hg <- process_list(bs_filt$Human)
bs_norm_mm <- process_list(bs_filt$Mouse)
bs_norm_ortho <- process_list(bs_filt$Ortho)


# Get the rank of each genes mean bind score across the different normalization methods
#-------------------------------------------------------------------------------

# return a data frame where rows represent genes and columns the different
# processing schemes in bind_list, where elements are the rank of that gene's
# average bind score across experiments. lower rank = higher average score

rank_df <- function(bind_list) {

  rank_list <- lapply(bind_list, function(x) {
    rank(rowMeans(-x), ties.method = "min")
  })
  do.call(cbind, rank_list)
}


rank_hg <- rank_df(bs_norm_hg)
rank_mm <- rank_df(bs_norm_mm)
rank_ortho <- rank_df(bs_norm_ortho)

cor_hg <- cor(rank_hg)
cor_mm <- cor(rank_mm)
cor_ortho <- cor(rank_ortho)


# Save out binary+processed together

bs_norm_hg$Binary <- bin_filt$Human
bs_norm_mm$Binary <- bin_filt$Mouse
bs_norm_ortho$Binary <- bin_filt$Ortho


saveRDS(bs_norm_hg, 
        file = paste0(cmat_dir, "Human_refseq_", date, "_processed_bindmat_list_minpeak=", min_peaks, "_ouyang_dc=5000_intensity=FALSE_binary=", binary_dist/1e3, "kb.RDS"))
saveRDS(bs_norm_mm, 
        file = paste0(cmat_dir, "Mouse_refseq_", date, "_processed_bindmat_list_minpeak=", min_peaks, "_ouyang_dc=5000_intensity=FALSE_binary=", binary_dist/1e3, "kb.RDS"))
saveRDS(bs_norm_ortho, 
        file = paste0(cmat_dir, "Ortho_refseq_", date, "_processed_bindmat_list_minpeak=", min_peaks, "_ouyang_dc=5000_intensity=FALSE_binary=", binary_dist/1e3, "kb.RDS"))


# Plots
# ------------------------------------------------------------------------------


# Hist of peak counts with min peak as marker
# NOTE: In earlier version this was used on all experiments - peak filtering
# now done earlier.

p1 <- 
  ggplot(meta, aes(x = log10(N_peaks))) +
  geom_histogram(bins = 50, fill = "slategray4") +
  geom_vline(xintercept = log10(min_peaks)) +
  theme_classic() +
  xlab(paste0("Log10 count of peaks")) +
  ylab("Count") +
  theme(axis.title = element_text(size = 15),
        axis.text = element_text(size = 15))

ggsave(p1, device = "png", dpi = 300, height = 6, width = 8,
       filename = paste0(plot_dir, "hist_count_minpeak=", min_peaks, ".png"))


# Scatter of peak counts with count of binary gene assignments

plot_df2 <- data.frame(
  Experiment_ID = colnames(bin_list$Ortho),
  Count_genes = colSums(bin_list$Ortho)
) %>%
  left_join(meta[, c("Experiment_ID", "N_peaks")], by = "Experiment_ID")


p2 <- 
  ggplot(plot_df2, aes(x = log10(Count_genes), y = log10(N_peaks))) +
  geom_point() +
  geom_hline(yintercept = log10(min_peaks), col = "red") +
  theme_classic() +
  xlab(paste0("Log 10 count of genes (binary +/- ", binary_dist/1e3, "kb)")) +
  ylab(paste0("Log 10 count of peaks")) +
  theme(axis.title = element_text(size = 15),
        axis.text = element_text(size = 15))

ggsave(p2, device = "png", dpi = 300, height = 6, width = 8,
       filename = paste0(plot_dir, "ortho_peaks_vs_binary_genes_", binary_dist/1e3, "kb_minpeak=", min_peaks, ".png"))


# Density plots of gene/row sums of raw+log+binary binding scores before filtering


# Raw rowsums
plist1 <- lapply(bs_sums, function(x) {
  ggplot(data.frame(Rowsum = x), aes(x = Rowsum)) +
    geom_density() +
    xlab("Rowsum of raw continuous scores") +
    theme_classic() +
    scale_x_continuous(labels = scales::scientific) +
    theme(axis.title = element_text(size = 12),
          axis.text = element_text(size = 12))
})



# Log10 rowsums
plist2 <- lapply(bs_list, function(x) {
  data.frame(Rowsum = rowSums(log10(x + 1))) %>% 
  ggplot(., aes(x = Rowsum)) +
    geom_density() +
    xlab("Rowsum of log10+1 continuous scores") +
    theme_classic() +
    theme(axis.title = element_text(size = 12),
          axis.text = element_text(size = 12))
})


# Binary rowsums
plist3 <- lapply(bin_sums, function(x) {
  ggplot(data.frame(Rowsum = x), aes(x = Rowsum)) +
    geom_density() +
    xlab("Rowsum of binary scores") +
    theme_classic() +
    theme(axis.title = element_text(size = 12),
          axis.text = element_text(size = 12))
})


p3 <- plot_grid(plist1$Human, plist2$Human, plist3$Human,
                plist1$Mouse, plist2$Mouse, plist3$Mouse,
                plist1$Ortho, plist2$Ortho, plist3$Ortho,
                nrow = 3)

ggsave(p3, device = "png", dpi = 300, height = 9, width = 12,
       filename = paste0(plot_dir, "genewise_rowsum_bind_density_minpeak=", min_peaks, "_ouyang_dc=5000_intensity=FALSE_binary=", binary_dist/1e3, "kb.png"))



# Sample wise density over each processing scheme



plot_density <- function(mat, ptitle) {
  
  data.frame(mat) %>% 
    rownames_to_column(var = "ID") %>%
    pivot_longer(cols = 1:ncol(mat) + 1) %>% 
    ggplot(aes(value)) +
    geom_density(aes(group = ID), size = 0.3) +
    xlab("Bind scores") +
    ggtitle(ptitle) +
    theme_classic() +
    theme(axis.text = element_text(size = 15),
          title = element_text(size = 15),
          plot.title = element_text(hjust = 0.5))
}


plist4 <- lapply(1:length(bs_norm_hg), function(x) {
  plot_density(t(bs_norm_hg[[x]]), ptitle = paste("Human", names(bs_norm_hg)[x]))
})
names(plist4) <- names(bs_norm_hg)


plist5 <- lapply(1:length(bs_norm_mm), function(x) {
  plot_density(t(bs_norm_mm[[x]]), ptitle = paste("Mouse", names(bs_norm_mm)[x]))
})
names(plist5) <- names(bs_norm_mm)


plist6 <- lapply(1:length(bs_norm_ortho), function(x) {
  plot_density(t(bs_norm_ortho[[x]]), ptitle = paste("Ortho", names(bs_norm_ortho)[x]))
})
names(plist6) <- names(bs_norm_ortho)


for (i in 1:length(plist4)) {
  pi <- plot_grid(plist4[[i]], plist5[[i]], plist6[[i]], nrow = 1)
  ggsave(pi, device = "png", dpi = 300, height = 6, width = 12,
         filename = paste0(plot_dir, "experimentwise_bind_density_norm=", names(bs_norm_ortho)[i], "_minpeak=", min_peaks, "_ouyang_dc=5000_intensity=FALSE_binary=", binary_dist/1e3, "kb.png"))
}



# For inspecting density of sampled genes


# set.seed(7)
# samp_genes <- sample(rownames(bs_norm_hg$Raw), 300)
# 
# plot(density(bs_norm_hg$QN_log[samp_genes[1],]), xlim = c(0, 3), ylim = c(0, 10))
# for (i in samp_genes[2:length(samp_genes)]) {
#   lines(density(bs_norm_hg$QN_log[i,]))
# }


