## Plot different processing strategies for ChIP-seq bind matrix, and save
## a final processed object for downstream analysis
## -----------------------------------------------------------------------------

library(tidyverse)
library(preprocessCore)
library(scales)
library(cowplot)

date <- "Apr2022"  # most recent data freeze
intensity_flag <- FALSE  # use ouyang scores with or without macs2 intensity
binary_dist <- 25e3  # distance threshold used for binary gene scores
peakset <- "idr"  # idr for TF (but overlap still for mecp2) or all overlap
plot_dir <- "~/Plots/Chipseq/Preprocess/"  # plot export dir
mat_dir <- "~/Data/Annotated_objects/Bind_matrices/Encpipe/" # bind matrix path

# batch 1 ChIP-seq meta and QC info
meta <- read.delim(paste0("~/Data/Metadata/Chipseq/batch1_chip_meta_final_", date, ".tsv"), stringsAsFactors = FALSE)

# batch 1 binding matrices
# binmat = binary matrix based on distance thresholds.
# bsmat = binding score matrix
binmat_hg <- readRDS(paste0(mat_dir, "binary_refseq_human_batch1_", date, "_", "distance=", binary_dist/1e3, "kb_uniquesymbol=TRUE_peakset=", peakset, ".RDS"))
binmat_mm <- readRDS(paste0(mat_dir, "binary_refseq_mouse_batch1_", date, "_", "distance=", binary_dist/1e3, "kb_uniquesymbol=TRUE_peakset=", peakset, ".RDS"))
binmat_ortho <- readRDS(paste0(mat_dir, "binary_refseq_ortho_batch1_", date, "_", "distance=", binary_dist/1e3, "kb_uniquesymbol=TRUE_peakset=", peakset, ".RDS"))
bsmat_hg <- readRDS(paste0(mat_dir, "ouyang_refseq_human_batch1_", date, "_", "dc=5000_intensity=", intensity_flag, "_peakset=", peakset, ".RDS"))
bsmat_mm <- readRDS(paste0(mat_dir, "ouyang_refseq_mouse_batch1_", date, "_", "dc=5000_intensity=", intensity_flag, "_peakset=", peakset, ".RDS"))
bsmat_ortho <- readRDS(paste0(mat_dir, "ouyang_refseq_ortho_batch1_", date, "_", "dc=5000_intensity=", intensity_flag, "_peakset=", peakset, ".RDS"))


bs_list <- list("Human" = bsmat_hg,
                "Mouse" = bsmat_mm,
                "Ortho" = bsmat_ortho)

bin_list <- list("Human" = binmat_hg,
                 "Mouse" = binmat_mm,
                 "Ortho" = binmat_ortho)


if (peakset == "idr") {
  min_peaks <- 100
} else if (peakset == "overlap") {
  min_peaks <- 1000
  meta <- mutate(meta, N_peaks = N_overlap_peaks)
}


# Remove genes with no binding across all experiments - as some of the raw bind
# scores have row sum to ^-230, round to coerce to 0.
#-------------------------------------------------------------------------------


rowsum_filt <- 1
round_var = 10


filter_nobind <- function(mat, rowsum_arg, round_arg) {
  # remove rows that have no binding (all entries are 0)
  mat <- round(mat, round_arg)
  no_bind <- which(rowSums(mat) < rowsum_arg)
  if (length(no_bind) > 0) mat <- mat[-no_bind,]
  return(mat)
}


bs_sums <- lapply(bs_list, rowSums)
bin_sums <- lapply(bin_list, rowSums)

plot(density(bs_sums$Human))
abline(v = rowsum_filt)
bs_low <- lapply(bs_sums, function(x) sum(x < rowsum_filt))

# Remove no bind
bs_filt <- lapply(bs_list, filter_nobind, rowsum_filt, round_var)
bin_filt <- lapply(bin_list, filter_nobind, rowsum_filt, round_var)

# lapply(bs_filt, dim)
# lapply(bin_filt, dim)


# Process continuous binding scores
# 1) Raw scores 
# 2) Log10+1 
# 3) Quantile norm log10+1 **works on samples, not features**
# 4) Standardize log10+1 (features mean=0 sd=1)
# 5) Min max norm log10+1 (features bound [0, 1])
# 6) Scale after Quantile norm log10+1
#-------------------------------------------------------------------------------


qn_and_name <- function(mat) {
  # performs quantile and norm and label mat with original row/col names
  qn_mat <- preprocessCore::normalize.quantiles(mat) # default columns as samples
  rownames(qn_mat) <- rownames(mat)
  colnames(qn_mat) <- colnames(mat)
  return(qn_mat)
}


minmax_norm <- function(mat) {
  # normalize rows to be bound [0,1]
  t(apply(mat, 1, function(x) (x - min(x)) / (max(x) - min(x))))
}


process_list <- function(mat) {
  # take in a matrix of features x samples and return a list of different
  # processing steps
  list(
    Raw = mat,
    Log = log10(mat+1),
    QN_log = qn_and_name(log10(mat + 1)),
    Scale_log = t(scale(t(log10(mat + 1)))),
    Minmax = minmax_norm(log10(mat + 1)),
    Scale_QN_log = t(scale(t(qn_and_name(log10(mat + 1)))))
  )
}


bs_norm_hg <- process_list(bs_filt$Human)
bs_norm_mm <- process_list(bs_filt$Mouse)
bs_norm_ortho <- process_list(bs_filt$Ortho)


# Get the rank of each genes mean bind score across the different normalization methods
#-------------------------------------------------------------------------------


rank_df <- function(bind_list) {
  # return a data frame where rows represent genes and columns the different
  # processing schemes in bind_list, where elements are the rank of that gene's
  # average bind score across experiments. lower rank = higher average score
  rank_list <- lapply(bind_list, function(x) {
    rank(rowMeans(-x), ties.method = "min")
  })
  do.call(cbind, rank_list)
}


rank_hg <- rank_df(bs_norm_hg)
rank_mm <- rank_df(bs_norm_mm)
rank_ortho <- rank_df(bs_norm_ortho)

cor(rank_hg)
cor(rank_mm)
cor(rank_ortho)


# Save out binary+processed together

bs_norm_hg$Binary <- bin_filt$Human
bs_norm_mm$Binary <- bin_filt$Mouse
bs_norm_ortho$Binary <- bin_filt$Ortho

saveRDS(bs_norm_hg, 
        file = paste0(mat_dir, "Human_refseq_", date, "_processed_bindmat_list_minpeak=", min_peaks, "_ouyang_dc=5000_intensity=", intensity_flag, "_binary=", binary_dist/1e3, "kb_peakset=", peakset, ".RDS"))
saveRDS(bs_norm_mm, 
        file = paste0(mat_dir, "Mouse_refseq_", date, "_processed_bindmat_list_minpeak=", min_peaks, "_ouyang_dc=5000_intensity=", intensity_flag, "_binary=", binary_dist/1e3, "kb_peakset=", peakset, ".RDS"))
saveRDS(bs_norm_ortho, 
        file = paste0(mat_dir, "Ortho_refseq_", date, "_processed_bindmat_list_minpeak=", min_peaks, "_ouyang_dc=5000_intensity=", intensity_flag, "_binary=", binary_dist/1e3, "kb_peakset=", peakset, ".RDS"))


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
  xlab(paste0("Log 10 count of ", peakset, " peaks")) +
  ylab("Count") +
  theme(axis.title = element_text(size = 15),
        axis.text = element_text(size = 15))

ggsave(p1, device = "png", dpi = 300, height = 6, width = 8,
       filename = paste0(plot_dir, "hist_count_", peakset, "_minpeak=", min_peaks, ".png"))


# Scatter of peak counts with count of binary gene assignments

if (peakset == "idr") {
  plot_df2 <- data.frame(
    Experiment_ID = colnames(bin_list$Ortho),
    Count_genes = colSums(bin_list$Ortho)
  ) %>%
    left_join(meta[, c("Experiment_ID", "N_peaks")], by = "Experiment_ID")
} else if (peakset == "overlap") {
  plot_df2 <- data.frame(
    Experiment_ID = colnames(bin_list$Ortho),
    Count_genes = colSums(bin_list$Ortho)
  ) %>%
    left_join(meta[, c("Experiment_ID", "N_overlap_peaks")], by = "Experiment_ID") %>%
    dplyr::rename(N_peaks = N_overlap_peaks)
}


p2 <- 
  ggplot(plot_df2, aes(x = log10(Count_genes), y = log10(N_peaks))) +
  geom_point() +
  geom_hline(yintercept = log10(min_peaks), col = "red") +
  theme_classic() +
  xlab(paste0("Log 10 count of genes (binary +/- ", binary_dist/1e3, "kb)")) +
  ylab(paste0("Log 10 count of ", peakset, " peaks")) +
  theme(axis.title = element_text(size = 15),
        axis.text = element_text(size = 15))

ggsave(p2, device = "png", dpi = 300, height = 6, width = 8,
       filename = paste0(plot_dir, "ortho_", peakset, "_peaks_vs_binary_genes_", binary_dist/1e3, "kb_minpeak=", min_peaks, ".png"))



# Density plots of gene/row sums of raw+log+binary binding scores before filtering


plist1 <- lapply(bs_sums, function(x) {
  ggplot(data.frame(Rowsum = x), aes(x = Rowsum)) +
    geom_density() +
    xlab("Rowsum of raw continuous scores") +
    theme_classic() +
    scale_x_continuous(labels = scales::scientific) +
    theme(axis.title = element_text(size = 12),
          axis.text = element_text(size = 12))
})


rowSums(log10(bs_list$Human+1)) %>% density %>% plot


plist2 <- lapply(bs_list, function(x) {
  data.frame(Rowsum = rowSums(log10(x+1))) %>% 
  ggplot(., aes(x = Rowsum)) +
    geom_density() +
    xlab("Rowsum of log10+1 continuous scores") +
    theme_classic() +
    theme(axis.title = element_text(size = 12),
          axis.text = element_text(size = 12))
})


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
       filename = paste0(plot_dir, "genewise_rowsum_bind_density_minpeak=", min_peaks, "_ouyang_dc=5000_intensity=", intensity_flag, "_binary=", binary_dist/1e3, "kb_peakset=", peakset, ".png"))


# density of bind scores of sampled genes

set.seed(7)
samp_genes <- sample(rownames(bs_norm_hg$Raw), 300)

plot(density(bs_norm_hg$QN_log[samp_genes[1],]), xlim = c(0, 3), ylim = c(0, 10))
for (i in samp_genes[2:length(samp_genes)]) {
  lines(density(bs_norm_hg$QN_log[i,]))
}

# Sample wise density over each processing scheme

plot_density <- function(mat, ptitle) {
  
  data.frame(mat) %>% 
    rownames_to_column(var = "ID") %>%
    pivot_longer(cols = 1:ncol(mat)+1) %>% 
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
         filename = paste0(plot_dir, "experimentwise_bind_density_norm=", names(bs_norm_ortho)[i], "_minpeak=", min_peaks, "_ouyang_dc=5000_intensity=", intensity_flag, "_binary=", binary_dist/1e3, "kb_peakset=", peakset, ".png"))
}


# density of samples top scores

top_n <- 1000

plot(density(head(sort(bs_norm_hg$Scale_QN_log[,1], decreasing = TRUE), top_n)),
     xlim = c(0, 10), ylim = c(0, 10))
for (i in 2:ncol(bs_norm_hg$Scale_QN_log)) {
  lines(density(head(sort(bs_norm_hg$Scale_QN_log[,i], decreasing = TRUE), top_n)))
}
