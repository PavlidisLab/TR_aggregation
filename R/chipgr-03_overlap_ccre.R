## Explore overlap of ChIP-seq GR objects among ENCODE candidate cis reg elements
## TODO: re-consider design: script has redundancy in count calculations
## TODO: organize functions into GR utils
## -----------------------------------------------------------------------------

library(tidyverse)
library(GenomicRanges)
source("~/regnetR/R/utils/range_table_functions.R")

date <- "Apr2022"
peakset <- "idr"
pipeline_dir <- "/cosmos/data/pipeline-output/chipseq-encode-pipeline/chip"
plot_dir <- "~/Plots/Chipseq/GRanges/"

# GRanges objects
gr_hg <- readRDS(paste0("~/Data/Annotated_objects/GRanges/human_batch1_grlist_peakset=", peakset, "_", date, ".RDS"))
gr_mm <- readRDS(paste0("~/Data/Annotated_objects/GRanges/mouse_batch1_grlist_peakset=", peakset, "_", date, ".RDS"))

# cCRE tables -> GR objects
ccre_hg <- read.delim("~/Data/Chromosome_info/cCREs_V3_hg38.bed", stringsAsFactors = FALSE)
ccre_hg <- makeGRangesFromDataFrame(ccre_hg, keep.extra.columns = TRUE)
ccre_hg$Group <- str_replace_all(ccre_hg$Group, ",|-", "_")

ccre_mm <- read.delim("~/Data/Chromosome_info/cCREs_V3_mm10.bed", stringsAsFactors = FALSE)
ccre_mm <- makeGRangesFromDataFrame(ccre_mm, keep.extra.columns = TRUE)
ccre_mm$Group <- str_replace_all(ccre_mm$Group, ",|-", "_")

# batch 1 ChIP-seq meta
meta <- read.delim(paste0("~/Data/Metadata/Chipseq/batch1_chip_meta_final_", date, ".tsv"), stringsAsFactors = FALSE)
stopifnot(all(meta$Experiment_ID %in% c(names(gr_hg), names(gr_mm))))

meta_hg <- meta %>% 
  filter(Experiment_ID %in% names(gr_hg)) %>% 
  distinct(Experiment_ID, .keep_all = TRUE) %>% 
  arrange(match(Experiment_ID, names(gr_hg)))

meta_mm <- meta %>% 
  filter(Experiment_ID %in% names(gr_mm)) %>%
  distinct(Experiment_ID, .keep_all = TRUE) %>% 
  arrange(match(Experiment_ID, names(gr_mm)))

tfs_hg <- unique(meta_hg$Symbol)
tfs_mm <- unique(meta_mm$Symbol)

stopifnot(identical(meta_hg$Experiment_ID, names(gr_hg)))
stopifnot(identical(meta_mm$Experiment_ID, names(gr_mm)))


# Fix to 1bp summit for overlaps
# ------------------------------------------------------------------------------


window_size <- 0

summit_window <- function(gr, window_size) {
  # Fix range to the summit position and pad with a fixed window
  start(gr) <- end(gr) <- start(gr) + gr$Summit_from_start
  gr <- gr + window_size
  return(gr)
}

gr_rs_hg <- lapply(gr_hg, summit_window, window_size)
gr_rs_mm <- lapply(gr_mm, summit_window, window_size)


# Summarize how many cCRE elements are covered by each peak for each experiment 
# ------------------------------------------------------------------------------


count_peaks_with_ccre <- function(peaklist_gr, ccre_gr) {
  # Returns a named matrix with a row for each element in peaklist_gr and
  # columns corresponding to counts. Each element is the number of peaks that
  # overlap a cCRE
  
  count_list <- lapply(peaklist_gr, function(x) {
    countOverlaps(x, ccre_gr)
  })
  
  max_count <- max(unlist(lapply(count_list, max)))
  
  # count table of the overlaps, convert to matrix and tack on total count
  table_list <- lapply(count_list, function(x) {
    table(factor(x, levels = c(0:max_count)))  
  })
  
  mat <- as.matrix(do.call(rbind, table_list))
  colnames(mat) <- paste0("Count_", colnames(mat), "_overlapped")
  mat <- cbind(mat, Count_peaks = rowSums(mat))
  return(mat)
}


frac_peaks_with_ccre <- function(count_mat) {
  frac_mat <- count_mat
  frac_mat[, 1:ncol(frac_mat)-1] <- frac_mat[, 1:ncol(frac_mat)-1] / frac_mat[, "Count_peaks"]
  return(frac_mat)
}


# Human: All experiments overlap at least one cCRE. On average, 33% of peak 
# summits did not overlap a cCRE. GSE60006_HES1_Human_FLAG-HES1 has the least
# overlap (95% summits no overlap with a cCRE), but only has 40 peaks. The next
# highest was GSE125659_MECP2_Human_MECP2-KO_HISTONE (94%). As expected,
# MECP2 on average has the most peaks not overlapping a cCRE (84%)


ccre_counts_hg <- count_peaks_with_ccre(gr_rs_hg, ccre_hg)
ccre_frac_hg <- frac_peaks_with_ccre(ccre_counts_hg)

any_zero_hg <- which(ccre_frac_hg[, "Count_0_overlapped"] == 0)
max_zero_hg <- ccre_frac_hg[which.max(ccre_frac_hg[, "Count_0_overlapped"]), , drop = FALSE]
max_multi_hg <- names(which.max(ccre_frac_hg[, "Count_2_overlapped"]))

zero_by_tf_hg <- data.frame(ccre_frac_hg) %>%
  rownames_to_column(var = "Experiment_ID") %>%
  left_join(meta[, c("Experiment_ID", "Symbol")], by = "Experiment_ID") %>%
  group_by(Symbol) %>%
  summarize(Mean_0 = mean(Count_0_overlapped))

summary(ccre_frac_hg)


# Mouse: All experiments overlap at least one cCRE. On average, 50% of peak 
# summits did not overlap a cCRE. GSE139509_Mecp2_Mouse_Ab1-Mecp2-KO_HISTONE
# had the least overlap (all 209 peaks no cCRE overlap - note its a KO experiment).

# GSE60006_HES1_Human_FLAG-HES1 has the least
# overlap (95% summits no overlap with a cCRE), but only has 40 peaks. The next
# highest was GSE125659_MECP2_Human_MECP2-KO_HISTONE (94%). As expected,
# Mecp2 on average has the most peaks not overlapping a cCRE (96.5%)


ccre_counts_mm <- count_peaks_with_ccre(gr_rs_mm, ccre_mm)
ccre_frac_mm <- frac_peaks_with_ccre(ccre_counts_mm)

any_zero_mm <- which(ccre_frac_mm[, "Count_0_overlapped"] == 0)
max_zero_mm <- ccre_frac_mm[which.max(ccre_frac_mm[, "Count_0_overlapped"]), , drop = FALSE]
max_multi_mm <- names(which.max(ccre_frac_mm[, "Count_2_overlapped"]))

zero_by_tf_mm <- data.frame(ccre_frac_mm) %>% 
  rownames_to_column(var = "Experiment_ID") %>% 
  left_join(meta[, c("Experiment_ID", "Symbol")], by = "Experiment_ID") %>% 
  group_by(Symbol) %>% 
  summarize(Mean_0 = mean(Count_0_overlapped))

summary(ccre_frac_mm)


# Get counts by TF that each cCRE overlaps
# ------------------------------------------------------------------------------


count_ccre_with_tf <- function(gr_query, gr_subject, meta, tfs) {
  # Return a matrix with nrow = length of gr_query and ncol equal to length
  # of TFs. each row corresponds to a peak in gr_query and elements are the 
  # count of TF data sets from gr_subject that the peak overlaps
  
  # get the indices of the subject hits - requires coercion to GRList
  gr_ol <- findOverlaps(gr_query, GRangesList(gr_subject))
  ol_by_ix <- as(gr_ol, "List")
  
  # get the corresponding TF symbol of the subject and count
  ol_by_tf <- extractList(meta$Symbol, ol_by_ix)
  
  ol_by_count <- mclapply(ol_by_tf, function(x) {
    table(factor(x, levels = tfs))
  }, mc.cores = 8)
  
  # return the matrix of counts with rows named after the range
  count_mat <- as.matrix(do.call(rbind, ol_by_count))
  rownames(count_mat) <- 
    paste0(seqnames(gr_query), ":", start(gr_query), "-", end(gr_query))
  return(count_mat)
}


# Human
ccre_by_tf_hg <- count_ccre_with_tf(ccre_hg, gr_rs_hg, meta_hg, tfs_hg)
# ccre_by_tf_hg[order(ccre_by_tf_hg[, "ASCL1"], decreasing = TRUE),][1:30,]


# Mouse
ccre_by_tf_mm <- count_ccre_with_tf(ccre_mm, gr_mm, meta_mm, tfs_mm)
# ccre_by_tf_mm[order(ccre_by_tf_mm[, "Ascl1"], decreasing = TRUE),][1:30,]


# For each experiment, get the count of cCRE groups that were overlapped
# Note: The frac calculation is done by total peaks, not total assignments. A
# small amount of experiments have some peaks that are assigned to 2 cCREs due
# to being on the border
# ------------------------------------------------------------------------------


count_mat_by_group <- function(gr_query, gr_subject) {
  # gr_query assumed to be GRanges list of peaks, gr_subject to be cCRE GR.
  # returns a matrix of nrow = length(gr_query) and ncol = cCRE counts. each
  # element of the matrix is the number of cCRE elements for each group that was
  # overlapped by each experiment
  
  count_list <- mclapply(gr_query, function(x) {
    # get the indices of the subject hits
    gr_ol <- findOverlaps(x, gr_subject)
    ol_by_ix <- gr_ol@from
    # get the corresponding cCRE group of the subject and count
    ol_by_group <- unlist(extractList(gr_subject$Group, ol_by_ix))
    ol_count <- table(factor(ol_by_group, levels = unique(gr_subject$Group)))
  }, mc.cores = 8)
  # return the matrix of counts 
  count_mat <- as.matrix(do.call(rbind, count_list))
  return(count_mat)
}


# Human:

group_count_hg <- count_mat_by_group(gr_rs_hg, ccre_hg)

group_count_hg <- 
  cbind(None = ccre_counts_hg[, "Count_0_overlapped"], group_count_hg)

group_frac_hg <- group_count_hg / ccre_counts_hg[, "Count_peaks"]

group_by_tf_hg <- data.frame(group_frac_hg) %>% 
  rownames_to_column(var = "Experiment_ID") %>% 
  left_join(meta[, c("Experiment_ID", "Symbol")], by = "Experiment_ID") %>% 
  group_by(Symbol) %>% summarise(across(where(is.double), list(mean)))


# get summary of all TFs, separately considering Mecp2
summary(1 - group_frac_hg[filter(meta_hg, Symbol != "MECP2")$Experiment_ID, "None"])
summary(1 - group_frac_hg[filter(meta_hg, Symbol == "MECP2")$Experiment_ID, "None"])


# Mouse:

group_count_mm <- count_mat_by_group(gr_rs_mm, ccre_mm)

group_count_mm <- 
  cbind(None = ccre_counts_mm[, "Count_0_overlapped"], group_count_mm)

group_frac_mm <- group_count_mm / ccre_counts_mm[, "Count_peaks"]
  
group_by_tf_mm <- data.frame(group_frac_mm) %>% 
  rownames_to_column(var = "Experiment_ID") %>% 
  left_join(meta[, c("Experiment_ID", "Symbol")], by = "Experiment_ID") %>% 
  group_by(Symbol) %>% summarise(across(where(is.double), list(mean)))

# get summary of all TFs, separately considering Mecp2
summary(1 - group_frac_mm[filter(meta_mm, Symbol != "Mecp2")$Experiment_ID, "None"])
summary(1 - group_frac_mm[filter(meta_mm, Symbol == "Mecp2")$Experiment_ID, "None"])


# Plots
# ------------------------------------------------------------------------------


# boxplot of the fraction of peaks per experiment that did not overlap a cCRE,
# grouped by TF

boxplot(ccre_frac_hg[, "Count_0_overlapped"] ~ meta_hg$Symbol)

boxplot(ccre_frac_mm[, "Count_0_overlapped"] ~ meta_mm$Symbol)


# Stacked barchart of cCRE membership per experiment. Include total fraction
# of cCRE groups as a comparison


bin_groups <- function(str_vec) {
  # helper to group common cCREs regulatory element types, returning a factor
  
  str_vec <- str_replace_all(str_vec, "\\.|,|-", "_")
  
  group_vec <- vapply(as.character(str_vec), function(x) {
    if (x %in% c("pELS_CTCF_bound",
                 "dELS_CTCF_bound",
                 "dELS",
                 "pELS")) {
      return ("Enhancer-like")
    } else if (x %in% c("PLS_CTCF_bound",
                        "PLS")) {
      return ("Promoter-like")
    } else if (x %in% c("CTCF_only_CTCF_bound",
                        "DNase_H3K4me3_CTCF_bound",
                        "DNase_H3K4me3")) {
      return ("Other")
    } else {
      return ("None")
    }
    
  }, character(1))
  
  group_vec <- factor(
    group_vec, 
    levels = c("None", "Other", "Enhancer-like", "Promoter-like")
  )
}


all_ccre <- factor(
  c(unique(ccre_hg$Group), "None"),
  levels = c(
    "None",
    "PLS",
    "PLS_CTCF_bound",
    "pELS",
    "pELS_CTCF_bound",
    "dELS_CTCF_bound",
    "dELS",
    "DNase_H3K4me3",
    "DNase_H3K4me3_CTCF_bound",
    "CTCF_only_CTCF_bound"
  )
)


all_cols <- c('#d9d9d9', '#006d2c','#74c476','#7bccc4','#9ecae1','#4292c6', '#08519c', '#ff7f00','#cab2d6', '#fdbf6f')
bin_cols <- c('#d9d9d9','#6a3d9a','#1f78b4','#33a02c')


# Human

ccre_hg_df <- data.frame(Group = ccre_hg$Group, stringsAsFactors = FALSE) %>% 
  group_by(Group) %>% 
  summarize(Frac = n()/length(ccre_hg$Group)) %>% 
  ungroup() %>% 
  mutate(Plot_var = "1")
ccre_hg_df$Group <- factor(ccre_hg_df$Group, levels = levels(all_ccre))


plot_df_hg <- data.frame(group_frac_hg) %>% 
  rownames_to_column(var = "Experiment_ID") %>% 
  left_join(meta[, c("Symbol", "Experiment_ID", "N_peaks")], by = "Experiment_ID") %>%
  pivot_longer(cols = c(unique(ccre_hg$Group), "None")) %>% 
  mutate(name = factor(name, levels = levels(all_ccre)),
         Bin_group = bin_groups(name),
         Experiment_ID = factor(Experiment_ID, levels = unique(Experiment_ID)))

# All cCRE groups

p1a <- 
  plot_df_hg %>% 
  ggplot(aes(x = Experiment_ID, y = value, fill = name)) +
  facet_grid(~ Symbol, space = "free", scales = "free") +
  geom_bar(position = "stack", stat = "identity", width = 1) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.01))) +
  xlab("ChIP-seq experiment") +
  ylab("Proportion of peaks") +
  scale_fill_manual(values = all_cols, name = "cCRE group") +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 14),
        axis.text = element_text(size = 25),
        axis.title = element_text(size = 25),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20),
        legend.position = "bottom",
        plot.margin = margin(10, 10, 10, 10))  # add padding for removing legend


p1b <- p1a + theme(legend.position = "none",
                   strip.text = element_blank(),
                   axis.title.x = element_blank())

p1c <- 
  ccre_hg_df %>% 
  ggplot(aes(x = Plot_var, y = Frac, fill = Group)) +
  geom_bar(position = "stack", stat = "identity", width = 0.6) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.01))) +
  ylab("Proportion of all cCREs") +
  scale_fill_manual(values = all_cols[2:length(all_cols)], name = "cCRE group") +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text = element_text(size = 25),
        axis.title = element_text(size = 25),
        legend.position = "none",
        plot.margin = margin(10, 10, 10, 10))  # add padding for removing legend


ggsave(plot = p1a, width = 20, height = 12, dpi = 300, device = "png",
  filename = paste0(plot_dir, "cCRE_overlap_by_exp_human_legend_", date, "_peakset=", peakset, ".png"))

ggsave(plot = p1b, width = 20, height = 12, dpi = 300, device = "png",
       filename = paste0(plot_dir, "cCRE_overlap_by_exp_human_nolegend_", date, "_peakset=", peakset, ".png"))

ggsave(plot = p1c, width = 5, height = 12, dpi = 300, device = "png",
       filename = paste0(plot_dir, "cCRE_allfrac_human_", date, "_peakset=", peakset, ".png"))


# Binned cCRE groups

p2 <- 
  plot_df_hg %>% 
  ggplot(aes(x = Experiment_ID, y = value, fill = Bin_group)) +
  facet_grid(~ Symbol, space = "free", scales = "free") +
  geom_bar(position = "stack", stat = "identity") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.01))) +
  xlab("ChIP-seq experiment") +
  ylab("Proportion of peaks") +
  scale_fill_manual(values = bin_cols, name = "cCRE group") +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 12),
        axis.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15))


ggsave(plot = p2, width = 16, height = 9, dpi = 300, device = "png",
       filename = paste0(plot_dir, "bin_cCRE_overlap_by_exp_human_", date, "_peakset=", peakset, ".png"))


p3_list <- 
  lapply(tfs_hg, function(x) {
    p <- 
      plot_df_hg %>% 
      filter(Symbol == x) %>% 
      ggplot(aes(x = Experiment_ID, y = value, fill = Bin_group)) +
      geom_bar(position = "stack", stat = "identity") +
      scale_y_continuous(expand = expansion(mult = c(0, 0.01))) +
      xlab("ChIP-seq experiment") +
      ylab("Proportion of peaks") +
      ggtitle(x) +
      scale_fill_manual(values = bin_cols) +
      theme_classic() +
      theme(axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            plot.title = element_text(size = 12),
            axis.text = element_text(size = 15),
            axis.title = element_text(size = 15),
            axis.title.x = element_blank())
    if (x != tfs_hg[length(tfs_hg)]) {
      p <- p + theme(legend.position = "none")
    }
    return(p)
})

p3 <- cowplot::plot_grid(plotlist = p3_list)


# Mouse

ccre_mm_df <- data.frame(Group = ccre_mm$Group, stringsAsFactors = FALSE) %>% 
  group_by(Group) %>% 
  summarize(Frac = n()/length(ccre_mm$Group)) %>% 
  ungroup() %>% 
  mutate(Plot_var = "1")
ccre_mm_df$Group <- factor(ccre_mm_df$Group, levels = levels(all_ccre))


plot_df_mm <- data.frame(group_frac_mm) %>% 
  rownames_to_column(var = "Experiment_ID") %>% 
  left_join(meta[, c("Symbol", "Experiment_ID", "N_peaks")], by = "Experiment_ID") %>%
  pivot_longer(cols = c(unique(ccre_mm$Group), "None")) %>% 
  mutate(name = factor(name, levels = levels(all_ccre)),
         Bin_group = bin_groups(name),
         Experiment_ID = factor(Experiment_ID, levels = unique(Experiment_ID)))

# All cCRE groups

p4a <- 
  plot_df_mm %>% 
  ggplot(aes(x = Experiment_ID, y = value, fill = name)) +
  facet_grid(~ Symbol, space = "free", scales = "free") +
  geom_bar(position = "stack", stat = "identity", width = 1) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.01))) +
  xlab("ChIP-seq experiment") +
  ylab("Proportion of peaks") +
  scale_fill_manual(values = all_cols, name = "cCRE group") +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 14),
        axis.text = element_text(size = 25),
        axis.title = element_text(size = 25),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20),
        plot.margin = margin(10, 10, 10, 10))  # add padding for removing legend


p4b <- p4a + theme(legend.position = "none",
                   strip.text = element_blank(),
                   axis.title.x = element_blank())

p4c <- 
  ccre_mm_df %>% 
  ggplot(aes(x = Plot_var, y = Frac, fill = Group)) +
  geom_bar(position = "stack", stat = "identity", width = 0.6) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.01))) +
  ylab("Proportion of all cCREs") +
  scale_fill_manual(values = all_cols[2:length(all_cols)], name = "cCRE group") +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text = element_text(size = 25),
        axis.title = element_text(size = 25),
        legend.position = "none",
        plot.margin = margin(10, 10, 10, 10))  # add padding for removing legend



ggsave(plot = p4a, width = 20, height = 12, dpi = 300, device = "png",
       filename = paste0(plot_dir, "cCRE_overlap_by_exp_mouse_legend_", date, "_peakset=", peakset, ".png"))

ggsave(plot = p4b, width = 20, height = 12, dpi = 300, device = "png",
       filename = paste0(plot_dir, "cCRE_overlap_by_exp_mouse_nolegend_", date, "_peakset=", peakset, ".png"))

ggsave(plot = p4c, width = 5, height = 12, dpi = 300, device = "png",
       filename = paste0(plot_dir, "cCRE_allfrac_mouse_", date, "_peakset=", peakset, ".png"))


# Binned cCRE groups

p5 <- 
  plot_df_mm %>% 
  ggplot(aes(x = Experiment_ID, y = value, fill = Bin_group)) +
  facet_grid(~ Symbol, space = "free", scales = "free") +
  geom_bar(position = "stack", stat = "identity") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.01))) +
  xlab("ChIP-seq experiment") +
  ylab("Proportion of peaks") +
  scale_fill_manual(values = bin_cols, name = "cCRE group") +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 12),
        axis.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15))


ggsave(plot = p5, width = 16, height = 9, dpi = 300, device = "png",
       filename = paste0(plot_dir, "bin_cCRE_overlap_by_exp_mouse_", date, "_peakset=", peakset, ".png"))

p6_list <- 
  lapply(tfs_mm, function(x) {
    p <- 
      plot_df_mm %>% 
      filter(Symbol == x) %>% 
      ggplot(aes(x = Experiment_ID, y = value, fill = Bin_group)) +
      geom_bar(position = "stack", stat = "identity") +
      scale_y_continuous(expand = expansion(mult = c(0, 0.01))) +
      xlab("ChIP-seq experiment") +
      ylab("Proportion of peaks") +
      ggtitle(x) +
      scale_fill_manual(values = bin_cols, name = "cCRE group") +
      theme_classic() +
      theme(axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            plot.title = element_text(size = 12),
            axis.text = element_text(size = 15),
            axis.title = element_text(size = 15),
            axis.title.x = element_blank())
    if (x != tfs_mm[length(tfs_mm)]) {
      p <- p + theme(legend.position = "none")
    }
    return(p)
  })

p6 <- cowplot::plot_grid(plotlist = p6_list)
