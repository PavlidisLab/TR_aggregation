## Explore overlap of ChIP-seq GR objects among ENCODE candidate cis reg elements
## -----------------------------------------------------------------------------

library(tidyverse)
library(GenomicRanges)
library(cowplot)
source("R/setup-01_config.R")
source("R/utils/range_table_functions.R")

pipeline_dir <- paste0(pipeout_dir, "chip/")
plot_dir <- paste0(cplot_dir, "GRanges/")

# GRanges objects
gr_hg <- readRDS(paste0(gr_dir, "human_batch1_grlist_", date, ".RDS"))
gr_mm <- readRDS(paste0(gr_dir, "mouse_batch1_grlist_", date, ".RDS"))

# Region by TF counts
count_list <- readRDS(paste0(scratch_dir, date, "_count_mat_list.RDS"))

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


# Combine cCRE subgroups (promoter, enhancer, other) as factor
# ------------------------------------------------------------------------------


bin_groups <- function(str_vec) {
  
  str_vec <- str_replace_all(str_vec, "\\.|,|-", "_")
  
  group_vec <- vapply(as.character(str_vec), function(x) {
    if (x %in% c("pELS_CTCF_bound",
                 "dELS_CTCF_bound",
                 "dELS",
                 "pELS")) {
      return("Enhancer-like")
    } else if (x %in% c("PLS_CTCF_bound",
                        "PLS")) {
      return("Promoter-like")
    } else if (x %in% c("CTCF_only_CTCF_bound",
                        "DNase_H3K4me3_CTCF_bound",
                        "DNase_H3K4me3")) {
      return("Other")
    } else {
      return("None")
    }
  }, character(1))
  
  group_vec <- factor(
    group_vec, 
    levels = c("None", "Other", "Enhancer-like", "Promoter-like")
  )
}


# Create binned groups
ccre_hg$Bin_group <- as.character(bin_groups(ccre_hg$Group))
ccre_mm$Bin_group <- as.character(bin_groups(ccre_mm$Group))



# levels for all cCRE groups

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



# Fix ChIP-seq peaks to 1bp summit for overlaps with cCRE. 
# This is done to normalize peak sizes across experiments. Some experiments
# have monster peaks that can overlap upwards of 105 elements. Also explored
# adding +/- 150bp to peak summit and see similar results.
# ------------------------------------------------------------------------------


window_size <- 0

gr_rs_hg <- lapply(gr_hg, summit_window, window_size)
gr_rs_mm <- lapply(gr_mm, summit_window, window_size)


# Summarize how many cCRE elements are covered by each peak for each experiment.
# This was motivated by wanting to see how many peaks per data set had no overlap, 
# and to see how many cCRE elements could be covered if peaks were not resized.
# ------------------------------------------------------------------------------


# Returns a named matrix with a row for each element in exp_gr and
# columns tracking the integer count. Each element is the number of peaks from 
# the given experiments that overlapped 0, 1, ... n cCRE overlaps

count_peaks_with_ccre <- function(exp_gr, ccre_gr) {
 
  # for each experiment, get the count of elements that each peak overlaps
  count_list <- lapply(exp_gr, function(x) {
    countOverlaps(x, ccre_gr)
  })
  
  # keep track of the highest count so table can be bound into matrix
  max_count <- max(unlist(lapply(count_list, max)))
  
  table_list <- lapply(count_list, function(x) {
    table(factor(x, levels = c(0:max_count)))  
  })
  
  mat <- as.matrix(do.call(rbind, table_list))
  colnames(mat) <- paste0("Count_", colnames(mat), "_overlapped")
  mat <- cbind(mat, Count_peaks = rowSums(mat))
  
  return(mat)
}


# Convert count matrix to proportion

count_to_prop <- function(count_mat) {
  
  prop_mat <- count_mat
  prop_mat[, 1:ncol(prop_mat) - 1] <- prop_mat[, 1:ncol(prop_mat) - 1] / prop_mat[, "Count_peaks"]
  
  return(prop_mat)
}


# Human: All experiments overlap at least one cCRE. On average, 21% of peaks from
# original set did not overlap with a peak, 28% for resized at 1bp, 22% for 150bp.
# GSE125659_MECP2_Human_MECP2-OE-4x_HISTONE has the least overlap (94% summits 
# no overlap with a cCRE). As expected, MECP2 on average has the most peaks not 
# overlapping a cCRE (81%)


peak_count_hg <- count_peaks_with_ccre(gr_rs_hg, ccre_hg)
peak_count_nors_hg <- count_peaks_with_ccre(gr_hg, ccre_hg)
peak_prop_hg <- count_to_prop(peak_count_hg)
peak_prop_nors_hg <- count_to_prop(peak_count_nors_hg)

any_zero_hg <- which(peak_prop_hg[, "Count_0_overlapped"] == 1)
max_zero_hg <- peak_prop_hg[which.max(peak_prop_hg[, "Count_0_overlapped"]), , drop = FALSE]
max_multi_hg <- peak_prop_hg[which.max(peak_prop_hg[, "Count_2_overlapped"]), , drop = FALSE]

zero_by_tf_hg <- data.frame(peak_prop_hg) %>%
  rownames_to_column(var = "Experiment_ID") %>%
  left_join(meta[, c("Experiment_ID", "Symbol")], by = "Experiment_ID") %>%
  group_by(Symbol) %>%
  summarize(Mean_0 = mean(Count_0_overlapped))


summary(peak_prop_hg[, "Count_0_overlapped"])
summary(peak_prop_nors_hg[, "Count_0_overlapped"])
summary(peak_prop_hg[filter(meta_hg, Symbol != "MECP2")$Experiment_ID, ])
summary(peak_prop_hg[filter(meta_hg, Symbol == "MECP2")$Experiment_ID, ])


# Mouse: GSE139509_Mecp2_Mouse_Ab1-Mecp2-KO_HISTONE has no cCRE overlap - note
# that this is a KO study and has 209 peaks. On average, 40% of peak 
# summits did not overlap a cCRE from original set, 48% with 1bp, 42% 150bp.
# Mecp2 experiments on average 97% no overlap from 1bp set.


peak_count_mm <- count_peaks_with_ccre(gr_rs_mm, ccre_mm)
peak_count_nors_mm <- count_peaks_with_ccre(gr_mm, ccre_mm)
peak_prop_mm <- count_to_prop(peak_count_mm)
peak_prop_nors_mm <- count_to_prop(peak_count_nors_mm)

any_zero_mm <- which(peak_prop_mm[, "Count_0_overlapped"] == 1)
max_zero_mm <- peak_prop_mm[which.max(peak_prop_mm[, "Count_0_overlapped"]), , drop = FALSE]
max_multi_mm <- peak_prop_mm[which.max(peak_prop_mm[, "Count_2_overlapped"]), , drop = FALSE]

zero_by_tf_mm <- data.frame(peak_prop_mm) %>%
  rownames_to_column(var = "Experiment_ID") %>%
  left_join(meta[, c("Experiment_ID", "Symbol")], by = "Experiment_ID") %>%
  group_by(Symbol) %>%
  summarize(Mean_0 = mean(Count_0_overlapped))

summary(peak_prop_mm[, "Count_0_overlapped"])
summary(peak_prop_nors_mm[, "Count_0_overlapped"])
summary(peak_prop_mm[filter(meta_mm, Symbol != "Mecp2")$Experiment_ID, ])
summary(peak_prop_mm[filter(meta_mm, Symbol == "Mecp2")$Experiment_ID, ])


# Find out how many data sets for each TR overlap each cCRE
# ------------------------------------------------------------------------------


# Return a matrix with nrow = length of gr_query and ncol equal to unique TFs in
# meta. Each row corresponds to a peak in gr_query and elements are the 
# count of TF data sets from gr_subject that the peak overlaps

count_ccre_with_tf <- function(gr_query, gr_subject, meta, cores) {
  
  tfs <- unique(meta$Symbol)
  
  # get the indices of the subject hits - requires coercion to GRList
  gr_ol <- findOverlaps(gr_query, GRangesList(gr_subject))
  ol_by_ix <- as(gr_ol, "List")
  
  # get the corresponding TF symbol of the subject and count
  ol_by_tf <- extractList(meta$Symbol, ol_by_ix)
  
  ol_by_count <- mclapply(ol_by_tf, function(x) {
    table(factor(x, levels = tfs))
  }, mc.cores = cores)
  
  # return the matrix of counts with rows named after the range
  count_mat <- as.matrix(do.call(rbind, ol_by_count))
  rownames(count_mat) <- paste0(seqnames(gr_query), ":", start(gr_query), "-", end(gr_query))
  
  return(count_mat)
}


# Human
ccre_by_tf_hg <- count_ccre_with_tf(ccre_hg, gr_rs_hg, meta_hg, cores)

# Mouse
ccre_by_tf_mm <- count_ccre_with_tf(ccre_mm, gr_mm, meta_mm, cores)


# For each experiment, get the breakdown of peak overlap over each cCRE group.
# ------------------------------------------------------------------------------


# gr_query assumed to be GRanges list of peaks, gr_subject to be cCRE GR.
# returns a matrix of nrow = length(gr_query) and ncol = cCRE counts. each
# element of the matrix is the number of cCRE elements for each group that was
# overlapped by each experiment

count_mat_by_group <- function(gr_query, gr_subject, cores) {
 
  count_list <- mclapply(gr_query, function(x) {
    
    # get the indices of the subject hits
    gr_ol <- findOverlaps(x, gr_subject)
    ol_by_ix <- gr_ol@from
    
    # get the corresponding cCRE group of the subject and count
    ol_by_group <- unlist(extractList(gr_subject$Group, ol_by_ix))
    ol_count <- table(factor(ol_by_group, levels = unique(gr_subject$Group)))
  
    }, mc.cores = cores)
  
  count_mat <- as.matrix(do.call(rbind, count_list))
  
  return(count_mat)
}


# Human:

group_count_hg <- count_mat_by_group(gr_rs_hg, ccre_hg, cores)

group_count_hg <- 
  cbind(None = peak_count_hg[, "Count_0_overlapped"], group_count_hg)

group_prop_hg <- group_count_hg / rowSums(group_count_hg)

group_by_tf_hg <- data.frame(group_prop_hg) %>% 
  rownames_to_column(var = "Experiment_ID") %>% 
  left_join(meta[, c("Experiment_ID", "Symbol")], by = "Experiment_ID") %>% 
  group_by(Symbol) %>% 
  summarise(across(where(is.double), list(mean)))


# Using the ENCODE experiments targeting RUNX1 (but diff antibodies and seq
# tech) to compare group breakdown - find quite similar
group_prop_hg[c("GSE91747_RUNX1_Human_K562-ENCODE-Ab1", "GSE96253_RUNX1_Human_K562-ENCODE-Ab2"), ]

# Using RUNX1 Kasumi-1 (most represented cell type) to show cCRE group spread
summary(group_prop_hg[filter(meta, Symbol == "RUNX1")$Experiment_ID, ])
group_prop_hg[filter(meta, Symbol == "RUNX1")$Experiment_ID, ]

kasumi <- filter(meta, Symbol == "RUNX1" & Cell_Type == "Kasumi-1")$Experiment_ID
non_kasumi <- filter(meta, Symbol == "RUNX1" & Cell_Type != "Kasumi-1")$Experiment_ID


pxa <- 
  data.frame(group_prop_hg[kasumi, ]) %>% 
  rownames_to_column(var = "Experiment_ID") %>% 
  reshape::melt(id = "Experiment_ID") %>% 
  ggplot(., aes(x = variable, y = value)) +
  geom_boxplot() +
  ylab("Proportion of overlap") +
  ggtitle("RUNX1 Kasumi-1 experiments") +
  ylim(c(0, 0.6)) +
  theme_classic() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 20),
        axis.text = element_text(size = 20),
        axis.text.x = element_text(size = 15, angle = 90, vjust = 0.5, hjust = 1))

pxb <- 
  data.frame(group_prop_hg[non_kasumi, ]) %>% 
  rownames_to_column(var = "Experiment_ID") %>% 
  reshape::melt(id = "Experiment_ID") %>% 
  ggplot(., aes(x = variable, y = value)) +
  geom_boxplot() +
  ylab("Proportion of overlap") +
  ggtitle("RUNX1 non-Kasumi-1 experiments") +
  ylim(c(0, 0.6)) +
  theme_classic() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 20),
        axis.text = element_text(size = 20),
        axis.text.x = element_text(size = 15, angle = 90, vjust = 0.5, hjust = 1))


plot_grid(pxa, pxb)


summary(group_prop_hg[filter(meta, Symbol == "RUNX1" & Cell_Type != "Kasumi-1")$Experiment_ID, ])
boxplot(group_prop_hg[filter(meta, Symbol == "RUNX1" & Cell_Type != "Kasumi-1")$Experiment_ID, ])

summary(group_prop_hg[filter(meta, Symbol == "RUNX1" & Cell_Type == "Kasumi-1")$Experiment_ID, ])
boxplot(group_prop_hg[filter(meta, Symbol == "RUNX1" & Cell_Type == "Kasumi-1")$Experiment_ID, ])


# Mouse:

group_count_mm <- count_mat_by_group(gr_rs_mm, ccre_mm, cores)

group_count_mm <- 
  cbind(None = peak_count_mm[, "Count_0_overlapped"], group_count_mm)

group_prop_mm <- group_count_mm / rowSums(group_count_mm)
  
group_by_tf_mm <- data.frame(group_prop_mm) %>% 
  rownames_to_column(var = "Experiment_ID") %>% 
  left_join(meta[, c("Experiment_ID", "Symbol")], by = "Experiment_ID") %>% 
  group_by(Symbol) %>% 
  summarise(across(where(is.double), list(mean)))


# Explore region x TF matrix and look at cCRE status of top bound regions
# ------------------------------------------------------------------------------


# Sort count_mat by tf column while minimizing the sum of the rest of the cols.
# Return the df of in vs out counts, as well as the ordered count mat

top_tf_count <- function(count_mat, tf) {
  
  count_df <- data.frame(
    Intra_TF = count_mat[, tf],
    Inter_TF = rowSums(count_mat[, setdiff(colnames(count_mat), tf)]),
    row.names = NULL)  # otherwise names may get preserved - need int for sort
  
  count_df <- count_df[order(count_df$Intra_TF, -count_df$Inter_TF, decreasing = TRUE), ]
  colnames(count_df) <- c("Intra_TF", "Inter_TF")
  
  return(list(Count_mat = count_mat[as.integer(rownames(count_df)),],
              Count_df = count_df))
}


# Sort matrix of proportions by TF with the highest relative overlap

top_tf_prop <- function(count_mat, tf, meta) {
  
  tf_count <- c(In = sum(meta$Symbol == tf), Out = sum(meta$Symbol != tf))
  
  count_df <- data.frame(
    Intra_TF = (count_mat[, tf, drop = FALSE] / tf_count["In"]),
    Inter_TF = (rowSums(count_mat[, setdiff(colnames(count_mat), tf)]) / tf_count["Out"]),
    row.names = rownames(count_mat))
  
  colnames(count_df) <- c("Intra_TF", "Inter_TF")
  
  count_df[order(count_df$Intra_TF, -count_df$Inter_TF, decreasing = TRUE), ]
  
}


# Given regions as character vector of form "9:38022988-38023319" find overlaps
# with provded cCRE object and return the associated cCRE group. If no overlap, 
# return "None", if multiple exist, return "Mixed"


get_ccre_group <- function(regions, ccre_gr) {
  
  ol <- findOverlaps(GRanges(regions), ccre_gr)

  ccre_group <- vapply(1:length(regions), function(x) {
    
    group <- unique(ccre_gr$Bin_group[ol[ol@from == x]@to])
    
    if (length(group) == 0) {
      group <- "None"
    } else if (length(group) > 1) {
      group <- "Mixed"
    }
    return(group)
    
  }, FUN.VALUE = character(1))

  ccre_group <- factor(
    ccre_group,
    levels = c("None", "Other", "Mixed", "Enhancer-like", "Promoter-like"))
  
  return(ccre_group)
}



# Human - focused on ASCL1. Specifically find region located in intron of SHB
# (9, 38022988:38023319) where every ASCL1 experiment had a binding event, while
# only 2 other experiments (TCF4) were bound. Correpsonds to an enhancer cCRE


tf <- "ASCL1"

# Order count matrix by regions that have most common binding in TF relative
# to rest of the TFs
top_mat_hg <- top_tf_count(count_list$Human$Reduced_resize, tf)$Count_mat

# Regions with max binding
top_regions_hg <- top_mat_hg[top_mat_hg[, tf] == max(top_mat_hg[, tf]), ]
nrow(top_regions_hg)

# cCRE status of these top regions
top_ccre_hg <- get_ccre_group(rownames(top_regions_hg), ccre_hg)


# Mouse - focused on Neurod1. 10:49985040-49985394 highly conserved, no ccRE
# status. Closest gene is Gm48131


tf <- "Neurod1"
top_mat_mm <- top_tf_count(count_list$Mouse$Reduced_resize, tf)$Count_mat
top_regions_mm <- top_mat_mm[top_mat_mm[, tf] == max(top_mat_mm[, tf]), ]
nrow(top_regions_mm)
top_ccre_mm <- get_ccre_group(rownames(top_regions_mm), ccre_mm)



# Plots
# ------------------------------------------------------------------------------


# Colours for different cCRE groups
all_cols <- c('#d9d9d9', '#006d2c','#74c476','#7bccc4','#9ecae1','#4292c6', '#08519c', '#ff7f00','#cab2d6', '#fdbf6f')
names(all_cols) <- levels(all_ccre)

# Colours for binned group (+mixed for regions bordering 2 groups)
bin_cols <- list(None = '#d9d9d9',
                 `Promoter-like` = '#33a02c',
                 `Enhancer-like` = '#1f78b4',
                 Other = '#6a3d9a', 
                 Mixed = '#fdbf6f')


# boxplot of the proportion of peaks per experiment that did not overlap a cCRE,
# grouped by symbol
boxplot(peak_prop_hg[, "Count_0_overlapped"] ~ meta_hg$Symbol)
boxplot(peak_prop_mm[, "Count_0_overlapped"] ~ meta_mm$Symbol)


# for looking at distn of cCRE groups by symbol
boxplot(group_prop_hg[, "pELS"] ~ meta_hg$Symbol)
boxplot(group_prop_mm[, "pELS"] ~ meta_mm$Symbol)


# Stacked barchart of cCRE membership per experiment. Include total proportion
# of cCRE groups as a comparison. 


# data frame of the proportion of cCREs per group from the whole set

get_ccre_df <- function(gr_ccre) {
  
  ccre_df <- data.frame(Group = gr_ccre$Group, stringsAsFactors = FALSE) %>% 
    group_by(Group) %>% 
    summarize(Proportion = n()/length(gr_ccre$Group)) %>% 
    ungroup() %>% 
    mutate(Plot_var = "1")  # Var to supply plot call
  ccre_df$Group <- factor(ccre_df$Group, levels = levels(all_ccre))
  
  return(ccre_df)
}


# dataframe of each experiment's proportion of overlap with cCRE groups


get_prop_df <- function(prop_mat, meta, ccre_gr, group_levels) {
  
  prop_df <- data.frame(prop_mat) %>% 
    rownames_to_column(var = "Experiment_ID") %>% 
    left_join(meta[, c("Symbol", "Experiment_ID", "N_peaks")], 
              by = "Experiment_ID") %>%
    pivot_longer(cols = c(unique(ccre_gr$Group), "None"), 
                 names_to = "Group",
                 values_to = "Proportion") %>% 
    mutate(
      Group = factor(Group, levels = group_levels),
      Bin_group = bin_groups(Group),
      Experiment_ID = factor(Experiment_ID, levels = unique(Experiment_ID))
    )
  
  return(prop_df)
}


exp_barchart <- function(plot_df, group_fill, colours) {
  
  plot_df %>% 
    ggplot(aes(x = Experiment_ID, y = Proportion, fill = !!sym(group_fill))) +
    facet_grid(~ Symbol, space = "free", scales = "free") +
    geom_bar(position = "stack", stat = "identity", width = 1) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.01))) +
    xlab("ChIP-seq experiment") +
    ylab("Proportion of peaks") +
    scale_fill_manual(values = colours, name = "cCRE group") +
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
  
}


ccre_barchart <- function(ccre_df, colours) {
  
  ccre_df %>% 
    ggplot(aes(x = Plot_var, y = Proportion, fill = Group)) +
    geom_bar(position = "stack", stat = "identity", width = 0.6) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.01))) +
    ylab("Proportion of all cCREs") +
    scale_fill_manual(values = colours, name = "cCRE group") +
    theme_classic() +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.x = element_blank(),
          axis.text = element_text(size = 25),
          axis.title = element_text(size = 25),
          legend.position = "none",
          plot.margin = margin(10, 10, 10, 10))  # add padding for removing legend
}



# Human

ccre_df_hg <- get_ccre_df(ccre_hg)
prop_df_hg <- get_prop_df(group_prop_hg, meta_hg, ccre_hg, levels(all_ccre))

# +legend
p1a <- exp_barchart(prop_df_hg, "Group", all_cols)

# -legend and strip names
p1b <- p1a + theme(legend.position = "none",
                   strip.text = element_blank(),
                   axis.title.x = element_blank())

# all cCREs 
p1c <- ccre_barchart(ccre_df_hg, all_cols)

# bin groups
p1d <- exp_barchart(prop_df_hg, "Bin_group", bin_cols)

ggsave(plot = p1a, width = 20, height = 12, dpi = 300, device = "png",
  filename = paste0(plot_dir, "cCRE_overlap_by_exp_human_legend_", date, ".png"))

ggsave(plot = p1b, width = 20, height = 12, dpi = 300, device = "png",
       filename = paste0(plot_dir, "cCRE_overlap_by_exp_human_nolegend_", date, ".png"))

ggsave(plot = p1c, width = 5, height = 12, dpi = 300, device = "png",
       filename = paste0(plot_dir, "cCRE_allprop_human_", date, ".png"))

ggsave(plot = p1d, width = 16, height = 9, dpi = 300, device = "png",
       filename = paste0(plot_dir, "bin_cCRE_overlap_by_exp_human_", date, ".png"))

# Mouse

ccre_df_mm <- get_ccre_df(ccre_mm)
prop_df_mm <- get_prop_df(group_prop_mm, meta_mm, ccre_mm, levels(all_ccre))

# +legend
p2a <- exp_barchart(prop_df_mm, "Group", all_cols)

# -legend and strip names
p2b <- p2a + theme(legend.position = "none",
                   strip.text = element_blank(),
                   axis.title.x = element_blank())

# all cCREs 
p2c <- ccre_barchart(ccre_df_mm, all_cols)

# bin groups
p2d <- exp_barchart(prop_df_mm, "Bin_group", bin_cols)

ggsave(plot = p2a, width = 20, height = 12, dpi = 300, device = "png",
       filename = paste0(plot_dir, "cCRE_overlap_by_exp_mouse_legend_", date, ".png"))

ggsave(plot = p2b, width = 20, height = 12, dpi = 300, device = "png",
       filename = paste0(plot_dir, "cCRE_overlap_by_exp_mouse_nolegend_", date, ".png"))

ggsave(plot = p2c, width = 5, height = 12, dpi = 300, device = "png",
       filename = paste0(plot_dir, "cCRE_allprop_mouse_", date, ".png"))

ggsave(plot = p2d, width = 16, height = 9, dpi = 300, device = "png",
       filename = paste0(plot_dir, "bin_cCRE_overlap_by_exp_mouse_", date, ".png"))


# Scatter of high region overlap for in vs out, coloured by grouped cCRE status


plot_scatter <- function(plot_df, tf) {
  
  ggplot(plot_df, aes(x = Intra_TF, y = Inter_TF, fill = Group)) +
    geom_jitter(size = 3, shape = 21, height = 0.008, width = 0.01) +
    theme_classic() +
    ylab("Proportion overlap inter-TF") +
    xlab("Proportion overlap intra-TF") +
    ggtitle(tf) +
    scale_fill_manual(values = bin_cols) +
    theme(axis.text = element_text(size = 25),
          axis.title = element_text(size = 25),
          plot.title = element_text(size = 25),
          legend.text = element_text(size = 15),
          legend.title = element_text(size = 15),
          legend.position = c(.85, .9))
}



plot_all_scatter <- function(count_mat, meta, cutoff = 0.5) {
  
  tfs <- unique(meta$Symbol)
  
  plot_l <- lapply(tfs, function(tf) {
    
    top_df <-
      top_tf_prop(count_mat, tf, meta) %>%
      rownames_to_column(var = "Region") %>%
      filter(Intra_TF > cutoff)
    
    # in case there no regions at the cutoff
    
    if (nrow(top_df) == 0) {
      top_df <-
        top_tf_prop(count_mat, tf, meta_mm) %>%
        rownames_to_column(var = "Region") %>%
        filter(Intra_TF > 0.2)
    }
    
    top_df$Group <- get_ccre_group(top_df$Region, ccre_hg)
    plot_scatter(top_df, tf)
    
  })
  names(plot_l) <- tfs
  
  return(plot_l)
}


# Human

pl_hg <- plot_all_scatter(count_list$Human$Reduced_resize, meta_hg)

pdf(paste0(plot_dir, date, "_human_freq_bound_by_ccre.pdf"))
invisible(lapply(pl_hg, print))
graphics.off()


ggsave(pl_hg$ASCL1, dpi = 300, device = "png", height = 8, width = 10,
                filename = paste0(plot_dir, "ASCL1_freq_bound.png"))

# Mouse

pl_mm <- plot_all_scatter(count_list$Mouse$Reduced_resize, meta_mm)

pdf(paste0(plot_dir, date, "_mouse_freq_bound_by_ccre.pdf"))
invisible(lapply(pl_mm, print))
graphics.off()
