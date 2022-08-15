## Explore overlap of ChIP-seq GR objects among all and TF-specific experiments
## TODO: clean up top/count/frac mat helpers
## TODO: formalize overlap of frequent with cCREs (current plot scheme highly inefficient)
## TODO: save logic of region tables
## -----------------------------------------------------------------------------

library(tidyverse)
library(GenomicRanges)
library(ggbio)
source("~/regnetR/R/utils/range_table_functions.R")

date <- "Apr2022"
peakset <- "idr"
plot_dir <- "~/Plots/Chipseq/GRanges/"

# GRanges objects
gr_hg <- readRDS(paste0("~/Data/Annotated_objects/GRanges/human_batch1_grlist_peakset=", peakset, "_", date, ".RDS"))
gr_mm <- readRDS(paste0("~/Data/Annotated_objects/GRanges/mouse_batch1_grlist_peakset=", peakset, "_", date, ".RDS"))

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


# Look at summary of peak sizes across all experiments
# ------------------------------------------------------------------------------


# Human: GSE69394_NEUROD1_Human_H82 highest max peak size (123kb - 8:127638060-127761390)
# which is stretch 5' of MYC that contains non protein coding CASC11
# GSE58428_RUNX1_Human_VCaP highest mean peak size (1.67kb) - has high RSC

peak_size_hg <- data.frame(
  do.call(rbind, lapply(gr_hg, function(x) summary(width(x))))
) %>%
  rownames_to_column(var = "Experiment_ID") %>%
  left_join(meta[, c("Symbol", "Experiment_ID")])
rownames(peak_size_hg) <- peak_size_hg$Experiment_ID


max_hg <- peak_size_hg[which.max(peak_size_hg$Max.),]
max_peak_hg <- gr_hg[rownames(max_hg)][width(gr_hg[rownames(max_hg)]) == max_hg$Max.]
max_mean_hg <- peak_size_hg[which.max(peak_size_hg$Mean),]
meta[meta$Experiment_ID %in% c(rownames(max_mean_hg), names(max_peak_hg)), ]


# Mouse: Largest max peak size GSE90893_Runx1_Mouse_MEF-Control
# (73kb - 1:71582211-71655087) 
# Largest mean peak size belong to GSE126375_Runx1_Mouse_Pax5-KO (1.7kb)


peak_size_mm <- data.frame(
  do.call(rbind, lapply(gr_mm, function(x) summary(width(x))))
) %>%
  rownames_to_column(var = "Experiment_ID") %>%
  left_join(meta[, c("Symbol", "Experiment_ID")])
rownames(peak_size_mm) <- peak_size_mm$Experiment_ID


max_mm <- peak_size_mm[which.max(peak_size_mm$Max.),]
max_peak_mm <- gr_mm[rownames(max_mm)][width(gr_mm[rownames(max_mm)]) == max_mm$Max.]
max_mean_mm <- peak_size_mm[which.max(peak_size_mm$Mean),]
meta[meta$Experiment_ID %in% c(rownames(max_mean_mm), names(max_peak_mm)), ]


# Look at how many peaks overlap within data sets before fixing size.
# Find no intra-peak overlaps (which is expected based on pipeline output).
# Closest peaks across all datasets is 33bp for both human and mouse, all in
# Runx1 experiments
# ------------------------------------------------------------------------------


any_overlaps_hg <- vapply(gr_hg, function(gr) {
  return (length(gr) == length(reduce(gr)))
}, FUN.VALUE = logical(1))

stopifnot(all(any_overlaps_hg))


any_overlaps_mm <- vapply(gr_mm, function(gr) {
  return (length(gr) == length(reduce(gr)))
}, FUN.VALUE = logical(1))

stopifnot(all(any_overlaps_mm))


nearest_dist_hg <- unlist(mclapply(gr_hg, function(x) {
  min(distanceToNearest(x)@elementMetadata$distance)
}, mc.cores = 8))


nearest_dist_mm <- unlist(mclapply(gr_mm, function(x) {
  min(distanceToNearest(x)@elementMetadata$distance)
}, mc.cores = 8))


nearest_dist_hg[nearest_dist_hg == min(nearest_dist_hg)]
nearest_dist_mm[nearest_dist_mm == min(nearest_dist_mm)]


# Because peaks have variable sizes, fix to summit position and add padding to
# make all same size for overlap (1bp summit too stringent for overlap)
# ------------------------------------------------------------------------------


window_size <- 150

summit_window <- function(gr, window_size) {
  # Fix range to the summit position and pad with a fixed window
  start(gr) <- end(gr) <- start(gr) + gr$Summit_from_start
  gr <- gr + window_size
  return(gr)
}


# re-size
gr_rs_hg <- GRangesList(lapply(gr_hg, summit_window, window_size))
gr_rs_mm <- GRangesList(lapply(gr_mm, summit_window, window_size))


# Note that re-sizing 100bp+ results in overlapping peaks that were not
# overlapping in their original forms. In the "reduced" set these re-sized
# peaks are then merged, despite the fact that they originally represented two
# summits from a single data set. Following was explored using
# "GSE61197_ASCL1_Human_H1775-ASCL1-high"


plot_gr_summit <- function(gr) {
  # Make a standard Grange plot and add a line corresponding to the summit
  ggbio::autoplot(gr) + 
    geom_vline(xintercept = gr$Summit) +
    theme_clear()
}


# stash an unedited peak and directly add summit for demo plot
test_gr_before <- gr_hg[["GSE61197_ASCL1_Human_H1775-ASCL1-high"]]
test_gr_before$Summit <- start(test_gr_before) + test_gr_before$Summit_from_start
length(test_gr_before) == length(reduce(test_gr_before)) # expect T


# stash an edited peak
test_gr_after <- gr_rs_hg[["GSE61197_ASCL1_Human_H1775-ASCL1-high"]]
test_gr_after$Summit <- test_gr_before$Summit
length(test_gr_after) == length(reduce(test_gr_after))  # expect F


# split the re-sized peak and find the indices corresponding to self overlaps
self_overlap <- findOverlaps(
  split(test_gr_after, rep(1:length(test_gr_after), each=1)), test_gr_after)

dup_ix <- unique(self_overlap@to[duplicated(self_overlap@to)])

# get the closest distances for peaks that overlap in the re-sized set
self_dups_before <- test_gr_before[dup_ix]
self_dups_after <- test_gr_after[dup_ix]

dist_before <- distanceToNearest(self_dups_before)@elementMetadata$distance
dist_after <- distanceToNearest(self_dups_after)@elementMetadata$distance

# plot demonstrates overlapping ranges (plus summit position) before/after re-size
p_before <- plot_gr_summit(self_dups_before[1:2])
p_after <- plot_gr_summit(self_dups_after[1:2])

# plot that shows example of a large peak that originally does not overlap
# but does so after re-sizing
max_ix <- dup_ix[which.max(width(self_dups_before))]
p_max_before <- plot_gr_summit(test_gr_before[(max_ix-1):max_ix])
p_max_after <- plot_gr_summit(test_gr_after[(max_ix-1):max_ix])


# Describe the full and reduced (non-overlapping) set of peaks for all 
# experiments with and without the re-size.
# Note that when peaks have been re-sized, the reduced set can create large 
# peaks that were not-overlapping in their original forms (summits near the 
# edge close to another peak, and the window size padding causes overlap)
# ------------------------------------------------------------------------------


# Human: All and re-size ~3.5m peaks across 141 experiments which drops to ~1.22m
# in the re-sized and reduced set (-65%) and ~893k in the original reduced set (-75%)

# The longest reduced peak in the original reduced set is 138kb 
# (8: 127630738-127769153) which overlaps peaks in 67 data sets.
# The longest in the re-sized reduced set is 3.4kb (19: 16585409-16588814)
# in 40 data sets


all_hg <- list(
  All = unlist(gr_hg),
  Resize = unlist(gr_rs_hg),
  Reduced = reduce(unlist(gr_hg)),
  Reduced_resize = reduce(unlist(gr_rs_hg))
)

count_hg <- sapply(all_hg, length)
change_hg <- (count_hg["Reduced"] - count_hg["All"]) / count_hg["All"]
change_rs_hg <- (count_hg["Reduced_resize"] - count_hg["Resize"]) / count_hg["Resize"]

maxpeak_hg <- all_hg$Reduced[which.max(width(all_hg$Reduced))]
maxpeak_rs_hg <- all_hg$Reduced_resize[which.max(width(all_hg$Reduced_resize))]

table(meta_hg$Symbol[unique(findOverlaps(maxpeak_hg, gr_hg)@to)])
table(meta_hg$Symbol[unique(findOverlaps(maxpeak_rs_hg, gr_rs_hg)@to)])


# Mouse: All and re-size ~1.95m peaks across 136 experiments which drops to ~686k
# in the re-sized and reduced set (-65%) and ~588k in the original reduced set (-70%)

# The longest reduced peak in the original reduced set is 73kb 
# (1: 71582211-71655460 ) which overlaps peaks in 36 data sets.
# The longest in the re-sized reduced set is 3.3kb (7: 127024978-127028247)
# in 57 data sets


all_mm <- list(
  All = unlist(gr_mm),
  Resize = unlist(gr_rs_mm),
  Reduced = reduce(unlist(gr_mm)),
  Reduced_resize = reduce(unlist(gr_rs_mm))
)

count_mm <- sapply(all_mm, length)
change_mm <- (count_mm["Reduced"] - count_mm["All"]) / count_mm["All"]
change_rs_mm <- (count_mm["Reduced_resize"] - count_mm["Resize"]) / count_mm["Resize"]

maxpeak_mm <- all_mm$Reduced[which.max(width(all_mm$Reduced))]
maxpeak_rs_mm <- all_mm$Reduced_resize[which.max(width(all_mm$Reduced_resize))]

table(meta_mm$Symbol[unique(findOverlaps(maxpeak_mm, gr_mm)@to)])
table(meta_mm$Symbol[unique(findOverlaps(maxpeak_rs_mm, gr_rs_mm)@to)])


# Look at the sites that are most commonly overlapped across the all set
# ------------------------------------------------------------------------------


all_overlap <- function(all_gr, exp_gr, meta) {
  # all_gr is a single gr, exp_gr is a gr list of experiments
  
  ol <- findOverlaps(all_gr, exp_gr)
  counts <- table(ol@from)
  most_ol <- all_gr[which(counts == max(counts))]
  most_ol_tf <- table(meta$Symbol[unique(findOverlaps(most_ol, exp_gr)@to)])
  width_most_ol <- width(range(most_ol))
  
  return(list(Overlap = ol,
              Count_table = counts,
              Most_overlapped = most_ol,
              Count_exp_most_overlapped = sum(most_ol_tf),
              Most_overlapped_TF = most_ol_tf,
              Most_overlapped_width = width_most_ol))
}


# Human: 
# Most overlapping regions 16:85604925-85645912, 11: 65496125-65508279, and
# 19: 44753325-44753587

all_ol_hg <- list(
  All = all_overlap(all_hg$All, gr_hg, meta_hg),
  Resize = all_overlap(all_hg$Resize, gr_rs_hg, meta_hg),
  Reduced = all_overlap(all_hg$Reduced, gr_hg, meta_hg),
  Reduced_resize = all_overlap(all_hg$Reduced_resize, gr_rs_hg, meta_hg)
)

lapply(all_ol_hg, `[[`, "Most_overlapped")
lapply(all_ol_hg, `[[`, "Count_exp_most_overlapped")
lapply(all_ol_hg, `[[`, "Most_overlapped_width")


# Mouse: 
# Most overlapping regions 3:88509241-88512113

all_ol_mm <- list(
  All = all_overlap(all_mm$All, gr_mm, meta_mm),
  Resize = all_overlap(all_mm$Resize, gr_rs_mm, meta_mm),
  Reduced = all_overlap(all_mm$Reduced, gr_mm, meta_mm),
  Reduced_resize = all_overlap(all_mm$Reduced_resize, gr_rs_mm, meta_mm)
)

lapply(all_ol_mm, `[[`, "Most_overlapped")
lapply(all_ol_mm, `[[`, "Count_exp_most_overlapped")
lapply(all_ol_mm, `[[`, "Most_overlapped_width")


# Create TF-specific all and reduced GR objects
# ------------------------------------------------------------------------------


gr_by_tf <- function(tfs, exp_gr, meta) {
  # Grlist of experiments grouped into list by TF
  tf_gr <- lapply(tfs, function(tf) {
    runs <- meta[meta$Symbol == tf, "Experiment_ID"]
    unlist(exp_gr[runs])
  })
  names(tf_gr) <- tfs
  return(tf_gr)
}


tf_gr_hg <- list(
  All = gr_by_tf(tfs_hg, gr_hg, meta_hg),
  Resize = gr_by_tf(tfs_hg, gr_rs_hg, meta_hg))
tf_gr_hg$Reduced <- lapply(tf_gr_hg$All, reduce)
tf_gr_hg$Reduced_resize <- lapply(tf_gr_hg$Resize, reduce)


tf_gr_mm <- list(
  All = gr_by_tf(tfs_mm, gr_mm, meta_mm),
  Resize = gr_by_tf(tfs_mm, gr_rs_mm, meta_mm))
tf_gr_mm$Reduced <- lapply(tf_gr_mm$All, reduce)
tf_gr_mm$Reduced_resize <- lapply(tf_gr_mm$Resize, reduce)



# Decrease in size when reduced. a small % change suggests the TF experiments
# were different, since the reduced set retains unique peaks. a large % change
# suggests the experiments were similar - overlapping peaks get reduced
# ------------------------------------------------------------------------------


# Human: RUNX1 most reduced size, HES1 least

tf_counts_hg <- sapply(tf_gr_hg, function(peakset) sapply(peakset, function(tf) length(tf)))
change_hg <- (tf_counts_hg[, "Reduced"] - tf_counts_hg[, "All"]) / tf_counts_hg[, "All"]
change_rs_hg <- (tf_counts_hg[, "Reduced_resize"] - tf_counts_hg[, "Resize"]) / tf_counts_hg[, "Resize"]

# Mouse: Runx1 most reduced size, Mecp2 least

tf_counts_mm <- sapply(tf_gr_mm, function(peakset) sapply(peakset, function(tf) length(tf)))
change_mm <- (tf_counts_mm[, "Reduced"] - tf_counts_mm[, "All"]) / tf_counts_mm[, "All"]
change_rs_mm <- (tf_counts_mm[, "Reduced_resize"] - tf_counts_mm[, "Resize"]) / tf_counts_mm[, "Resize"]


# Get the count of TF data sets that each range overlaps. Arrange by counts 
# for each TF while minimizing sum of counts for the other TFs. Suggestive of
# TF-specific binding sites
# NOTE: exploration does on reduced, re-sized
# ------------------------------------------------------------------------------


get_count_mat <- function(gr_query, gr_subject, meta, tfs) {
  # Return a matrix with nrow = length of gr_query and ncol equal to length
  # of TFs. each row corresponds to a peak in gr_query and elements are the 
  # number of TF data sets from gr_subject that the peak overlaps
  
  # get the indices of the subject hits
  gr_ol <- findOverlaps(gr_query, gr_subject)
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


count_to_frac <- function(count_mat, meta) {
  
  tf_count <- filter(meta) %>% 
    count(Symbol) %>% 
    arrange(Symbol, colnames(count_mat)) %>% 
    pull(n)
  
  t(t(count_mat) / tf_count)
}


top_tf <- function(count_mat, tf) {
  
  # Sort count_mat by tf column while minimizing the sum of the rest of the cols.
  # Return the df of in vs out counts, as well as the ordered count mat
  
  count_df <- data.frame(
    Intra_TF = count_mat[, tf],
    Inter_TF <- rowSums(count_mat[, setdiff(colnames(count_mat), tf)]),
    row.names = NULL)  # otherwise names may get preserved - need int for sort
  
  count_df <- 
    count_df[order(count_df$Intra_TF, -count_df$Inter_TF, decreasing = TRUE), ]
  colnames(count_df) <- c("Intra_TF", "Inter_TF")
  
  return(list(Count_mat = count_mat[as.integer(rownames(count_df)),],
              Count_df = count_df))
}


top_tf_frac <- function(count_mat, tf, meta) {
  
  tf_count <- c(In = sum(meta$Symbol == tf), Out = sum(meta$Symbol != tf))
  
  count_df <- data.frame(
    Intra_TF = (count_mat[, tf, drop = FALSE] / tf_count["In"]),
    Inter_TF <- (rowSums(count_mat[, setdiff(colnames(count_mat), tf)]) / tf_count["Out"]),
    row.names = rownames(count_mat))
  
  colnames(count_df) <- c("Intra_TF", "Inter_TF")
  
  count_df[order(count_df$Intra_TF, -count_df$Inter_TF, decreasing = TRUE), ]
  
}


# Human: 

count_mat_hg <- list(
  All = get_count_mat(all_hg$All, gr_hg, meta_hg, tfs_hg),
  Resize = get_count_mat(all_hg$Resize, gr_rs_hg, meta_hg, tfs_hg),
  Reduced = get_count_mat(all_hg$Reduced, gr_hg, meta_hg, tfs_hg),
  Reduced_resize = get_count_mat(all_hg$Reduced_resize, gr_rs_hg, meta_hg, tfs_hg)
)


tf <- "ASCL1"

top_mat <- top_tf(count_mat_hg$Reduced_resize, tf)$Count_mat
sum(top_mat[, tf] == max(top_mat[, tf]))


# names(gr_hg)[findOverlaps(GRanges(9, 38022988:38023319), gr_rs_hg)@to]


# Mouse: 


count_mat_mm <- list(
  All = get_count_mat(all_mm$All, gr_mm, meta_mm, tfs_mm),
  Resize = get_count_mat(all_mm$Resize, gr_rs_mm, meta_mm, tfs_mm),
  Reduced = get_count_mat(all_mm$Reduced, gr_mm, meta_mm, tfs_mm),
  Reduced_resize = get_count_mat(all_mm$Reduced_resize, gr_rs_mm, meta_mm, tfs_mm)
)


tf <- "Neurod1"

top_mat <- top_tf(count_mat_mm$Reduced_resize, tf)$Count_mat


# Save out
saveRDS(list(Human = count_mat_hg, Mouse = count_mat_mm),
        file = paste0("~/scratch/R_objects/", date, "_count_mat_list.RDS"))


# Plots
# ------------------------------------------------------------------------------


### DEMO overlap with cCRE

# cCRE tables -> GR objects
ccre_hg <- read.delim("~/Data/Chromosome_info/cCREs_V3_hg38.bed", stringsAsFactors = FALSE)
ccre_hg <- makeGRangesFromDataFrame(ccre_hg, keep.extra.columns = TRUE)
ccre_hg$Group <- str_replace_all(ccre_hg$Group, ",|-", "_")

ccre_mm <- read.delim("~/Data/Chromosome_info/cCREs_V3_mm10.bed", stringsAsFactors = FALSE)
ccre_mm <- makeGRangesFromDataFrame(ccre_mm, keep.extra.columns = TRUE)
ccre_mm$Group <- str_replace_all(ccre_mm$Group, ",|-", "_")


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


ccre_hg$Bin_group <- as.character(bin_groups(ccre_hg$Group))
ccre_mm$Bin_group <- as.character(bin_groups(ccre_mm$Group))


get_ccre_group <- function(top_df, ccre_gr) {
  # Given top_df with a column of regions, find overlaps with provded cCRE
  # object and return the associated cCRE group. If no overlap, return "None", 
  # if multiple exist, return "Mixed"
  
  top_gr <- GRanges(top_df$Region)
  ol <- findOverlaps(top_gr, ccre_gr)
  # out_vec <- rep("None", nrow(top_df))
  
  ccre_group <- unlist(lapply(1:nrow(top_df), function(x) {
    group <- unique(ccre_gr$Bin_group[ol[ol@from == x]@to])
    if (length(group) == 0) {
      group <- "None"
    } else if (length(group) > 1) {
      group <- "Mixed"
    }
    return (group)
  }))
  
  ccre_group <- factor(
    ccre_group,
    levels = c("None", "Other", "Mixed", "Enhancer-like", "Promoter-like"))
  
  return(ccre_group)

}


# ccre_hg[findOverlaps(GRanges("9:137441196-137441859"), ccre_hg)@to]


bin_cols <- c('#d9d9d9','#6a3d9a','#fdbf6f','#1f78b4', "#33a02c")


# scatter of high region overlap for in vs out


plist_hg <- lapply(tfs_hg, function(tf) {
  
  top_df <- 
    top_tf_frac(count_mat_hg$Reduced_resize, tf, meta_hg) %>% 
    rownames_to_column(var = "Region") %>% 
    filter(Intra_TF > 0.5)
  
  top_df$Group <- get_ccre_group(top_df, ccre_hg)
  
  p <- 
    ggplot(top_df, aes(x = Intra_TF, y = Inter_TF, fill = Group)) +
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
  
})
names(plist_hg) <- tfs_hg


pdf(paste0(plot_dir, date, "_human_freq_bound_by_ccre.pdf"))
invisible(lapply(plist_hg, print))
graphics.off()


ggplot2::ggsave(plist_hg$ASCL1, dpi = 300, device = "png", height = 8, width = 10,
       filename = paste0(plot_dir, "ASCL1_freq_bound.png"))


plist_mm <- lapply(tfs_mm, function(tf) {
  
  top_df <- 
    top_tf_frac(count_mat_mm$Reduced_resize, tf, meta_mm) %>% 
    rownames_to_column(var = "Region") %>% 
    filter(Intra_TF > 0.5)
  
  if (nrow(top_df) == 0) {
    top_df <- 
      top_tf_frac(count_mat_mm$Reduced_resize, tf, meta_mm) %>% 
      rownames_to_column(var = "Region") %>% 
      filter(Intra_TF > 0.2)
  }
  
  
  top_df$Group <- get_ccre_group(top_df, ccre_mm)
  
  p <- 
    ggplot(top_df, aes(x = Intra_TF, y = Inter_TF, fill = Group)) +
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
  
})
names(plist_mm) <- tfs_mm


pdf(paste0(plot_dir, date, "_mouse_freq_bound_by_ccre.pdf"))
invisible(lapply(plist_mm, print))
graphics.off()
