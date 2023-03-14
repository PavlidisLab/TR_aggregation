## Explore overlap of ChIP-seq GR objects within and between TR groups. Create
## an "all" GR object that collapses peaks from every experiment, +/- resizing
## peaks relative to summit and +/- reduction (collapsing overlapping peaks).
## Also exports a matrix of these "all" regions and how many data sets from each
## TR overlaps with these regions
## -----------------------------------------------------------------------------

library(tidyverse)
library(GenomicRanges)
library(parallel)
library(ggbio)
source("R/setup-01_config.R")
source("R/utils/range_table_functions.R")

window_size <- 150  # padding to add to either direction of peak summit

# GRanges objects
gr_hg <- readRDS(grlist_hg_path)
gr_mm <- readRDS(grlist_mm_path)

# batch 1 ChIP-seq meta
meta <- read.delim(chip_meta_path, stringsAsFactors = FALSE)

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


# Look at summary of peak sizes across all experiments
# ------------------------------------------------------------------------------


get_width_df <- function(gr_list, meta) {
  
  width_df <- data.frame(do.call(
    rbind, lapply(gr_list, function(x) summary(width(x)))
    )) %>% 
    rownames_to_column(var = "Experiment_ID") %>%
    left_join(meta[, c("Symbol", "Experiment_ID")], by = "Experiment_ID")
  
  rownames(width_df) <- meta$Experiment_ID
  return(width_df)
}



# Human: GSE69394_NEUROD1_Human_H82 highest max peak size (123kb - 8:127638060-127761390)
# which is stretch 5' of MYC that contains non protein coding CASC11
# GSE58428_RUNX1_Human_VCaP highest mean peak size (1.67kb) - has high RSC

width_hg <- get_width_df(gr_hg, meta_hg)

max_hg <- width_hg[which.max(width_hg$Max.),]
max_peak_hg <- gr_hg[rownames(max_hg)][width(gr_hg[rownames(max_hg)]) == max_hg$Max.]

max_mean_hg <- width_hg[which.max(width_hg$Mean),]
meta[meta$Experiment_ID %in% c(rownames(max_mean_hg), names(max_peak_hg)), ]



# Mouse: Largest max peak size GSE90893_Runx1_Mouse_MEF-Control
# (73kb - 1:71582211-71655087) 
# Largest mean peak size belong to GSE126375_Runx1_Mouse_Pax5-KO (1.7kb)


width_mm <- get_width_df(gr_mm, meta_mm)

max_mm <- width_mm[which.max(width_mm$Max.),]
max_peak_mm <- gr_mm[rownames(max_mm)][width(gr_mm[rownames(max_mm)]) == max_mm$Max.]

max_mean_mm <- width_mm[which.max(width_mm$Mean),]
meta[meta$Experiment_ID %in% c(rownames(max_mean_mm), names(max_peak_mm)), ]


# Look at how many peaks overlap within data sets before fixing size.
# ------------------------------------------------------------------------------


# Check if reduced set (collapses overlapping peaks) is same length as original

any_overlaps <- function(gr_list) {
  
  vapply(gr_list, function(gr) {
    return(length(gr) == length(reduce(gr)))
  }, FUN.VALUE = logical(1))

}


# Find no intra-peak overlaps (which is expected based on pipeline output).

stopifnot(all(any_overlaps(gr_hg)))
stopifnot(all(any_overlaps(gr_mm)))


# For each data set, find the smallest distance between peaks


nearest_dist <- function(gr_list, cores) {
  
  unlist(mclapply(gr_list, function(x) {
    min(distanceToNearest(x)@elementMetadata$distance)
  }, mc.cores = cores))
  
}


# Closest peaks across all datasets is 33bp for both human and mouse, all in
# Runx1 experiments

nearest_hg <- nearest_dist(gr_hg, cores = cores)
nearest_hg[nearest_hg == min(nearest_hg)]

nearest_mm <- nearest_dist(gr_mm, cores = cores)
nearest_mm[nearest_mm == min(nearest_mm)]



# Because peaks have variable sizes, fix to summit position and add padding to
# make all same size for overlap (1bp summit too stringent for overlap)
# ------------------------------------------------------------------------------


# re-size
gr_rs_hg <- GRangesList(lapply(gr_hg, summit_window, window_size))
gr_rs_mm <- GRangesList(lapply(gr_mm, summit_window, window_size))



# Note that re-sizing 100bp+ results in overlapping peaks that were not
# overlapping in their original forms. In the "reduced" set these re-sized
# peaks are then merged, despite the fact that they originally represented two
# summits from a single data set. Following was explored using
# "GSE61197_ASCL1_Human_H1775-ASCL1-high"
# ------------------------------------------------------------------------------



# Make a standard Grange plot and add a line corresponding to the summit

plot_gr_summit <- function(gr) {
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
  split(test_gr_after, rep(1:length(test_gr_after), each = 1)), test_gr_after)

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
p_max_before <- plot_gr_summit(test_gr_before[(max_ix - 1):max_ix])
p_max_after <- plot_gr_summit(test_gr_after[(max_ix - 1):max_ix])


# Describe the full and reduced (non-overlapping) set of peaks for all 
# experiments with and without the re-size.
# Note that when peaks have been re-sized, the reduced set can create large 
# peaks that were not-overlapping in their original forms (summits near the 
# edge close to another peak, and the window size padding causes overlap)
# ------------------------------------------------------------------------------


# Human: All and re-size ~3.5m peaks across 129 experiments which drops to ~1.22m
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

count_hg <- unlist(lapply(all_hg, length))
change_hg <- (count_hg["Reduced"] - count_hg["All"]) / count_hg["All"]
change_rs_hg <- (count_hg["Reduced_resize"] - count_hg["Resize"]) / count_hg["Resize"]

maxpeak_hg <- all_hg$Reduced[which.max(width(all_hg$Reduced))]
maxpeak_rs_hg <- all_hg$Reduced_resize[which.max(width(all_hg$Reduced_resize))]

table(meta_hg$Symbol[unique(findOverlaps(maxpeak_hg, gr_hg)@to)])
table(meta_hg$Symbol[unique(findOverlaps(maxpeak_rs_hg, gr_rs_hg)@to)])


# Mouse: All and re-size ~1.95m peaks across 126 experiments which drops to ~686k
# in the re-sized and reduced set (-65%) and ~588k in the original reduced set (-70%)

# The longest reduced peak in the original reduced set is 73kb 
# (1: 71582211-71655460) which overlaps peaks in 33 data sets.
# The longest in the re-sized reduced set is 3.3kb (7: 127024978-127028247)
# in 56 data sets


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


# all_gr is a single gr of all experiments, exp_gr are the same experiments
# grouped in a list. 

all_overlap <- function(all_gr, exp_gr, meta) {
  
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
# Most frequently overlapped regions 
# Original/All - 11: 65496125-65509326 (13 kbp) overlapped in 94 data sets
# Resize - 19:44753276-44753763 (488 bp) overlapped in 90 data sets
# Reduced - 16:85604925-85645912 (41 kb) overlapped in 96 data sets
# Reduced resize - 19:44753276-44753763 (488 bp) overlapped in 90 data sets


all_ol_hg <- list(
  All = all_overlap(all_hg$All, gr_hg, meta_hg),
  Resize = all_overlap(all_hg$Resize, gr_rs_hg, meta_hg),
  Reduced = all_overlap(all_hg$Reduced, gr_hg, meta_hg),
  Reduced_resize = all_overlap(all_hg$Reduced_resize, gr_rs_hg, meta_hg)
)

lapply(all_ol_hg, `[[`, "Most_overlapped")
lapply(all_ol_hg, `[[`, "Most_overlapped_TF")
lapply(all_ol_hg, `[[`, "Count_exp_most_overlapped")
lapply(all_ol_hg, `[[`, "Most_overlapped_width")


# Mouse: 
# Most frequently overlapped regions 
# Original/All - 3: 88507136-88515082 (8 kbp) overlapped in 82 data sets
# Resize - 3: 88509931-88510302 (372 bp) overlapped in 79 data sets
# Reduced - 3: 88507004-88515082 (8 kb) overlapped in 82 data sets
# Reduced resize - 3: 88509708-88510703 (996 bp) overlapped in 79 data sets



all_ol_mm <- list(
  All = all_overlap(all_mm$All, gr_mm, meta_mm),
  Resize = all_overlap(all_mm$Resize, gr_rs_mm, meta_mm),
  Reduced = all_overlap(all_mm$Reduced, gr_mm, meta_mm),
  Reduced_resize = all_overlap(all_mm$Reduced_resize, gr_rs_mm, meta_mm)
)

lapply(all_ol_mm, `[[`, "Most_overlapped")
lapply(all_ol_mm, `[[`, "Most_overlapped_TF")
lapply(all_ol_mm, `[[`, "Count_exp_most_overlapped")
lapply(all_ol_mm, `[[`, "Most_overlapped_width")



# Create TF-specific all and reduced GR objects
# ------------------------------------------------------------------------------


gr_by_tf <- function(gr_list, meta) {
  
  tfs <- unique(meta$Symbol)
  
  tf_gr <- lapply(tfs, function(tf) {
    runs <- meta[meta$Symbol == tf, "Experiment_ID"]
    unlist(gr_list[runs])
  })
  names(tf_gr) <- tfs
  
  return(tf_gr)
}


tf_gr_hg <- list(
  All = gr_by_tf(gr_hg, meta_hg),
  Resize = gr_by_tf(gr_rs_hg, meta_hg))
tf_gr_hg$Reduced <- lapply(tf_gr_hg$All, reduce)
tf_gr_hg$Reduced_resize <- lapply(tf_gr_hg$Resize, reduce)


tf_gr_mm <- list(
  All = gr_by_tf(gr_mm, meta_mm),
  Resize = gr_by_tf(gr_rs_mm, meta_mm))
tf_gr_mm$Reduced <- lapply(tf_gr_mm$All, reduce)
tf_gr_mm$Reduced_resize <- lapply(tf_gr_mm$Resize, reduce)



# Decrease in size when reduced. a small % change suggests the TF experiments
# were different, since the reduced set retains unique peaks. a large % change
# suggests the experiments were similar - overlapping peaks get reduced
# ------------------------------------------------------------------------------


tf_count_df <- function(tf_gr_list) {
  data.frame(
    do.call(
      cbind, lapply(tf_gr_list, function(x) unlist(lapply(x, length)))
    )
  )
}


count_change <- function(count_df) {
  data.frame(
    Symbol = rownames(count_df),
    Reduced = (count_df[, "Reduced"] - count_df[, "All"]) / count_df[, "All"],
    Reduced_resize = (count_df[, "Reduced_resize"] - count_df[, "Resize"]) / count_df[, "Resize"]
  )
}



# Human: RUNX1 most reduced size, HES1 least

tf_count_hg <- tf_count_df(tf_gr_hg)
tf_change_hg <- count_change(tf_count_hg)

# Mouse: Runx1 most reduced size, Mecp2 least (Tcf4 only has 1 data set)

tf_count_mm <- tf_count_df(tf_gr_mm)
tf_change_mm <- count_change(tf_count_mm)



# Get the count of TF data sets that each range overlaps. Arrange by counts 
# for each TF while minimizing sum of counts for the other TFs. Suggestive of
# TF-specific binding sites
# ------------------------------------------------------------------------------


# Return a matrix where rows are every region of all_gr (all ChIP-seq experiments
# as one GRange object), and columns are the unique TRs in meta. Each element
# of this matrix is the number of experiments for the given TR that had a peak
# overlapping the given range.


get_count_mat <- function(all_gr, gr_list, meta, cores) {
  
  # find which data sets of gr_list each range of all_gr overlaps
  gr_ol <- findOverlaps(all_gr, gr_list)
  ol_by_ix <- as(gr_ol, "List")
  
  # get the corresponding TF symbol of the overlapping hit
  ol_by_tf <- extractList(meta$Symbol, ol_by_ix)
  
  # counts to table - as factor to maintain tally of TFs with 0 counts
  ol_by_count <- mclapply(ol_by_tf, function(x) {
    table(factor(x, levels = tfs))
  }, mc.cores = cores)
  
  # return the matrix of counts with rows named after the range
  count_mat <- as.matrix(do.call(rbind, ol_by_count))
  rownames(count_mat) <- paste0(seqnames(all_gr), ":", start(all_gr), "-", end(all_gr))
  
  return(count_mat)
}


# NOTE: slow. limiting step appears to be tallying the symbols for each
# range in get_count_mat()


if (!file.exists(ol_count_path)) {
  
  count_mat_hg <- list(
    All = get_count_mat(all_hg$All, gr_hg, meta_hg, cores),
    Resize = get_count_mat(all_hg$Resize, gr_rs_hg, meta_hg, cores),
    Reduced = get_count_mat(all_hg$Reduced, gr_hg, meta_hg, cores),
    Reduced_resize = get_count_mat(all_hg$Reduced_resize, gr_rs_hg, meta_hg, cores)
  )
  
  count_mat_mm <- list(
    All = get_count_mat(all_mm$All, gr_mm, meta_mm, cores),
    Resize = get_count_mat(all_mm$Resize, gr_rs_mm, meta_mm, cores),
    Reduced = get_count_mat(all_mm$Reduced, gr_mm, meta_mm, cores),
    Reduced_resize = get_count_mat(all_mm$Reduced_resize, gr_rs_mm, meta_mm, cores)
  )
  
  saveRDS(list(Human = count_mat_hg, Mouse = count_mat_mm), ol_count_path) 
  
  } else {
  
  dat <- readRDS(ol_count_path)
  count_mat_hg <- dat$Human
  count_mat_mm <- dat$Mouse

}
