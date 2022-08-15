## This script plots an intronic region of the SHB gene using the trackplot
## package to demonstrate a region with specific ASCL1 binding
## -----------------------------------------------------------------------------

library(tidyverse)
library(GenomicRanges)
library(trackplot)
source("~/regnetR/R/utils/range_table_functions.R")
source("~/regnetR/R/utils/plot_functions.R")


# installing bwtools: https://gist.github.com/PoisonAlien/e19b482ac6146bfb03142a0de1c4fbc8
# saved this in my own /home/user/bin, then added this dir to R path:
bwtool_path <- "/home/amorin/bin/bwtool/"
if (!str_detect(Sys.getenv("PATH"), bwtool_path)) {
  Sys.setenv(PATH = paste(Sys.getenv("PATH"), bwtool_path, sep = ":"))
}

date <- "Apr2022"  # latest update when metadata is saved out

plot_dir <- "~/Plots/Chipseq/Trackplots/"
run_dirs <- read.delim(paste0("~/Data/Metadata/Chipseq/batch1_run_dirs_", date, ".tsv"), stringsAsFactors = FALSE)
meta <- read.delim(paste0("~/Data/Metadata/Chipseq/batch1_chip_meta_final_", date, ".tsv"), stringsAsFactors = FALSE)
chip_dir <- "/cosmos/data/pipeline-output/chipseq-encode-pipeline/chip/"

# For inspection of coordinates and how often a TR dataset overlappeds
region_list <- readRDS(paste0("~/scratch/R_objects/", date, "_count_mat_list.RDS"))
regions <- region_list$Human$Reduced_resize
regions <- regions[order(regions[, "ASCL1"], decreasing = TRUE), ]
locus <- "chr9:38022988-38023319"


tf <- "Ascl1"

ncores <- 8
pad_window <- 500  # how many bps should be added to start and end of region


get_bw_files <- function(run_id,
                         dir_location = run_dirs,
                         bw_track = "pval.signal.bigwig",
                         pooled_reps = "include") {
  
  # Returns the full path of bigwig files produced by the ENCODE pipeline. For 
  # the given Run ID. Assumes typical dir structure of a run after organize meta
  # has been run. 
  # run_id: the metadata Experiment_ID associated with an ENCODE run
  # dir_location: table that associates the run with its dir location
  # bw_track: either "pval.signal.bigwig" or "fc.signal.bigwig" (pval or foldchange)
  # pooled_reps: if multiple reps exist, return the pooled rep? one of exclude", "include", or "only"
  
  path <- paste0(filter(dir_location, Experiment_ID == run_id)$Dir, "/signal")
  reps <- list.files(path, full.names = TRUE)
  
  if (length(reps) == 0) stop("No files for the provided ID")
  stopifnot(pooled_reps %in% c("exclude", "include", "only"))
  
  if (pooled_reps == "exclude") {
    reps <- reps[!str_detect(reps, "pooled-rep")]
  } else if (pooled_reps == "only") {
    if (length(reps) > 1) {
      reps <- reps[str_detect(reps, "pooled-rep")]
    }
  }
  
  rep_bw <- unlist(lapply(reps, function(rep) {
    list.files(rep, pattern = bw_track, full.names = TRUE)
  }))
  
  return(rep_bw)
}


pad_locus <- function(locus, window_size) {
  # Increase the size of the locus by provided window size
  split_locus <- str_split(locus, ":|-", simplify = TRUE)
  start <- as.integer(split_locus[, 2]) - window_size
  end <- as.integer(split_locus[, 3]) + window_size
  return(paste0(split_locus[1], ":", start, "-", end))
}



# Get ASCL1 IDs, sample 2 IDs from rest of TRs, and tack on the 2 TCF4 experiments
# that also had a peak called

set.seed(13)

id_ascl1 <- meta %>%
  filter(Species == "Human" & Symbol == "ASCL1") %>%
  pull(Experiment_ID) %>%
  unique()

id_sample <- meta %>% 
  distinct(Experiment_ID, .keep_all = TRUE) %>% 
  filter(Species == "Human" & Symbol != "ASCL1") %>% 
  group_by(Symbol) %>% 
  slice_sample(n = 2) %>% 
  select(Symbol, Experiment_ID) %>% 
  ungroup()

tcf4 <- meta %>% 
  distinct(Experiment_ID, .keep_all = TRUE) %>% 
  filter(Experiment_ID %in% c("GSE119475_TCF4_Human_SU-DHL-2", "GSE119476_TCF4_Human_TMD8")) %>% 
  select(Symbol, Experiment_ID)


id_sample <- rbind(id_sample, tcf4)


# trackplot colours by TR

col_ascl1 <- tf_pal[tf]
col_sample <- tf_pal[str_to_title(id_sample$Symbol)]


# get full paths to big wig files

bw_ascl1 <- unlist(lapply(id_ascl1, function(x) get_bw_files(x, pooled_reps = "only")))
bw_sample <- unlist(lapply(id_sample$Experiment_ID, function(x) get_bw_files(x, pooled_reps = "only")))


# Extract track data

track_ascl1 <- track_extract(
  bigWigs = bw_ascl1,
  loci = pad_locus(locus, pad_window),
  binsize = 10,
  custom_names = id_ascl1,
  nthreads = ncores
)


track_sample <- track_extract(
  bigWigs = bw_sample,
  loci = pad_locus(locus, pad_window),
  binsize = 10,
  custom_names = id_sample$Experiment_ID,
  nthreads = ncores
)


# For plotting sample tracks, use the median max value of the ASCL1 experiments
# to emphasize the comparative lack of signal (save for the 2 TCF4 experiments)

track_max <- median(unlist(lapply(track_ascl1, function(x) max(data.frame(x)$max))))


# Plot tracks

# ASCL1
png(height = 20, width = 12, units = "in", res = 300,
    filename = paste0(plot_dir, "Human_ASCL1_", locus, ".png"))

track_plot(
  summary_list = track_ascl1,
  col = col_ascl1,
  build = "hg38",
  collapse_txs = FALSE,
  groupAutoScale = FALSE,  # each track has its own ylim scale
  track_names = id_ascl1,
  track_names_pos = 1
)

graphics.off()


# Sample
png(height = 20, width = 12, units = "in", res = 300,
    filename = paste0(plot_dir, "Human_sample-TR_", locus, ".png"))

track_plot(
  summary_list = track_sample,
  col = col_sample,
  build = "hg38",
  collapse_txs = FALSE,
  y_max = track_max,
  track_names = id_sample$Experiment_ID,
  track_names_pos = 1
)
graphics.off()
