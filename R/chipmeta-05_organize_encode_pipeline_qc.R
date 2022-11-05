## This script organizes and saves two tables from the ENCODE pipeline qc2tsv 
## output.The first table cleans/organizes selected stats which are output by 
## sample. The second averages these stats for samples within a run, and adds 
## additional information about that run/experiment (count of peaks, number of 
## samples) and merges it with the curated metadata.
## -----------------------------------------------------------------------------

library(tidyverse)
source("R/setup-01_config.R")

# The raw output of the qc2tsv tool
qc_input <- paste0(pipeout_dir, "qc_reports/", date, "_chip_qc_report.tsv")

# Cleaned QC info for all samples
qc_sample_output <- paste0("~/Data/Metadata/Chipseq/batch1_clean_sample_qc_", date, ".tsv")

# Cleaned QC info at the experiment level (samples averaged)
qc_run_output <- paste0("~/Data/Metadata/Chipseq/batch1_chip_meta_completed_withqc_", date, ".tsv")

# The final metadata that includes completed experiments with their QC metrics
meta_final_output <- paste0("~/Data/Metadata/Chipseq/batch1_chip_meta_final_", date, ".tsv")

# The metadata tables before adding the QC metrics
meta_all <- read.delim(paste0("~/Data/Metadata/Chipseq/batch1_chip_meta_all_", date, ".tsv"), stringsAsFactors = FALSE)
meta_distinct <- read.delim(paste0("~/Data/Metadata/Chipseq/batch1_chip_meta_completed_runs_", date, ".tsv"), stringsAsFactors = FALSE)


# Reminder - current imp doesn't auto run qc2tsv after chipmeta-04
if (!file.exists(qc_input)) {
  stop("File does not exist - ensure 'qc2tsv' has been run")
} else {
  qc <- read.delim(qc_input, skip = 4, stringsAsFactors = FALSE)
}


# Only want to keep most immediately relevant QC fields. For some fields like
# title (corresponding to Experiment_ID) only the first element is named and the 
# rest of that set are blank. Need to fill out. Also must be aware of 
# redundancies in the the naming of QC output columns. R adds a digit suffix
# ------------------------------------------------------------------------------


keep_cols <- c(
  "total_reads",  # samstat
  "mapped_reads", # samstat
  "pct_mapped_reads", # samstat
  "mapped_reads.1", # no_dup samstat
  "pct_duplicate_reads",
  "reproducibility",  # overlap reproducibility
  "reproducibility.1",  # idr reproducibility
  "num_peaks",  # MACS2 pval 0.01 peaks for reproducible analysis
  "N_opt", # number of optimal overlap reproducible peaks
  "N_opt.1", # number of optimal IDR reproducible peaks
  "NRF",
  "NSC",
  "RSC"
)

stopifnot(all(keep_cols %in% colnames(qc)))


# Helper to fill out a character vector with blanks - assumes that the blanks
# should be filled with the last non blank character

fill_chr <- function(chr_vec) {

  if (chr_vec[1] == "") {
    stop("Can't fill if the first element is blank")
  }
  
  for (i in 1:length(chr_vec)) {
    chr <- chr_vec[i]
    if (chr != "" ) {
      save_chr <- chr
    } else {
      chr_vec[i] <- save_chr
    }
  }
  return(chr_vec)
}


# Fill out N IDR/Overlap peaks and reproducible peak status for all non-input
# samples of a run (since only the first sample of a run has this info, and it
# may get lost when splitting samples)

fill_repr <- function(qc_df) {
  
  qc_list <- split(qc_df, qc_df$Experiment_ID)
  
  qc_list <- lapply(qc_list, function(x) {
    x$N_opt[which(!x$Is_input)] <- x$N_opt[1]
    x$N_opt.1[which(!x$Is_input)] <- x$N_opt.1[1]
    x$reproducibility[which(!x$Is_input)] <- x$reproducibility[1]
    x$reproducibility.1[which(!x$Is_input)] <- x$reproducibility.1[1]
    return(x)
  })
  
  qc_df <- do.call(rbind, qc_list)
  rownames(qc_df) <- NULL
  return(qc_df)
}



# Given a qc df and the desired columns to keep, return a df of equal nrow
# to qc_df and only with the requested columns, as well as the SRX ID,
# Run title, and whether or not the sample is an input control

organize_qc <- function(qc_df, keep_cols) {
  
  ids <- unlist(str_split(qc_df$description, ", "))
  ids <- ids[ids != ""]
  
  is_input <- str_detect(qc_df$replicate, "ctl[:digit:]")
  
  runs <- fill_chr(qc_df$title)
  
  qc_df <- data.frame(ID = ids,
                      Is_input  = is_input,
                      Experiment_ID = runs,
                      qc_df[, keep_cols],
                      stringsAsFactors = FALSE)
  
  qc_df <- fill_repr(qc_df)
  
  qc_df <- dplyr::rename(
    qc_df,
    Overlap_reproducibility = reproducibility,
    IDR_reproducibility = reproducibility.1,
    mapped_reads_nodup = mapped_reads.1,
    N_relaxed_peaks = num_peaks,
    N_overlap_peaks = N_opt,
    N_IDR_peaks = N_opt.1
  )
  
  return(qc_df)
  
}


qc_sub <- organize_qc(qc, keep_cols)



# Now average selected QC info over samples within a run/experiment to summarize.
# NOTE: Histone mode didn't output this information: can match the completed
# Mecp2 histone runs with the equivalent TF run IF the TF run also completed.
# Because not all Mecp2 samples succeeded the TF mode, they will be absent
# ------------------------------------------------------------------------------


qc_samples <- filter(qc_sub, !Is_input)
qc_input <- filter(qc_sub, Is_input) # this will have duplicates b/c some inputs shared across runs


# Fill in information from mecp2 TF runs, if available 

mecp2_na <- qc_samples$ID[is.na(qc_samples$mapped_reads_nodup)]

qc_na <- qc_sub %>% 
  filter(ID %in% mecp2_na & !is.na(mapped_reads_nodup)) %>% 
  dplyr::select(c(ID, total_reads, mapped_reads, pct_mapped_reads, 
                  mapped_reads_nodup, pct_duplicate_reads, NRF))


qc_samples <-
  full_join(qc_samples, qc_na, by = "ID", suffix = c("", ".y")) %>%
  dplyr::mutate(
    total_reads = ifelse(is.na(total_reads), total_reads.y, total_reads),
    mapped_reads = ifelse(is.na(mapped_reads), mapped_reads.y, mapped_reads),
    pct_mapped_reads = ifelse(is.na(pct_mapped_reads), pct_mapped_reads.y, pct_mapped_reads),
    mapped_reads_nodup = ifelse(is.na(mapped_reads_nodup), mapped_reads_nodup.y, mapped_reads_nodup),
    pct_duplicate_reads = ifelse(is.na(pct_duplicate_reads), pct_duplicate_reads.y, pct_duplicate_reads),
    NRF = ifelse(is.na(NRF), NRF.y, NRF)
  ) %>% 
  dplyr::select(!ends_with('.y'))


# Sample counts within an experiment 

samps_per_run <- meta_all %>% 
  group_by(Experiment_ID) %>% 
  dplyr::summarise(Count_samples = n())

input_per_run <- meta_all %>%
  group_by(Experiment_ID) %>%
  dplyr::summarise(Count_input = sum(!is.na(unique(unlist(str_split(Input_ID, ", "))))))


# For peak counts, consider average relaxed peak counts (applied to each sample),
# overlap and IDR counts (which are per run/experiment), and a count that uses 
# IDR for TFs and overlap for Mecp2 (as this is what is used in analysis)


peaks <- qc_samples %>% 
  group_by(Experiment_ID) %>% 
  dplyr::mutate(Avg_relaxed_peaks = round(mean(N_relaxed_peaks, na.rm = TRUE), 1)) %>% 
  distinct(Experiment_ID, .keep_all = TRUE) %>%
  dplyr::mutate(N_peaks = ifelse(str_detect(str_to_lower(Experiment_ID), "mecp2"),
                               N_overlap_peaks, 
                               N_IDR_peaks)) %>% 
  dplyr::select(Experiment_ID, Avg_relaxed_peaks, N_overlap_peaks, N_IDR_peaks, N_peaks)


# Then average rest of statistics of interest

exp_reads <- qc_samples %>% 
  group_by(Experiment_ID) %>% 
  dplyr::mutate(Avg_exp_mapped_reads = mean(mapped_reads),
                Avg_exp_mapped_reads_nodup = mean(mapped_reads_nodup)) %>% 
  distinct(Experiment_ID, .keep_all = TRUE) %>%
  dplyr::select(Experiment_ID, Avg_exp_mapped_reads, Avg_exp_mapped_reads_nodup)


input_reads <- qc_input %>% 
  group_by(Experiment_ID) %>% 
  dplyr::mutate(Avg_input_mapped_reads = mean(mapped_reads),
                Avg_input_mapped_reads_nodup = mean(mapped_reads_nodup)) %>% 
  distinct(Experiment_ID, .keep_all = TRUE) %>%
  dplyr::select(Experiment_ID, Avg_input_mapped_reads, Avg_input_mapped_reads_nodup)


exp_nrf <-  qc_samples %>% 
  group_by(Experiment_ID) %>% 
  dplyr::mutate(Avg_exp_NRF = mean(NRF)) %>% 
  distinct(Experiment_ID, .keep_all = TRUE) %>%
  dplyr::select(Experiment_ID, Avg_exp_NRF)


input_nrf <-  qc_input %>% 
  group_by(Experiment_ID) %>% 
  dplyr::mutate(Avg_input_NRF = mean(NRF)) %>% 
  distinct(Experiment_ID, .keep_all = TRUE) %>%
  dplyr::select(Experiment_ID, Avg_input_NRF)


nsc <-  qc_samples %>% 
  group_by(Experiment_ID) %>% 
  dplyr::mutate(Avg_NSC = mean(NSC)) %>% 
  distinct(Experiment_ID, .keep_all = TRUE) %>%
  dplyr::select(Experiment_ID, Avg_NSC)


rsc <-  qc_samples %>% 
  group_by(Experiment_ID) %>% 
  dplyr::mutate(Avg_RSC = mean(RSC)) %>% 
  distinct(Experiment_ID, .keep_all = TRUE) %>%
  dplyr::select(Experiment_ID, Avg_RSC)


# join back together for run-level QC
qc_run <- Reduce(function(...) merge(..., by = "Experiment_ID", all.x = TRUE), 
                 list(meta_distinct,
                      peaks,
                      exp_reads, 
                      input_reads,
                      samps_per_run,
                      input_per_run,
                      exp_nrf,
                      input_nrf,
                      nsc,
                      rsc))


# Save out copy of only experiments considered for analysis (passed IDR/overlap
# and have more than min peaks)

# NOTE: Uncertain why ENCODE pipeline will spit out IDR failure in the terminal
# for most cases, but in others the run will complete and IDR failure is only
# reflected in the QC report.


fail_idr <- qc_samples %>% 
  filter(!(str_detect(str_to_lower(Experiment_ID), "mecp2")) & 
           IDR_reproducibility == "fail") %>% 
  pull(Experiment_ID) %>% 
  unique()


fail_overlap <- qc_samples %>%
  filter(str_detect(str_to_lower(Experiment_ID), "mecp2") &
           Overlap_reproducibility == "fail") %>%
  pull(Experiment_ID) %>% 
  unique()


fail_npeak <- unique(filter(qc_run, N_peaks < min_peaks)$Experiment_ID)


meta_final <- qc_run %>% 
  filter(!(Experiment_ID %in% c(fail_idr, fail_overlap, fail_npeak))) %>% 
  arrange(Symbol, Species)


# Save out

write.table(qc_sub,
            sep = "\t",
            quote = FALSE,
            file = qc_sample_output)

write.table(qc_run,
            sep = "\t",
            quote = FALSE,
            file = qc_run_output)

write.table(meta_final,
            sep = "\t",
            quote = FALSE,
            file = meta_final_output)
