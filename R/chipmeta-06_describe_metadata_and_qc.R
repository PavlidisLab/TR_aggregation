## This script summarizes and plots information from the curated ChIP-seq metadata
## and resuting QC report. Summary plots are saved out.
## -----------------------------------------------------------------------------

library(tidyverse)
library(cowplot)
source("R/setup-01_config.R")
source("R/utils/plot_functions.R")

plot_dir <- paste0(cplot_dir, "Describe_meta/")

# meta final is what is used for analysis - filtered for min peaks and passing IDR/overlap
# meta runs is all completed runs (including those under min peak) - used for plotting
meta_final <- read.delim(paste0(meta_dir, "Chipseq/batch1_chip_meta_final_", date, ".tsv"), stringsAsFactors = FALSE)
meta_runs <- read.delim(paste0(meta_dir, "Chipseq/batch1_chip_meta_completed_withqc_", date, ".tsv"), stringsAsFactors = FALSE)
meta_all <- read.delim(paste0(meta_dir, "Chipseq/batch1_chip_meta_all_", date, ".tsv"), stringsAsFactors = FALSE)


# Describe metadata
# ------------------------------------------------------------------------------

# Number of samples relative to unique completed runs

n_runs <- nrow(meta_final)

summary(meta_final$Count_samples)

frac_runs_1sample <- sum(meta_final$Count_samples == 1)/n_runs

# Count of samples by TF

n_samps_species <- meta_final %>% 
  mutate(Symbol = str_to_title(Symbol)) %>% 
  group_by(Symbol, Species) %>% 
  dplyr::summarise(n = sum(Count_samples))

n_samps_all <- meta_final %>% 
  mutate(Symbol = str_to_title(Symbol)) %>% 
  group_by(Symbol) %>% 
  dplyr::summarise(n = sum(Count_samples))


# Sample/runs by species

table(meta_final$Species)

meta_final %>% 
  group_by(Species) %>% 
  summarise(n = sum(Count_samples))

# How many different runs did each GSE accession have

runs_per_geo <- meta_final %>% 
  group_by(GSE) %>% 
  summarise(Count_runs = n_distinct(Experiment_ID))

summary(runs_per_geo)


# Failures for TF/non-Mecp2

meta_all %>% 
  filter(Complete != "1" & str_to_title(Symbol) != "Mecp2") %>% 
  distinct(Experiment_ID, .keep_all = TRUE) %>% 
  dplyr::count(Complete)


# Failures for Mecp2 histone mode

meta_all %>% 
  filter(Complete != "1" & str_to_title(Symbol) == "Mecp2") %>% 
  filter(str_detect(Experiment_ID, "HISTONE")) %>% 
  distinct(Experiment_ID, .keep_all = TRUE) %>% 
  dplyr::count(Complete)

# Runs with input

n_inputs <- sum(!is.na(
  unique(unlist(str_split(meta_final$Input_ID, ", ")))))

summary(meta_final$Count_input)

frac_runs_noinput <- sum(meta_final$Count_input == 0)/n_runs

# Experiments with condition (1: TF is directly perturbed, 2: Some treatment)

table(meta_final$Condition)

# Sum of condition across samples broken down by TF

meta_final %>% 
  dplyr::mutate(Symbol = str_to_title(Symbol)) %>% 
  group_by(Symbol, Condition) %>% 
  summarise(n())

# Ratio of condition by TF for experiments/runs

meta_final %>%
  dplyr::mutate(Symbol = str_to_title(Symbol)) %>% 
  group_by(Symbol, Condition) %>% 
  summarise(N = n()) %>% 
  dplyr::mutate(Ratio = round(N / sum(N), 2))

# Cell types

sort(table(meta_final$Cell_Type), decreasing = TRUE)

# Summary of read counts for ChIP'ed samples and inputs. 

summary(meta_final[, c("Avg_exp_mapped_reads_nodup", "Avg_input_mapped_reads_nodup")] / 1e6)

# Summary of QC metrisc

summary(meta_final[, c("Avg_exp_NRF", "Avg_input_NRF", "Avg_NSC", "Avg_RSC")])


# Correlation and testing of peak counts with technical/experiment design factors
# NOTE - N_peaks is IDR for TF, overlap for mecp2, which is what is used in paper
# ------------------------------------------------------------------------------


round(cor(dplyr::select_if(meta_final, is.numeric), use = "pairwise.complete.obs"), 3)

# Peaks vs de-duplicated mapped read counts (averaged over replicates)

cor.test(meta_final$Avg_exp_mapped_reads_nodup, meta_final$N_peaks, 
         use = "pairwise.complete.obs",
         method = "spearman")

# difference in peak counts for +/- input control

# Wilcox reported in paper
wilcox.test(meta_final$N_peaks ~ (meta_final$Count_input > 0))
t.test(log10(meta_final$N_peaks+1) ~ (meta_final$Count_input > 0))

# difference in peak counts for +/- reps

# Wilcox reported in paper
wilcox.test(meta_final$N_peaks ~ (meta_final$Count_samples > 1))
t.test(log10(meta_final$N_peaks+1) ~ (meta_final$Count_samples > 1))


# Plot
# ------------------------------------------------------------------------------


# Stacked bar chart of count of samples

p1 <- n_samps_species %>% 
  ggplot(., aes(x = reorder(Symbol, n), y = n, fill = Species)) +
  geom_bar(stat = "identity", colour = "black", width = 0.8) +
  ylab("Count of samples") +
  scale_y_continuous(limits = c(0, 200), breaks = seq(0, 175, 25)) +
  theme_classic() +
  scale_fill_manual(values = c("royalblue", "goldenrod")) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 30),
        axis.text.y = element_text(size = 30),
        axis.text.x = element_text(size = 30, angle = 60, vjust = 1, hjust = 1),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15),
        legend.position = "bottom")


ggsave(p1,
       dpi = 300,
       device = "png",
       height = 10,
       width = 12,
       file = paste0(plot_dir, "batch1_chip_sample_counts_all_", date, ".png"))


# Count of experiments as stacked bar

p2 <- meta_final %>% 
  dplyr::mutate(Symbol = str_to_title(Symbol)) %>% 
  distinct(Experiment_ID, .keep_all = TRUE) %>% 
  count(Symbol, Species) %>% 
  ggplot(., aes(x = reorder(Symbol, n), y = n, fill = Species)) +
  geom_bar(stat = "identity", colour = "black", width = 0.8) +
  ylab("Count of experiments") +
  scale_y_continuous(limits = c(0, 150), breaks = seq(0, 150, 15)) +
  theme_classic() +
  scale_fill_manual(values = c("royalblue", "goldenrod")) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 30),
        axis.text.y = element_text(size = 30),
        axis.text.x = element_text(size = 30, angle = 60, vjust = 1, hjust = 1),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15),
        legend.position = "bottom")



p2 <- meta_final %>%
  dplyr::mutate(Symbol = str_to_title(Symbol)) %>%
  distinct(Experiment_ID, .keep_all = TRUE) %>%
  count(Symbol, Species) %>%
  ggplot(., aes(y = reorder(Symbol, n), x = n, fill = Species)) +
  geom_bar(stat = "identity", colour = "black", width = 0.8) +
  xlab("Count of experiments") +
  scale_x_continuous(limits = c(0, 150), breaks = seq(0, 150, 15)) +
  theme_classic() +
  scale_fill_manual(values = c("royalblue", "goldenrod")) +
  theme(axis.title.y = element_blank(),
        axis.title.x = element_text(size = 30),
        axis.text.y = element_text(size = 30),
        axis.text.x = element_text(size = 30, angle = 60, vjust = 1, hjust = 1),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15),
        legend.position = "bottom")



ggsave(p2,
       dpi = 300,
       device = "png",
       height = 10,
       width = 12,
       file = paste0(plot_dir, "batch1_chip_experiment_counts_all_", date, ".png"))


# Counts of peaks as jitter plots with median bars, grouped by TR or presence of input.
# Note use of meta_runs, which includes experiments that were ultimately filtered


plot_jitter_median <- function(df, xvar, yvar, ylab, hline = NULL) {
  
  ggplot(df) +
    geom_jitter(aes(x = !!sym(xvar), y = !!sym(yvar), fill = Species),
                color = "black", shape = 21, width = 0.1, size = 4) +
    stat_summary(aes(x = !!sym(xvar), y = !!sym(yvar)),
                 fun = median, fun.min = median, fun.max = median, 
                 geom = "crossbar", width = 0.5, inherit.aes = FALSE) +
    geom_hline(yintercept = hline, linetype = "dashed", color = "black") +
    ylab(ylab) +
    scale_fill_manual(values = c("royalblue", "goldenrod")) +
    theme_classic() +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_text(size = 25),
          axis.text = element_text(size = 20),
          legend.text = element_text(size = 20),
          legend.title = element_text(size = 20))
}


# overlap counts ~  symbol 

p3a <- meta_runs %>% 
  dplyr::mutate(Symbol = str_to_title(Symbol),
                Peaks = log10(N_overlap_peaks)) %>% 
  plot_jitter_median(ylab = "Log10 count of overlap peaks", yvar = "Peaks", xvar = "Symbol")


ggsave(p3a,
       dpi = 300,
       device = "png",
       height = 8,
       width = 12,
       file = paste0(plot_dir, "batch1_chip_count_overlap_peaks_by_symbol_", date, ".png"))


# IDR counts

p3b <- meta_runs %>% 
  dplyr::mutate(Symbol = str_to_title(Symbol),
         Peaks = log10(N_IDR_peaks)) %>% 
  plot_jitter_median(ylab = "Log10 count of IDR peaks", yvar = "Peaks", xvar = "Symbol")

ggsave(p3b,
       dpi = 300,
       device = "png",
       height = 8,
       width = 12,
       file = paste0(plot_dir, "batch1_chip_count_IDR_peaks_by_symbol_", date, ".png"))

# overlap for mecp2, IDR for TF (this is what is used in the paper)

p3c <- meta_runs %>% 
  dplyr::mutate(Symbol = str_to_title(Symbol),
         Peaks = log10(N_peaks)) %>% 
  plot_jitter_median(ylab = "Log10 count of reproducible peaks",
                     yvar = "Peaks", xvar = "Symbol",
                     hline = log10(min_peaks))

ggsave(p3c,
       dpi = 300,
       device = "png",
       height = 8,
       width = 12,
       file = paste0(plot_dir, "batch1_chip_count_tf-idr_mecp2_overlap_peaks_by_symbol_", date, ".png"))


# counts ~ binary has input as jitter+median - here overlap for mecp2, idr for tf


p4 <- meta_final %>% 
  dplyr::mutate(Symbol = str_to_title(Symbol),
         Has_input = Count_input > 0) %>% 
  ggplot(.) +
  geom_jitter(aes(x = Has_input, y = log10(N_peaks), fill = Symbol),
              color = "black", shape = 21, size = 2.5, width = 0.1) +
  stat_summary(aes(x = Has_input, y = log10(N_peaks)), 
               fun = median, fun.min = median, fun.max = median, 
               geom = "crossbar", width = 0.3, inherit.aes = FALSE) +
  ylab("Log10 count of reproducible peaks") + 
  xlab("Has input control") +
  scale_fill_manual(values = tf_pal) +
  theme_classic() +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 20))

ggsave(p4,
       dpi = 300,
       device = "png",
       height = 6,
       width = 4,
       file = paste0(plot_dir, "batch1_chip_count_tf-idr_mecp2_overlap_peaks_by_input_", date, ".png"))


# Correlation of read counts and peaks
# Note that this will be missing mecp2 samples that did not complete TF pipeline
# (histone pipeline did not output the qc.json)


# IDR peaks ~ average non-dup mapped reads as scatter

p5a <- meta_final %>% 
  filter(!is.na(Avg_exp_mapped_reads_nodup) & !is.na(N_IDR_peaks)) %>% 
  dplyr::mutate(Peaks = log10(N_IDR_peaks), Reads = log10(Avg_exp_mapped_reads_nodup)) %>%
  ggplot(aes(y = Peaks, x = Reads)) +
  geom_point() +
  geom_smooth(method = lm) +
  ylab("Log10 count of reproducible peaks") +
  xlab("Log10 count of average mapped reads") +
  theme_classic() +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 20))

ggsave(p5a,
       dpi = 300,
       device = "png",
       height = 6,
       width = 6,
       file = paste0(plot_dir, "batch1_chip_IDR_peaks_vs_avg_exp_mreads", date, ".png"))


# overlap peaks vs mapped reads

p5b <- meta_final %>% 
  filter(!is.na(Avg_exp_mapped_reads_nodup) & !is.na(N_overlap_peaks)) %>% 
  dplyr::mutate(Peaks = log10(N_overlap_peaks), Reads = log10(Avg_exp_mapped_reads_nodup)) %>% 
  ggplot(aes(y = Peaks, x = Reads)) +
  geom_point() +
  geom_smooth(method = lm) +
  ylab("Log10 count of reproducible peaks") +
  xlab("Log10 count of average mapped reads") +
  theme_classic() +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 20))

ggsave(p5b,
       dpi = 300,
       device = "png",
       height = 6,
       width = 6,
       file = paste0(plot_dir, "batch1_chip_overlap_overlap_peaks_vs_avg_exp_mreads", date, ".png"))


# overlap for mecp2, IDR for TF (used in paper)

p5c <- meta_final %>%
  filter(!is.na(Avg_exp_mapped_reads_nodup) & !is.na(N_peaks)) %>%
  dplyr::mutate(Symbol = str_to_title(Symbol)) %>% 
  ggplot(aes(y = log10(N_peaks), x = log10(Avg_exp_mapped_reads_nodup))) +
  geom_point() +
  geom_smooth(method = lm) +
  ylab("Log10 count of reproducible peaks") +
  xlab("Log10 count of average mapped reads") +
  theme_classic() +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 20))

ggsave(p5c,
       dpi = 300,
       device = "png",
       height = 6,
       width = 6,
       file = paste0(plot_dir, "batch1_chip_count_tf-idr_mecp2_overlap_peaks_vs_avg_exp_mreads", date, ".png"))
