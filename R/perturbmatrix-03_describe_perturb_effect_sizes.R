## This script summarizes and plots various effect sizes of perturbation
## experiments: the perturbed TF's fold change, experiments with unexpected
## changes in expression for the perturbed TF, and the count of DEGs at FDR05
## -----------------------------------------------------------------------------

library(tidyverse)
library(cowplot)
library(googlesheets4)
source("R/setup-01_config.R")
source("R/utils/gemma_functions.R")
source("R/utils/plot_functions.R")


# FC cutoff for demonstration, not used in final analysis.
# fdr (false discovery rate control for DE) declared in config
fc <- 1

plot_dir <- paste0(pplot_dir, "Effect_size/")

# Load meta and list of DE result sets +/- filtering
meta <- read.delim(file =  paste0(meta_dir, "batch1_tfperturb_meta_final_", date, ".tsv"), stringsAsFactors = FALSE)
results <- readRDS(paste0(expr_dir, "TF_perturb_batch1_rslist_", date, ".RDS"))
results_unfilt <-  readRDS(paste0(expr_dir, "TF_perturb_batch1_unfiltered_rslist_", date, ".RDS"))

stopifnot(identical(names(results), meta$Experiment_ID))


# Extract from the result sets the rows corresponding to the perturbed TF of 
# interest. Also consider the unfiltered results to see if there are multiple 
# mapped elements to the TF in divergent directions
# ------------------------------------------------------------------------------


# init blank row to sub in for TFs missing data
blank_result <- results[[1]][1, ]
blank_result[1,] <- NA

# Isolate stats for the perturbed TF from the filtered result sets

pert_list <- lapply(1:length(results), function(x) {
  df <- filter_symbol(results[[x]], meta$Symbol[x])
  if (nrow(df) == 0) df <- blank_result
  df$Element_Name <- as.character(df$Element_Name) # otherwise mix of int/chr that messes list bind
  df$NCBI_ID <- as.character(df$NCBI_ID)
  return(df)
})
names(pert_list) <- names(results)


# Ditto but for the unfiltered results

pert_list_unfilt <- lapply(1:length(results_unfilt), function(x) {
  filter_symbol(results_unfilt[[x]], meta$Symbol[x])
})
names(pert_list) <- names(results)


# Look at unfiltered results to find if multiple probes match to TF
which_multi <- unlist(lapply(pert_list_unfilt, function(x) (nrow(x) > 1)))
which_multi_accession <- meta[which_multi, c("GSE", "Symbol", "Experiment_ID", "Platform")]

# which experiments do not have any data for the TF of interest
which_none <- unlist(lapply(pert_list, function(x) is.na(x$FoldChange)))
which_none_accession <- meta[which_none, c("GSE", "Symbol", "Experiment_ID", "Platform")]


# Inspect the experiments with no associated elements using unfiltered results
# and fuzzy string matching to see if elements are mapped to multiple genes
# ------------------------------------------------------------------------------


# GSE64264 has no probes mapping to PAX6

filter_symbol(results_unfilt$GSE64264_PAX6_Human_Knockdown, "PAX6", exact = FALSE)

# GSE79598 has an ambiguously mapped probe to RUNX1 (this was found by manually
# searching the associated platform  on Gemma). Coerce to RUNX1

pert_list$GSE79598_RUNX1_Human_Knockdown <- 
  filter(results_unfilt$GSE79598_RUNX1_Human_Knockdown, Element_Name == "TC21000423.hg.1") %>% 
  mutate(Symbol = "RUNX1", Adj_pval = as.numeric(NA), PercRankFC = as.numeric(NA))
  
pert_list$`GSE79598_RUNX1_Human_Knockdown-1` <- 
  filter(results_unfilt$`GSE79598_RUNX1_Human_Knockdown-1`, Element_Name == "TC21000423.hg.1") %>% 
  mutate(Symbol = "RUNX1", Adj_pval = as.numeric(NA), PercRankFC = as.numeric(NA))

# GSE54295 has an ambiguously mapped probe, but ENSEMBL lumps to RUNX1. Going
# to keep this multi mapped probe and coerce name to just RUNX1

pert_list$GSE54295_RUNX1_Human_Mutant <- 
  filter(results_unfilt$GSE54295_RUNX1_Human_Mutant, Element_Name == "16925365") %>% 
  mutate(Symbol = "RUNX1", Adj_pval = as.numeric(NA), PercRankFC = as.numeric(NA))

# GSE151385 no match to RUNX1 - knockdown resulted in gene being filtered?

filter_symbol(results_unfilt$GSE151385_RUNX1_Human_Knockdown, "RUNX1", exact = FALSE)

# Consistency for downstream binding
pert_list <- lapply(pert_list, as.data.frame)  


# Constructing a summary data frame of the effect sizes of each perturbation.
# ------------------------------------------------------------------------------


pert_df <- data.frame(
  do.call(rbind, pert_list)) %>% 
  rownames_to_column(var = "Experiment_ID") %>% 
  mutate(PercRankFC = as.double(PercRankFC),
         Species = meta$Species,
         Symbol = meta$Symbol,
         Perturbation = meta$Perturbation)

pert_df$Count_DEG <- vapply(results, function(x) {
  sum(x$Adj_pval < fdr)
}, FUN.VALUE = double(1))

pert_df$Count_DEG_FC <- vapply(results, function(x) {
  sum(x$Adj_pval < fdr & abs(x$FoldChange) > fc)
}, FUN.VALUE = double(1))


# summary of DEGs
summary(pert_df$Count_DEG)
summary(pert_df$Count_DEG_FC)

# Count of 0 DEGs, and check by species and perturbation type
sum(pert_df$Count_DEG == 0)/nrow(pert_df)
table(pert_df[pert_df$Count_DEG == 0, "Perturbation"])
table(pert_df[pert_df$Count_DEG == 0, "Species"])

# Test perturbation type for count of DEGs +/- requiring at least 1 DEG
pert_df_min <- filter(pert_df, Count_DEG > 0)
kruskal.test(pert_df$Count_DEG ~ pert_df$Perturbation)
kruskal.test(pert_df_min$Count_DEG ~ pert_df_min$Perturbation)

# Note that knockdowns have the most DEGs
pert_df %>% group_by(Perturbation) %>% summarize(Median = median(Count_DEG))
pert_df_min %>% group_by(Perturbation) %>% summarize(Median = median(Count_DEG))

pert_df %>% 
  mutate(Is_KD = (Perturbation == "Knockdown")) %>% 
  group_by(Is_KD) %>% 
  summarize(Median = median(Count_DEG))


write.table(
  x = pert_df,  
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  file = paste0(expr_dir, date, "_perturbation_effect_sizes.tsv")
)


# Which TFs have unexpected perturbations (upregulated after KD/KO or downreg
# after OE)? If multiple probes match are they all opposite of expectations?
# Split by single/multi probes assigned to the perturbed TF 
# ------------------------------------------------------------------------------


# Returns a 1x5 df that summarizes the probe values for perturbed TFs that are
# not in the expected direction of change. This is fed into an lapply call and 
# assumes that meta_df is a 1 row slice that corresponds to the supplied 
# result set table

summarize_multi_probes <- function(meta_df, resultset) {
  
  perturb <- meta_df$Perturbation
  mean_fc <- mean(resultset$FoldChange)
  above_zero <- sum(resultset$FoldChange > 0)
  below_zero <- sum(resultset$FoldChange < 0)
  
  data.frame(
    Perturbation = perturb,
    Experiment_ID = meta_df$Experiment_ID,
    Mean_FC = mean_fc,
    Above_zero = above_zero,
    Below_zero = below_zero
  )
}


which_unexpected <- which(
  (
    pert_df$Perturbation %in% c("Knockdown", "Knockout") & 
      pert_df$FoldChange > 0
  ) |
    (
      pert_df$Perturbation == "Overexpression" &
        pert_df$FoldChange < 0
    )
)

# Isolate unexpected experiments in summary df and unfiltered results list

unexpected_df <- pert_df[which_unexpected, c("Experiment_ID", "FoldChange", "PercRankFC")]
pert_list_unexpected <- pert_list_unfilt[which_unexpected]
names(pert_list_unexpected) <- names(results)[which_unexpected]

# Isolate those with multiple probes and summarize their direction/mean FC

which_multi_unexpected <- intersect(which_multi, which_unexpected)

multi_unexpected <- lapply(which_multi_unexpected, function(x) {
  summarize_multi_probes(meta[x, ], pert_list_unfilt[[x]])
})

multi_unexpected_df <- do.call(rbind, multi_unexpected)

# look at the percentile rank fold change of the perturbed TF - if unexpected
# in direction, is it still highly ranked?

summary(unexpected_df$PercRankFC)  # median 0.84


# save out meta of unexpected for curation/inspection
# ------------------------------------------------------------------------------


# write_sheet(data = filter(unexpected_df, !Experiment_ID %in% multi_unexpected_df$Experiment_ID),
#             ss = gsheets_perturb,
#             sheet = paste0("Unexpected_single_probe_", date)
# )
# 
# write_sheet(data = multi_unexpected_df,
#             ss = gsheets_perturb,
#             sheet = paste0("Unexpected_multi_probe_", date)
# )


# Get a df of fold changes for each experiment for plotting
# -----------------------------------------------------------------------------


get_fc_df <- function(results_list, meta) {
  
  fc_l <- lapply(1:length(results_list), function(x) {
    
    data.frame(
      FoldChange = results_list[[x]][, "FoldChange"],
      Experiment_ID = meta$Experiment_ID[x],
      Symbol = str_to_upper(meta$Symbol[x])
    )
  })
  
  return(do.call(rbind, fc_l))
}


fc_df <- get_fc_df(results, meta)


# Correlations of perturb effect size ~ count DEGs
# -----------------------------------------------------------------------------

pert_df_nona <- pert_df[!is.na(pert_df$FoldChange), ]

# percentile rank fold change
cor.test(pert_df_nona$PercRankFC, pert_df_nona$Count_DEG, method = "spearman")

# fold change +/- abs
cor.test(pert_df_nona$FoldChange, pert_df_nona$Count_DEG, method = "spearman")
cor.test(abs(pert_df_nona$FoldChange), pert_df_nona$Count_DEG, method = "spearman")

# require count DEG > 0

# abs fold change - all
cor.test(
  abs(pert_df_nona[pert_df_nona$Count_DEG > 0, "FoldChange"]),
  pert_df_nona[pert_df_nona$Count_DEG > 0, "Count_DEG"],
  method = "spearman"
)

# abs fold change - exclude mutants
cor.test(
  abs(pert_df_nona[pert_df_nona$Count_DEG > 0 & pert_df_nona$Perturbation != "Mutant", "FoldChange"]),
  pert_df_nona[pert_df_nona$Count_DEG > 0 & pert_df_nona$Perturbation != "Mutant", "Count_DEG"],
  method = "spearman"
)

# require at least log2FC >= 1 for the perturbed TF and exclude mutants

# percentile rank fold change
cor.test(
  pert_df_nona[abs(pert_df_nona$FoldChange) > 1 & pert_df_nona$Perturbation != "Mutant", "PercRankFC"],
  pert_df_nona[abs(pert_df_nona$FoldChange) > 1 & pert_df_nona$Perturbation != "Mutant", "Count_DEG"],
  method = "spearman"
)

# fold change
cor.test(
  abs(pert_df_nona[abs(pert_df_nona$FoldChange) > 1 & pert_df_nona$Perturbation != "Mutant", "FoldChange"]),
  pert_df_nona[abs(pert_df_nona$FoldChange) > 1 & pert_df_nona$Perturbation != "Mutant", "Count_DEG"],
  method = "spearman"
)


# Plot
# ------------------------------------------------------------------------------


# remove NAs and collapse mouse and human symbols

plot_df <-  mutate(pert_df, Symbol = str_to_upper(Symbol))
plot_df_nona <- filter(plot_df, !is.na(FoldChange))

species_shape <- c("Human" = 21, "Mouse" = 24)

# Scatter plot of the perturbed TF fold change, grouped by perturbation and 
# colored by TF.

p1 <- 
  ggplot(plot_df_nona, aes(x = Experiment_ID, y = FoldChange, fill = Symbol)) +
  geom_point(size = 4, color = "black", shape = 21) +
  scale_x_discrete(expand = c(0.1, 0.1)) +  # points were getting cut off
  geom_hline(yintercept = 0, linetype = "dashed") +
  facet_grid(. ~ Perturbation) +
  ylab("Log2 fold change of perturbed TR") +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 25),
        axis.title = element_text(size = 25),
        strip.text.x = element_text(size = 20),
        strip.background = element_rect(colour = "black", fill = "azure3"),
        panel.grid = element_blank(),  # to work with expand in scale x
        legend.text = element_text(size = 25),
        legend.title = element_blank()) +
  scale_fill_manual(values = tf_pal_hg)


ggsave(p1, height = 8, width = 13, dpi = 300, device = "png",
       filename = paste0(plot_dir, "Perturbed_TF_FC_", date, ".png"))


# Scatter plot of the perturbed TF PercRankFC, grouped by perturbation and 
# colored by TF.

p2 <- 
  ggplot(plot_df_nona, aes(x = Experiment_ID, y = PercRankFC, fill = Symbol)) +
  geom_point(size = 4, color = "black", shape = 21) +
  scale_x_discrete(expand = c(0.1, 0.1)) +  # points were getting cut off
  facet_grid(. ~ Perturbation) +
  ylab("Percentile rank of absolute fold change") +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 25),
        axis.title = element_text(size = 25),
        strip.text.x = element_text(size = 20),
        strip.background = element_rect(colour = "black", fill = "azure3"),
        panel.grid = element_blank(),  # to work with expand in scale x
        legend.text = element_text(size = 25),
        legend.title = element_blank()) +
  scale_fill_manual(values = tf_pal_hg)


ggsave(p2, height = 8, width = 12, dpi = 300, device = "png",
       filename = paste0(plot_dir, "Perturbed_TF_PercRankFC_", date, ".png"))


# plot the count of DEGs for each experiment

# col by TF
p3a <- 
  ggplot(plot_df, aes(x = Experiment_ID, y = log10(Count_DEG + 1), fill = Symbol)) +
  geom_point(size = 4.6, color = "black", shape = 21) +
  facet_grid(. ~ Symbol) +
  scale_x_discrete(expand = c(0.2, 0.2)) +  # points were getting cut off
  ylab(paste0("log10 Count of DEGs at FDR", fdr)) +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),        
        axis.text.y = element_text(size = 25),
        axis.ticks.x = element_blank(),
        axis.title = element_text(size = 25),
        strip.text = element_text(size = 18),
        strip.background = element_rect(colour = "black", fill = "cornsilk2"),
        panel.grid = element_blank(), # to work with expand in scale x
        legend.position = "none") +
  scale_fill_manual(values = tf_pal_hg)

# col by perturb
p3b <- 
  ggplot(plot_df, aes(x = Experiment_ID, y = log10(Count_DEG + 1), fill = Perturbation)) +
  geom_point(size = 5.5, color = "black", shape = 21) +
  facet_grid(. ~ Symbol) +
  scale_x_discrete(expand = c(0.2, 0.2)) +  # points were getting cut off
  ylab(paste0("log10 Count of DEGs at FDR", fdr)) +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),        
        axis.text.y = element_text(size = 25),
        axis.ticks.x = element_blank(),
        axis.title = element_text(size = 25),
        strip.text = element_text(size = 18),
        strip.background = element_rect(colour = "black", fill = "cornsilk2"),
        panel.grid = element_blank(), # to work with expand in scale x
        legend.text = element_text(size = 25),
        legend.title = element_blank()) +
  scale_fill_manual(values = pert_anno$Perturbation)

# col by perturb and shape by species (used in paper)
p3c <- 
  ggplot(plot_df, aes(x = Experiment_ID, y = log10(Count_DEG + 1), fill = Perturbation, shape = Species)) +
  geom_point(size = 5.5, color = "black") +
  facet_grid(. ~ Symbol) +
  scale_x_discrete(expand = c(0.2, 0.2)) +  # points were getting cut off
  ylab(paste0("log10 Count of DEGs at FDR", fdr)) +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),        
        axis.text.y = element_text(size = 25),
        axis.ticks.x = element_blank(),
        axis.title = element_text(size = 25),
        strip.text = element_text(size = 18),
        strip.background = element_rect(colour = "black", fill = "cornsilk2"),
        panel.grid = element_blank(),  # to work with expand in scale x
        legend.text = element_text(size = 25),
        legend.title = element_blank()) +
  scale_fill_manual(values = pert_anno$Perturbation) +
  scale_shape_manual(values = species_shape) +
  guides(fill = guide_legend(override.aes = list(shape = 21)),
         shape = guide_legend(override.aes = list(fill = "white")))

p3c_noleg <- p3c + theme(legend.position = "none")

ggsave(p3a, height = 8, width = 12, dpi = 300, device = "png",
       filename = paste0(plot_dir, "DEG_counts_byTF_", date, ".png"))

ggsave(p3b, height = 8, width = 12, dpi = 300, device = "png",
       filename = paste0(plot_dir, "DEG_counts_byPert_", date, ".png"))

ggsave(p3c, height = 6, width = 14, dpi = 300, device = "png",
       filename = paste0(plot_dir, "DEG_counts_byPertandSpecies_", date, ".png"))

ggsave(p3c_noleg, height = 6, width = 12, dpi = 300, device = "png",
       filename = paste0(plot_dir, "DEG_counts_byPertandSpecies_nolegend_", date, ".png"))



# Scatter plot of the relationships between perturb effect size and count DEGs


# PRFC (used in paper)
p4 <- 
  ggplot(plot_df_nona, aes(x = PercRankFC, y = log10(Count_DEG + 1), fill = Symbol)) +
  geom_jitter(size = 5, color = "black", shape = 21) +
  ylab(paste0("log10 Count of DEGs at FDR", fdr)) +
  xlab("Percentile rank fold change of perturbed TR") +
  theme_classic() +
  theme(axis.text = element_text(size = 25),
        axis.title = element_text(size = 25),
        legend.text = element_text(size = 25),
        legend.title = element_blank()) +
  scale_fill_manual(values = tf_pal_hg)

ggsave(p4, height = 9, width = 12, dpi = 300, device = "png",
       filename = paste0(plot_dir, "DEG_counts_vs_PercRankFC_", date, ".png"))

# fold change +/- abs
p5a <- 
  ggplot(plot_df_nona, aes(x = FoldChange, y = log10(Count_DEG + 1), fill = Symbol)) +
  geom_jitter(size = 5, color = "black", shape = 21) +
  ylab(paste0("log10 Count of DEGs at FDR", fdr)) +
  xlab("Log2 fold change of perturbed TR") +
  theme_classic() +
  theme(axis.text = element_text(size = 25),
        axis.title = element_text(size = 25),
        legend.text = element_text(size = 25),
        legend.title = element_blank()) +
  scale_fill_manual(values = tf_pal_hg)

p5b <- 
  ggplot(plot_df_nona, aes(x = abs(FoldChange), y = log10(Count_DEG + 1), fill = Symbol)) +
  geom_jitter(size = 5, color = "black", shape = 21) +
  ylab(paste0("log10 Count of DEGs at FDR", fdr)) +
  xlab("Absolute log2 fold change of perturbed TR") +
  theme_classic() +
  theme(axis.text = element_text(size = 25),
        axis.title = element_text(size = 25),
        legend.text = element_text(size = 25),
        legend.title = element_blank()) +
  scale_fill_manual(values = tf_pal_hg)

ggsave(p5a, height = 9, width = 12, dpi = 300, device = "png",
       filename = paste0(plot_dir, "DEG_counts_vs_FC_", date, ".png"))

ggsave(p5b, height = 9, width = 12, dpi = 300, device = "png",
       filename = paste0(plot_dir, "DEG_counts_vs_absFC_", date, ".png"))


# plot the dist of fold changes for each experiment


# order by median FC within TR group

fc_df <- fc_df %>% 
  group_by(Symbol, Experiment_ID) %>% 
  mutate(median = median(FoldChange)) %>% 
  arrange(Symbol, median)

fc_df$Experiment_ID <- factor(fc_df$Experiment_ID, levels = unique(fc_df$Experiment_ID))

p6a <- 
  ggplot(fc_df, aes(x = Experiment_ID, y = FoldChange, fill = Symbol)) +
  geom_boxplot(outlier.shape = NA, width = 1) +
  geom_hline(yintercept = c(-1, 1), color = "red", size = 1.5) +
  theme_classic() +
  ylim(-3, 3) +
  ylab("Log2 fold change") +
  xlab("TR perturbation experiments") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title = element_text(size = 25),
        axis.text.y = element_text(size = 25),
        legend.text = element_text(size = 25),
        legend.title = element_blank()) +
  scale_fill_manual(values = tf_pal_hg)

p6b <- p6a + theme(legend.position = "none")


ggsave(p6a, height = 6, width = 15, dpi = 300, device = "png",
       filename = paste0(plot_dir, "DEG_counts_boxplot_", date, ".png"))

ggsave(p6b, height = 6, width = 12, dpi = 300, device = "png",
       filename = paste0(plot_dir, "DEG_counts_boxplot_nolegend_", date, ".png"))


# boxplot of count of DEGs by perturbation type

p7a <- pert_df %>% 
  ggplot(., aes(x = Perturbation, y = log10(Count_DEG + 1))) +
  geom_boxplot(width = 0.4, fill = "turquoise4") +
  ylab(paste0("log10 Count of DEGs at FDR", fdr)) +
  # ggtitle("All experiments") +
  theme_classic() +
    theme(axis.title.x = element_blank(),
          axis.title = element_text(size = 25),
          axis.text = element_text(size = 20),
          plot.title = element_text(hjust = 0.5, size = 20))


p7b <- pert_df %>% 
  filter(Count_DEG > 0) %>% 
  ggplot(., aes(x = Perturbation, y = log10(Count_DEG + 1))) +
  geom_boxplot(width = 0.4, fill = "turquoise4") +
  ylab(paste0("log10 Count of DEGs at FDR", fdr)) +
  # ggtitle("Excluding experiments with 0 DEGs") +
  theme_classic() +
  theme(axis.title.x = element_blank(),
        axis.title = element_text(size = 25),
        axis.text = element_text(size = 20),
        plot.title = element_text(hjust = 0.5, size = 20))


ggsave(p7a, height = 6, width = 9, dpi = 300, device = "png",
       filename = paste0(plot_dir, "DEG_counts_boxplot_byperturb_allexp_", date, ".png"))

ggsave(p7b, height = 6, width = 9, dpi = 300, device = "png",
       filename = paste0(plot_dir, "DEG_counts_boxplot_byperturb_exclnodeg_", date, ".png"))
