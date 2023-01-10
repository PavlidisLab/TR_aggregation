## Loads matrices of perturb effect sizes, creates summary data frames of the
## count of times genes were DE across studies, plots/summarizes this info, and
## exports these summary dataframes to be used in the final rankings.
## -----------------------------------------------------------------------------

library(tidyverse)
library(parallel)
library(cowplot)
library(ggExtra)
library(ggrepel)
source("R/setup-01_config.R")
source("R/utils/perturbmatrix_functions.R")
source("R/utils/plot_functions.R")

plot_dir <- paste0(pplot_dir, "Describe_FDR_counts/")

# Output list of DE counts, grouped by TR or for all experiment combined
tf_outfile <- paste0(expr_dir, "TF_perturb_DE_counts_list_by_TF_FDR01_", date, ".RDS")
all_outfile <- paste0(expr_dir, "TF_perturb_DE_counts_list_all_FDR01_", date, ".RDS")

# Load meta and lists of perturb effect size matrices
meta <- read.delim(file =  paste0(meta_dir, "batch1_tfperturb_meta_final_", date, ".tsv"), stringsAsFactors = FALSE)
mlist_hg <- readRDS(paste0(pmat_dir, "human_list_perturb_matrix_", date, ".RDS"))
mlist_mm <- readRDS(paste0(pmat_dir, "mouse_list_perturb_matrix_", date, ".RDS"))
mlist_ortho <- readRDS(paste0(pmat_dir, "ortho_list_perturb_matrix_", date, ".RDS"))

# ortho protein coding genes
pc_ortho <- read.delim(paste0(meta_dir, "hg_mm_1to1_ortho_genes_DIOPT-v8.tsv"), stringsAsFactors = FALSE)

# DE prior tables
pdeg_hg <- read.delim(paste0(meta_dir, "DE_prior_hg.tsv"), stringsAsFactors = FALSE)
pdeg_mm <- read.delim(paste0(meta_dir, "DE_prior_mm.tsv"), stringsAsFactors = FALSE)

# Split meta by species
meta_mm <- meta[meta$Species == "Mouse", ]
meta_hg <- meta[meta$Species == "Human", ]


# List of data frames summarizing gene stats across perturbation experiments.
# Generate for all experiments and split by TF experiments. 
#
# Count_DE: Tally of gene being diff expr at FDR cutoff across experiments
# Proportion_DE_measured: Prop. of DE across experiments with non-NAs for gene
# Proportion_DE_all: Prop. of DE across all experiments
# Count_NA: Tally of NA measurements for the gene across experiments
# GoF_down: Tally when a gene had FC < 0 in overexpression experiments
# GoF_up: Tally when a gene had FC > 0 in overexpression experiments
# LoF_down: Tally when a gene had FC < 0 in KO and KD experiments
# LoF_up: Tally when a gene had FC > 0 in KO and KD experiments
# Avg_abs_FC: A gene's average absolute fold change across experiments
# FC_purity: Measure of consistency of FC direction within GoF/LoF experiments
# FC_signed_purity: FC_purity but make negative if GoF and LoF agree in direction
# DE_prior_rank: Ranks genes by how likely they are to be DE in diverse studies
# ------------------------------------------------------------------------------


all_de <- list(
  
  Human = process_all(fdr_mat = mlist_hg$FDR_mat, 
                      fdr = fdr, 
                      fc_mat = mlist_hg$FC_mat, 
                      meta = meta_hg, 
                      de_prior = pdeg_hg),
  
  Mouse = process_all(fdr_mat = mlist_mm$FDR_mat, 
                      fdr = fdr, 
                      fc_mat = mlist_mm$FC_mat, 
                      meta = meta_mm, 
                      de_prior = pdeg_mm),
  
  Ortho = process_all(fdr_mat = mlist_ortho$FDR_mat, 
                      fdr = fdr, 
                      fc_mat = mlist_ortho$FC_mat, 
                      meta = meta)
)


tf_de <- list(
  
  Human = process_tf(fdr_mat = mlist_hg$FDR_mat, 
                     fdr = fdr, 
                     fc_mat = mlist_hg$FC_mat, 
                     meta = meta_hg, 
                     de_prior = pdeg_hg),
  
  Mouse = process_tf(fdr_mat = mlist_mm$FDR_mat, 
                     fdr = fdr, 
                     fc_mat = mlist_mm$FC_mat, 
                     meta = meta_mm, 
                     de_prior = pdeg_mm),
  
  Ortho = process_tf(fdr_mat = mlist_ortho$FDR_mat, 
                     fdr = fdr, 
                     fc_mat = mlist_ortho$FC_mat, 
                     meta = meta, 
                     ortho = TRUE)
)

# Post-hoc addition of species-specific DE counts

tf_de$Ortho <- merge_ortho_counts(count_ortho = tf_de$Ortho, 
                                  count_hg = tf_de$Human, 
                                  count_mm = tf_de$Mouse, 
                                  pc_ortho = pc_ortho)



# Describe DE/NA counts across all experiments
# NOTE: While all NAs have been filtered, some genes still have mostly (eg, all
# but 1) NAs - be aware that this can influence trends when looking across exps.
# ------------------------------------------------------------------------------


# summary of times a gene was DE across all experiments for hg/mm/ortho
all_summ <- lapply(all_de, function(x) summary(x$Count_DE))

# count of genes that were DE at least once
all_min1 <- lapply(all_de, function(x) sum(x$Count_DE > 0))

# count of genes that were never measured as DE
all_none <- lapply(all_de, function(x) sum(x$Count_DE == 0))

# count of times a gene was not measured across experiments
all_na_summ <- lapply(all_de, function(x) summary(x$Count_NA))


# Inspect relationship of DE prior and DE counts
# ------------------------------------------------------------------------------


# Correlation of DE prior and DE count

cor.test(all_de$Human$DE_prior_rank,
         all_de$Human$Count_DE,
         use = "pairwise.complete.obs")


cor.test(all_de$Mouse$DE_prior_rank,
         all_de$Mouse$Count_DE,
         use = "pairwise.complete.obs")


# Genes with max DE counts across all experiments

de_max_hg <- slice_max(all_de$Human, Count_DE)$Symbol

de_max_mm <- slice_max(all_de$Mouse, Count_DE)$Symbol

max_df_hg <- 
  do.call(rbind, lapply(tf_de$Human, function(x) filter(x, Symbol %in% de_max_hg)))

max_df_mm <- 
  do.call(rbind, lapply(tf_de$Mouse, function(x) filter(x, Symbol %in% de_max_mm)))


# Inspect and provide examples of high DE count and low DE prior: Genes that are
# frequently DE in this collection, but not in the global DE prior

prior_lo_hg <- all_de$Human %>%
  filter(Count_DE > 20) %>%
  arrange(DE_prior_rank, desc(Count_DE))

prior_lo_mm <- all_de$Mouse %>%
  filter(Count_DE > 30) %>%
  arrange(DE_prior_rank, desc(Count_DE))


low_hg <- c("NPEPPS", "ATE1", "AAGAB", "INPPL1")
low_mm <- c("Irak1", "Pofut2")

low_hg_df <- 
  do.call(rbind, lapply(tf_de$Human, function(x) filter(x, Symbol %in% low_hg))) %>% 
  arrange(Symbol)

low_mm_df <- 
  do.call(rbind, lapply(tf_de$Mouse, function(x) filter(x, Symbol %in% low_mm))) %>% 
  arrange(Symbol)


# Inspect relationship of count DE and average absolute FC
# In general see weak cor, but in some instances (eg mouse Mecp2) see weak 
# negative cor. Appears to be driven by genes that are mostly NA (so will have
# low DE counts regardless) and high FC in the limited studies they are msrd.
# ------------------------------------------------------------------------------


cor(all_de$Human[, c("Count_DE", "Avg_abs_FC", "Count_NA")])
summary(lm(Count_DE ~ Avg_abs_FC + Count_NA, data = all_de$Human))
lapply(tf_de$Human, function(x) cor.test(x$Count_DE, x$Avg_abs_FC))

cor(all_de$Mouse[, c("Count_DE", "Avg_abs_FC", "Count_NA")])
summary(lm(Count_DE ~ Avg_abs_FC + Count_NA, data = all_de$Mouse))
lapply(tf_de$Mouse, function(x) cor.test(x$Count_DE, x$Avg_abs_FC))


# Genes with high average absolute FC and high DE count
# ------------------------------------------------------------------------------


# Filter df for var in top quantile (qtl)

qtl_filter <- function(df, var, qtl = 0.9) {
  filter(df, !!sym(var) >= quantile(!!sym(var), qtl, na.rm = TRUE))
}


# Heuristic - remove genes in the top quantile (qtl) of count NAs

na_filter <- function(df, qtl = 0.25) {
  filter(df, Count_NA <= quantile(df$Count_NA, qtl))
}


fc_hi_hg <- lapply(tf_de$Human, function(x) {
  qtl_filter(x, var = "Count_DE") %>% qtl_filter(var = "Avg_abs_FC")
})


fc_hi_mm <- lapply(tf_de$Mouse, function(x) {
  qtl_filter(x, var = "Count_DE") %>% qtl_filter(var = "Avg_abs_FC")
})


# High FC and high DE count tend to have extremely high DE prior. Example
# of low is MEF2C-CLBA1 (DEPR = 0.09, DE 2/4).

lapply(fc_hi_hg, function(x) summary(x$DE_prior_rank))
lapply(fc_hi_mm, function(x) summary(x$DE_prior_rank))
filter(fc_hi_hg$MEF2C, DE_prior_rank == min(DE_prior_rank, na.rm = TRUE))


# Examples of genes with high abs FC and low/no DE count across exps - largely
# genes with many NAs

tf_de$Mouse$Runx1 %>% 
  arrange(desc(Avg_abs_FC), Count_DE) %>% 
  head(10)


# Inspect FC purity
# ------------------------------------------------------------------------------


# Genes with high purity (GoF/LoF gives same FC direction) and elevated DE 
# counts. Remove genes with frequent NAs.

pur_hi_hg <- lapply(tf_de$Human, function(x) {
  na_filter(x) %>% 
  qtl_filter(var = "FC_purity") %>%
  arrange(desc(Avg_abs_FC))
})

pur_hi_mm <- lapply(tf_de$Mouse, function(x) {
  na_filter(x) %>% 
  qtl_filter(var = "FC_purity") %>%
  arrange(desc(Avg_abs_FC))
})

# ASCL1-CABP7: DE in 6/8 experiments with a Signed purity == 1

filter(pur_hi_hg$ASCL1, Symbol == "CABP7")


# Genes with high DE count and low purity

pur_lo_hg <- lapply(tf_de$Human, function(x) {
  filter(x, FC_purity < 0.6) %>% qtl_filter(var = "Count_DE")
})

pur_lo_mm <- lapply(tf_de$Mouse, function(x) {
  filter(x, FC_purity < 0.6) %>% qtl_filter(var = "Count_DE")
})

# RUNX1-TFPI example of high DE count (16/32 exps), avg_avs_FC (0.88), 
# DEPR (0.98), and low purity (0.55). Similar for Mecp2-Plxdc1 (altho lower FC)

filter(pur_lo_hg$RUNX1, Symbol == "TFPI")
filter(pur_lo_mm$Mecp2, Symbol == "Plxdc1")


# Correlation of DE counts and purity - filter genes with abundant NAs

cor_de_purity <- unlist(lapply(tf_de, function(species) {
  lapply(species, function(df) {
    df <- na_filter(df)
    cor(df$Count_DE, df$FC_purity, use = "pairwise.complete.obs")
  })
}))


# Min cor is Human PAX6 and max is Mouse Mecp2

cor_de_purity[which.min(abs(cor_de_purity))]
cor_de_purity[which.max(abs(cor_de_purity))]


# Mouse Neurod1 example of two genes with high count DE and purity, but one is 
# not consistent across GoF/LoF and so signed purity == -1). 

filter(tf_de$Mouse$Neurod1, Symbol == "Sord")
filter(tf_de$Mouse$Neurod1, Symbol == "Sft2d1")


# Find low-intermediate purity (0-0.8) is common for genes with Count_DE > 2 

pur_summ_hg <- lapply(tf_de$Human, function(x) {
  filter(x, Count_DE > 2) %>% 
    pull(FC_purity) %>% 
    summary()
})


pur_summ_mm <- lapply(tf_de$Mouse, function(x) {
  filter(x, Count_DE > 2) %>% 
    pull(FC_purity) %>% 
    summary()
})


# Inspect divergence of FC direction from contextually similar experiments
# Looking at GSE80745 b/c two types of anti-Tcf4 shRNA in same cell line
# ------------------------------------------------------------------------------


context_df <- data.frame(
  Exp1 = mlist_mm$FC_mat[, "GSE80745_Tcf4_Mouse_Knockdown"],
  Exp2 = mlist_mm$FC_mat[, "GSE80745_Tcf4_Mouse_Knockdown-1"]) %>% 
  mutate(Diff = Exp2 - Exp1) %>% 
  arrange(desc(abs(Diff)))

cor(context_df$Exp1, context_df$Exp2, use = "pairwise.complete.obs")
# plot(context_df$Exp1, context_df$Exp2)


# Save out

saveRDS(all_de, file = all_outfile)
saveRDS(tf_de, file = tf_outfile)


# Plots
#-------------------------------------------------------------------------------


# Hists for count DE + proportion measured across all experiments


hist_plot <- function(all_df, var, colour, fdr, title) {
  
  stopifnot(var %in% c("Count_DE", "Proportion_DE_measured"))
  
  if (var == "Count_DE") {
    xlab <- paste0("Count DE (FDR < ", fdr, ")")
    bins <- max(all_df$Count_DE + 1)
  } else {
    xlab <- paste0("Proportion DE of measured (FDR < ", fdr, ")")
    bins <- 51
  }

  vline <- median(all_df[, var])
  
  ggplot(all_df, aes(x = !!sym(var))) +
    geom_histogram(bins = bins, fill = colour) +
    geom_vline(xintercept = vline, col = "red") +
    theme_classic() +
    ylab("Count") +
    xlab(xlab) +
    ggtitle(title) +
    theme(axis.text = element_text(size = 25),
          axis.title = element_text(size = 25),
          plot.title = element_text(size = 25))
  
}


p1a <- hist_plot(all_df = all_de$Human, 
                 var = "Count_DE", 
                 colour = "royalblue", 
                 fdr = fdr,
                 title = "Human")

p1b <- hist_plot(all_df = all_de$Human, 
                 var = "Proportion_DE_measured", 
                 colour = "royalblue", 
                 fdr = fdr,
                 title = "Human")

p1c <- hist_plot(all_df = all_de$Mouse, 
                 var = "Count_DE", 
                 colour = "goldenrod", 
                 fdr = fdr,
                 title = "Mouse")

p1d <- hist_plot(all_df = all_de$Mouse, 
                 var = "Proportion_DE_measured", 
                 colour = "goldenrod", 
                 fdr = fdr,
                 title = "Mouse")


p1 <- plot_grid(p1a, p1b, p1c, p1d, nrow = 2)

ggsave(p1, height = 14, width = 20, dpi = 300, device = "png",
       filename = paste0(plot_dir, "hist_countde_propde_FDR=", fdr, "_", date, ".png"))



# plot relationship between total counts and DE prior ranking


de_scatter <- function(all_df, fdr, title) {
  
  filter(all_df, !is.na(DE_prior_rank)) %>% 
    ggplot(., aes(x = Count_DE, y = DE_prior_rank)) +
    geom_bin2d(bins = max(all_df$Count_DE + 1)) +
    geom_smooth(method = "lm", col = "black", se = FALSE) + 
    theme_classic() +
    ylab("DE prior rank") +
    xlab(paste0("Count DE (FDR < ", fdr, ")")) +
    ggtitle(title) +
    scale_fill_gradient(low = "lightgrey", high = "black", name = "Count per bin") +
    theme(axis.text = element_text(size = 25),
          axis.title = element_text(size = 25),
          plot.title = element_text(size = 25),
          legend.text = element_text(size = 25),
          legend.title = element_text(size = 25))
}


# Human
p2a <- de_scatter(all_df = all_de$Human, fdr = fdr, title = "Human")
p2a_noleg <- p2a + theme(legend.position = "none")

ggsave(p2a, height = 8, width = 10, dpi = 300, device = "png",
       filename = paste0(plot_dir, "prior_decount_bin_human_FDR=", fdr, "_", date, ".png"))

ggsave(p2a_noleg, height = 8, width = 10, dpi = 300, device = "png",
       filename = paste0(plot_dir, "prior_decount_bin_human_nolegend_FDR=", fdr, "_", date, ".png"))


# Mouse
p3a <- de_scatter(all_df = all_de$Mouse, fdr = fdr, title = "Mouse")
p3a_noleg <- p3a + theme(legend.position = "none")

ggsave(p3a, height = 8, width = 10, dpi = 300, device = "png",
       filename = paste0(plot_dir, "prior_decount_bin_mouse_FDR=", fdr, "_", date, ".png"))

ggsave(p3a_noleg, height = 8, width = 10, dpi = 300, device = "png",
       filename = paste0(plot_dir, "prior_decount_bin_mouse_nolegend_FDR=", fdr, "_", date, ".png"))


# Example of human genes with highest Count DE in relation to DE prior. Scatter
# of DE prior ~ Count DE, with top genes highlighted, and strip chart of FCs 
# with DE status paneled beside.


top_genes <- all_de$Human %>% 
  arrange(desc(Count_DE)) %>%
  slice_head(n = 2) %>% 
  pull(Symbol)


p4a <- all_de$Human %>% 
  mutate(Top_count = Symbol %in% top_genes) %>% 
  ggplot() +
  geom_point(aes(x = Count_DE, y = DE_prior_rank), 
             data = . %>% filter(Top_count),
             fill = "blue", size = 3.5, shape = 21) +
  geom_jitter(aes(x = Count_DE, y = DE_prior_rank), 
              data = . %>% filter(!Top_count),
              shape = 21, size = 1, alpha = 0.4, width = 0.1, height = 0.1) +
  geom_text_repel(aes(x = Count_DE, y = DE_prior_rank, label = Symbol),
                  data = . %>% filter(Top_count),
                  force = 0.5, force_pull = 0.5, size = 5) +
  theme_classic() +
  xlab("Count DE FDR < 0.1") +
  ylab("DE prior rank") +
  ggtitle(paste0("Human n=", nrow(meta_hg), " experiments")) +
  theme_classic() +
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size = 25),
        plot.title = element_text(size = 25))



top_list <- lapply(top_genes, function(x) {
  data.frame(
    FC = mlist_hg$FC_mat[x, ],
    DE = mlist_hg$FDR_mat[x, ] < fdr,
    TR = meta_hg$Symbol)
})
names(top_list) <- top_genes


p4_list <- lapply(names(top_list), function(x) {
  
  ggplot() +
    geom_jitter(aes(x = TR, y = FC, fill = DE),
                data = filter(top_list[[x]], !DE),
                shape = 21, size = 3, width = 0.1, height = 0.1) +
    geom_jitter(aes(x = TR, y = FC, fill = DE),
                data = filter(top_list[[x]], DE),
                shape = 21, size = 3, width = 0.1, height = 0.1) +
    scale_fill_manual(values = c("white", "red")) +
    ggtitle(x) +
    ylab("Log2 Fold Change") +
    theme_classic() +
    theme(axis.text = element_text(size = 15),
          axis.title = element_text(size = 20),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), 
          axis.title.x = element_blank(),
          plot.title = element_text(size = 20, hjust = 0.5))
})


p4b <- plot_grid(plotlist = p4_list, ncol = 1)

p4 <- plot_grid(p4a, p4b, ncol = 2, rel_widths = c(1.5, 1), scale = 0.9)

ggsave(p4, height = 9, width = 14, dpi = 300, device = "png", bg = "white",
       filename = paste0(plot_dir, "human_max_countde_vs_deprior_", date, ".png"))


# TF-specific barcharts of count DE


plot_tf_hist <- function(tf_list, meta, species, tf_pal) {
  
  tfs <- names(tf_list)
  
  p_l <- lapply(tfs, function(tf) {
    
    plot_df <- tf_list[[tf]]
    n_tf <- nrow(filter(meta, Species == species & Symbol == tf))
    bins <- max(plot_df$Count_DE + 1)
    
    ggplot(plot_df, aes(x = Count_DE)) +
      geom_histogram(bins = bins, fill = tf_pal[tf]) +
      theme_classic() +
      ylab("Count of genes") +
      xlab("Count DE") +
      ggtitle(paste0(tf, " n=", n_tf)) +
      scale_x_continuous(breaks = pretty_breaks) +
      theme(
        axis.text = element_text(size = 30),
        axis.title = element_text(size = 30),
        plot.title = element_text(size = 35),
        plot.margin = margin(10, 15, 10, 10))
  })
  names(p_l) <- tfs
  
  return(p_l)
}


# Plot ASCL1 separately for example in figure, and combine rest of TRs

p5a_list <- plot_tf_hist(tf_list = tf_de$Human, meta = meta, species = "Human", tf_pal = tf_pal_hg)
p5b_list <- plot_tf_hist(tf_list = tf_de$Mouse, meta = meta, species = "Mouse", tf_pal = tf_pal_mm)

p5c <- plot_grid(plot_grid(plotlist = p5a_list[-1], nrow = 1), 
                 plot_grid(plotlist = p5b_list[-1], nrow = 1), 
                 nrow = 2)

ggsave(p5a_list$ASCL1, dpi = 300, device = "png", height = 6, width = 6,
       filename = paste0(plot_dir, "Human_ASCL1_hist_decounts_FDR=", fdr, "_", date, ".png"))

ggsave(p5b_list$Ascl1, dpi = 300, device = "png", height = 6, width = 6,
       filename = paste0(plot_dir, "Mouse_Ascl1_hist_decounts_FDR=", fdr, "_", date, ".png"))

ggsave(p5c, dpi = 300, device = "png", height = 12, width = 40,
       filename = paste0(plot_dir, "TF_hist_decounts_FDR=", fdr, "_", date, ".png"))


# Boxplot Purity ~ DE Counts, with rect colour indicating hi purity (>0.8)


purity_de_bplot <- function(tf_list, 
                            species,
                            na_filter = 2,
                            hi_purity = 0.8) {
  
  tfs <- names(tf_list)
  
  plist <- lapply(tfs, function(x) {
    
    plot_df <- tf_list[[x]] %>% 
      filter(Count_NA <= na_filter)
    
    ggplot(plot_df, aes(x = as.factor(Count_DE), y = FC_purity)) +
      geom_rect(xmin = 0, 
                xmax = Inf, 
                ymin = 0, 
                ymax = hi_purity, 
                fill = "lightgrey",
                alpha = 0.3) +
      geom_rect(xmin = 0, 
                xmax = Inf, 
                ymin = hi_purity, 
                ymax = Inf, 
                fill = "forestgreen",
                alpha = 0.05) +
      geom_boxplot(width = 0.4) +
      xlab("Count DE FDR < 0.1") +
      ggtitle(paste(species, x)) +
      theme_classic() +
      theme(axis.text = element_text(size = 15),
            axis.title = element_text(size = 25),
            plot.title = element_text(size = 25))
  })
  names(plist) <- tfs
  return(plist)
}


pur_bplot_hg <- purity_de_bplot(tf_list = tf_de$Human, species = "Human")
pur_bplot_mm <- purity_de_bplot(tf_list = tf_de$Mouse, species = "Mouse")


# only want TRs with at least 8 experiments (arbitrary minimum to reduce small
# sample effect sizes on purity)

min_tr <- meta %>% 
  filter(Perturbation != "Mutant") %>% 
  count(Symbol) %>% 
  filter(n >= 8)

min_hg_plot <- pur_bplot_hg[intersect(min_tr$Symbol, names(pur_bplot_hg))]
min_mm_plot <- pur_bplot_mm[intersect(min_tr$Symbol, names(pur_bplot_mm))]

p6 <- plot_grid(plotlist = c(min_hg_plot, min_mm_plot), ncol = 2)

ggsave(p6,
       dpi = 300, device = "png", height = 20, width = 18, bg = "white",
       filename = paste0(plot_dir, "Purity_vs_decount_bplot_mincount_FDR=", fdr, "_", date, ".png"))


# Example heatmap of genes on opposite extremes of Purity


p7_df <- filter(tf_de$Human$ASCL1, Symbol %in% c("CABP7", "PDGFRA", "PCLAF"))

p7 <- fc_heatmap(tf = "ASCL1", tf_df = p7_df, fc_mat = mlist_hg$FC_mat, meta = meta)

png(width = 7, height = 6, units = "in", res = 300,
    filename = paste0(plot_dir, "Human_ASCL1_example_purity_heatmap_", date, ".png"))

p7

graphics.off()


# Example heatmaps of fold changes + DE prior (+/- binarized)

# Human

tf_hg <- "ASCL1"

# Sort by final perturbation rankings
p8_df <- tf_de$Human[[tf_hg]] %>% 
  arrange(desc(Count_DE), desc(Avg_abs_FC)) %>%
  slice_head(n = 15)

p8_list <- list(
  fc = fc_heatmap(tf = tf_hg, tf_df = p8_df, fc_mat = mlist_hg$FC_mat, meta = meta),
  depr = depr_heatmap(tf_df = p8_df, deprior = pdeg_hg),
  binary = depr_heatmap(tf_df = p8_df, deprior = pdeg_hg, deprior_binary = TRUE)
)


png(width = 6, height = 6, units = "in", res = 300,
    filename = paste0(plot_dir, "Human_", tf_hg, "_FC_heatmap_", date, ".png"))
p8_list[[1]]
graphics.off()

png(width = 6, height = 6, units = "in", res = 300,
    filename = paste0(plot_dir, "Human_", tf_hg, "_DEprior_heatmap_", date, ".png"))
p8_list[[2]]
graphics.off()

png(width = 6, height = 6, units = "in", res = 300,
    filename = paste0(plot_dir, "Human_", tf_hg, "_DEprior_binary_heatmap_", date, ".png"))
p8_list[[3]]
graphics.off()

# Mouse

tf_mm <- "Ascl1"

p9_df <- tf_de$Mouse[[tf_mm]] %>% 
  arrange(desc(Count_DE), desc(Avg_abs_FC)) %>% 
  slice_head(n = 15)

p9_list <- list(
  fc = fc_heatmap(tf = tf_mm, tf_df = p9_df, fc_mat = mlist_mm$FC_mat, meta = meta),
  depr = depr_heatmap(tf_df = p9_df, deprior = pdeg_mm),
  binary = depr_heatmap(tf_df = p9_df, deprior = pdeg_mm, deprior_binary = TRUE)
)


png(width = 7, height = 6, units = "in", res = 300,
    filename = paste0(plot_dir, "Mouse_", tf_mm, "_FC_heatmap_", date, ".png"))
p9_list[[1]]
graphics.off()

png(width = 7, height = 6, units = "in", res = 300,
    filename = paste0(plot_dir, "Mouse_", tf_mm, "_DEprior_heatmap_", date, ".png"))
p9_list[[2]]
graphics.off()

png(width = 7, height = 6, units = "in", res = 300,
    filename = paste0(plot_dir, "Mouse_", tf_mm, "_DEprior_binary_heatmap_", date, ".png"))
p9_list[[3]]
graphics.off()
