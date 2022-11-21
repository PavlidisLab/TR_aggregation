## Script that explores counts of differentially expressed genes (DEGs)
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


# Functions
# ------------------------------------------------------------------------------

# Filter df for var in top quantile (qtl)
qtl_filter <- function(df, var, qtl = 0.9) {
  filter(df, !!sym(var) >= quantile(!!sym(var), qtl, na.rm = TRUE))
}


# Heuristic - remove genes in the top quantile (qtl) of count NAs
na_filter <- function(df, qtl = 0.25) {
  filter(df, Count_NA <= quantile(df$Count_NA, qtl))
}


# For all experiments and for each TF, get a data frame of 
# 1) DE counts/fraction measured/NA counts; 
# 2) The count of times a gene was up/down in each of GoF and LoF experiments; 
# 3) Whether the counts were consistent in direction between LoF and GoFs;
# 4) The GoF/LoF class purity; 
# 5) The average absolute FC
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

# post-hoc addition of species-specific DE counts
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

cor.test(
  all_de$Human[!is.na(all_de$Human$DE_Prior_Rank), "DE_Prior_Rank"],
  all_de$Human[!is.na(all_de$Human$DE_Prior_Rank), "Count_DE"]
)

cor.test(
  all_de$Mouse[!is.na(all_de$Mouse$DE_Prior_Rank), "DE_Prior_Rank"],
  all_de$Mouse[!is.na(all_de$Mouse$DE_Prior_Rank), "Count_DE"]
)


# Genes with max DE counts across all experiments

de_max_hg <- filter(all_de$Human, Count_DE == max(all_de$Human$Count_DE))$Symbol
de_max_mm <- filter(all_de$Mouse, Count_DE == max(all_de$Mouse$Count_DE))$Symbol

max_df_hg <- 
  do.call(rbind, lapply(tf_de$Human, function(x) filter(x, Symbol %in% de_max_hg)))

max_df_mm <- 
  do.call(rbind, lapply(tf_de$Mouse, function(x) filter(x, Symbol %in% de_max_mm)))


# Inspect and provide examples of high DE count and low DE prior
# Eg, things that are frequently DE in this collection, but not in the global
# DE prior

prior_lo_hg <- all_de$Human %>%
  filter(Count_DE > 20) %>%
  arrange(DE_Prior_Rank, desc(Count_DE))

prior_lo_mm <- all_de$Mouse %>%
  filter(Count_DE > 30) %>%
  arrange(DE_Prior_Rank, desc(Count_DE))


low_hg <- c("NPEPPS", "ATE1", "AAGAB", "INPPL1")
low_mm <- c("Irak1", "Pofut2")

low_hg_df <- 
  do.call(rbind, lapply(tf_de$Human, function(x) filter(x, Symbol %in% low_hg))) %>% 
  arrange(Symbol)

low_mm_df <- 
  do.call(rbind, lapply(tf_de$Mouse, function(x) filter(x, Symbol %in% low_mm))) %>% 
  arrange(Symbol)


# Inspect relationship of DE prior and average absolute FC

# NOTE: In general see weak cor, but in some instances (eg mouse Mecp2) see weak 
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
fc_hi_hg <- lapply(tf_de$Human, function(x) {
  qtl_filter(x, var = "Count_DE") %>% qtl_filter(var = "Avg_abs_FC")
})

fc_hi_mm <- lapply(tf_de$Mouse, function(x) {
  qtl_filter(x, var = "Count_DE") %>% qtl_filter(var = "Avg_abs_FC")
})


# NOTE: High FC and high DE count tend to have extremely high DE prior. Example
# of low is MEF2C-CLBA1 (DEPR = 0.09, DE 2/4).
lapply(fc_hi_hg, function(x) summary(x$DE_Prior_Rank))
lapply(fc_hi_mm, function(x) summary(x$DE_Prior_Rank))
filter(fc_hi_hg$MEF2C, DE_Prior_Rank == min(DE_Prior_Rank, na.rm = TRUE))


# Examples of genes with high abs FC and low/no DE count across exps - largely
# genes with many NAs
tf_de$Mouse$Runx1 %>% 
  arrange(desc(Avg_abs_FC), Count_DE) %>% 
  head(10)


# Inspect purity
# ------------------------------------------------------------------------------


# Genes with high purity (GoF/LoF gives same FC direction). Remove genes with
# frequent NAs


pur_hi_hg <- lapply(tf_de$Human, function(x) {
  na_filter(x) %>% 
  qtl_filter(var = "Purity", qtl = 0.9) %>%
  arrange(desc(Avg_abs_FC))
})


pur_hi_mm <- lapply(tf_de$Mouse, function(x) {
  na_filter(x) %>% 
  qtl_filter(var = "Purity", qtl = 0.9) %>%
  arrange(desc(Avg_abs_FC))
})


# Inspect genes with high DE count and low purity
# NOTE: RUNX1-TFPI example of high DE count (16/32 exps), avg_avs_FC (0.88), 
# DEPR (0.98), and low purity (0.55). Similar for Mecp2-Plxdc1 (altho lower FC)

pur_lo_hg <- lapply(tf_de$Human, function(x) {
  filter(x, Purity < 0.6) %>% qtl_filter(var = "Count_DE")
})

pur_lo_mm <- lapply(tf_de$Mouse, function(x) {
  filter(x, Purity < 0.6) %>% qtl_filter(var = "Count_DE")
})

filter(pur_lo_hg$RUNX1, Symbol == "TFPI")
filter(pur_lo_mm$Mecp2, Symbol == "Plxdc1")


# Correlation of DE counts and purity - filter genes with abundant NAs
# Min cor is Human Hes1 (has minimal DE count) and max is Mouse Mecp2

cor_de_purity <- unlist(lapply(tf_de, function(species) {
  lapply(species, function(df) {
    df <- na_filter(df)
    cor(df$Count_DE, df$Purity, use = "pairwise.complete.obs")
  })
}))

cor_de_purity[which.min(abs(cor_de_purity))]
cor_de_purity[which.max(abs(cor_de_purity))]


# Mouse Neurod1 example of two genes with high count DE and purity, but not
# consistent (so signed purity == 1 and -1). But individual FCs are pretty weak

filter(tf_de$Mouse$Neurod1, Symbol == "Sord")
filter(tf_de$Mouse$Neurod1, Symbol == "Sft2d1")


# Find low-intermediate purity (0-0.8) is common for genes with Count_DE > 2 

pur_summ_hg <- lapply(tf_de$Human, function(x) {
  filter(x, Count_DE > 2) %>% 
    pull(Purity) %>% 
    summary()
})


pur_summ_mm <- lapply(tf_de$Mouse, function(x) {
  filter(x, Count_DE > 2) %>% 
    pull(Purity) %>% 
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
plot(context_df$Exp1, context_df$Exp2)


# Save out

saveRDS(all_de, file = tf_outfile)
saveRDS(tf_de, file = tf_outfile)


# Plots
#-------------------------------------------------------------------------------


# Hists for count DE + fraction measured

p1a <- 
  ggplot(all_de$Human, aes(x = Count_DE)) +
  geom_histogram(bins = max(all_de$Human$Count_DE + 1), fill = "royalblue") +
  geom_vline(xintercept = median(all_de$Human$Count_DE), col = "red") +
  theme_classic() +
  ylab("Count") +
  xlab(paste0("Count DE (FDR < ", fdr, ")")) +
  ggtitle("Human") +
  theme(axis.text = element_text(size = 25),
        axis.title = element_text(size = 25),
        plot.title = element_text(size = 25))


p1b <- 
  ggplot(all_de$Human, aes(x = Fraction_DE_measured)) +
  geom_histogram(bins = 51, fill = "royalblue") +
  geom_vline(xintercept = median(all_de$Human$Fraction_DE_measured), col = "red") +
  theme_classic() +
  ylab("Count") +
  xlab(paste0("Fraction DE of measured (FDR < ", fdr, ")")) +
  ggtitle("Human") +
  theme(axis.text = element_text(size = 25),
        axis.title = element_text(size = 25),
        plot.title = element_text(size = 25))

p1c <- 
  ggplot(all_de$Mouse, aes(x = Count_DE)) +
  geom_histogram(bins = max(all_de$Mouse$Count_DE+1), fill = "goldenrod") +
  geom_vline(xintercept = median(all_de$Mouse$Count_DE), col = "red") +
  theme_classic() +
  ylab("Count") +
  xlab(paste0("Count DE (FDR < ", fdr, ")")) +
  ggtitle("Mouse") +
  theme(axis.text = element_text(size = 25),
        axis.title = element_text(size = 25),
        plot.title = element_text(size = 25))


p1d <- 
  ggplot(all_de$Mouse, aes(x = Fraction_DE_measured)) +
  geom_histogram(bins = 51, fill = "goldenrod") +
  geom_vline(xintercept = median(all_de$Mouse$Fraction_DE_measured), col = "red") +
  theme_classic() +
  ylab("Count") +
  xlab(paste0("Fraction DE of measured (FDR < ", fdr, ")")) +
  ggtitle("Mouse") +
  theme(axis.text = element_text(size = 25),
        axis.title = element_text(size = 25),
        plot.title = element_text(size = 25))


p1 <- plot_grid(p1a, p1b, p1c, p1d, nrow = 2)

ggsave(p1, height = 14, width = 20, dpi = 300, device = "png",
       filename = paste0(plot_dir, "hist_countde_fracde_FDR=", fdr, "_", date, ".png"))


# plot relationship between total counts and DE prior ranking


# Human

p2a <- 
  filter(all_de$Human, !is.na(DE_Prior_Rank)) %>% 
  ggplot(., aes(x = Count_DE, y = DE_Prior_Rank)) +
  geom_jitter(shape = 19, alpha = 0.2) +
  theme_classic() +
  theme(axis.text = element_text(size = 25),
        axis.title = element_text(size = 25))

p2a <- ggMarginal(p2a, type = "histogram", fill = "royalblue")

p2b <- 
  filter(all_de$Human, !is.na(DE_Prior_Rank)) %>% 
  ggplot(., aes(x = Count_DE, y = DE_Prior_Rank)) +
  geom_bin2d(bins = max(all_de$Human$Count_DE+1)) +
  geom_smooth(method = "lm", col = "black", se = FALSE) + 
  theme_classic() +
  ylab("DE prior rank") +
  xlab(paste0("Count DE (FDR < ", fdr, ")")) +
  ggtitle("Human") +
  scale_fill_gradient(low = "lightgrey", high = "black", name = "Count per bin") +
  theme(axis.text = element_text(size = 25),
        axis.title = element_text(size = 25),
        plot.title = element_text(size = 25),
        legend.text = element_text(size = 25),
        legend.title = element_text(size = 25))

p2b_noleg <- p2b + theme(legend.position = "none")


ggsave(p2b, height = 8, width = 10, dpi = 300, device = "png",
       filename = paste0(plot_dir, "prior_decount_bin_human_FDR=", fdr, "_", date, ".png"))

ggsave(p2b_noleg, height = 8, width = 10, dpi = 300, device = "png",
       filename = paste0(plot_dir, "prior_decount_bin_human_nolegend_FDR=", fdr, "_", date, ".png"))


# mouse

p3a <- 
  filter(all_de$Mouse, !is.na(DE_Prior_Rank)) %>% 
  ggplot(., aes(x = Count_DE, y = DE_Prior_Rank)) +
  geom_jitter(shape = 19, alpha = 0.2) +
  theme_classic() +
  theme(axis.text = element_text(size = 25),
        axis.title = element_text(size = 25))

p3a <- ggMarginal(p3a, type = "histogram", fill = "goldenrod")

p3b <- 
  filter(all_de$Mouse, !is.na(DE_Prior_Rank)) %>% 
  ggplot(., aes(x = Count_DE, y = DE_Prior_Rank)) +
  geom_bin2d(bins = max(all_de$Mouse$Count_DE)+1) +
  geom_smooth(method = "lm", col = "black", se = FALSE) + 
  theme_classic() +
  ylab("DE prior rank") +
  xlab(paste0("Count DE (FDR < ", fdr, ")")) +
  ggtitle("Mouse") +
  # scale_fill_gradientn(colors = hex_colors(11), name = "Count per bin") +
  scale_fill_gradient(low = "lightgrey", high = "black", name = "Count per bin") +
  theme(axis.text = element_text(size = 25),
        axis.title = element_text(size = 25),
        plot.title = element_text(size = 25),
        legend.text = element_text(size = 25),
        legend.title = element_text(size = 25))

p3b_noleg <- p3b + theme(legend.position = "none")

ggsave(p3b, height = 8, width = 10, dpi = 300, device = "png",
       filename = paste0(plot_dir, "prior_decount_bin_mouse_FDR=", fdr, "_", date, ".png"))

ggsave(p3b_noleg, height = 8, width = 10, dpi = 300, device = "png",
       filename = paste0(plot_dir, "prior_decount_bin_mouse_nolegend_FDR=", fdr, "_", date, ".png"))


# Example of human genes with highest Count DE in relation to DE prior


p4a <- all_de$Human %>% 
  mutate(Max_count = Count_DE == max(Count_DE)) %>% 
  ggplot() +
  geom_point(aes(x = Count_DE, y = DE_Prior_Rank), 
             data = . %>% filter(Max_count),
             fill = "blue", size = 3.5, shape = 21) +
  geom_jitter(aes(x = Count_DE, y = DE_Prior_Rank), 
              data = . %>% filter(!Max_count),
              shape = 21, size = 1, alpha = 0.4, width = 0.1, height = 0.1) +
  geom_text_repel(aes(x = Count_DE, y = DE_Prior_Rank, label = Symbol),
                  data = . %>% filter(Max_count),
                  force = 0.5, force_pull = 0.5, size = 5) +
  theme_classic() +
  xlab("Count DE FDR < 0.1") +
  ylab("DE prior rank") +
  ggtitle(paste0("Human n=", nrow(meta_hg), " experiments")) +
  theme_classic() +
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size = 25),
        plot.title = element_text(size = 25))


max_count <- filter(all_de$Human, Count_DE == max(Count_DE))


max_list <- lapply(max_count$Symbol, function(x) {
  data.frame(
    FC = mlist_hg$FC_mat[x, ],
    DE = mlist_hg$FDR_mat[x, ] < fdr,
    TR = meta_hg$Symbol)
})
names(max_list) <- max_count$Symbol


p4_list <- lapply(names(max_list), function(x) {
  
  ggplot() +
    geom_jitter(aes(x = TR, y = FC, fill = DE),
                data = filter(max_list[[x]], !DE),
                shape = 21, size = 3, width = 0.1, height = 0.1) +
    geom_jitter(aes(x = TR, y = FC, fill = DE),
                data = filter(max_list[[x]], DE),
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


ggsave(p4, height = 9, width = 14, dpi = 300, device = "png",
       filename = paste0(plot_dir, "human_max_countde_vs_deprior_", date, ".png"))


# TF-specific barcharts of count DE


plot_tf_hist <- function(tf_list, meta, species) {
  
  tfs <- names(tf_list)
  
  for (tf in tfs) {
    
    plot_df <- tf_list[[tf]]
    n_tf <- nrow(filter(meta, Species == species & Symbol == tf))
    bins <- max(plot_df$Count_DE+1)
    
    p <- 
      ggplot(plot_df, aes(x = Count_DE)) +
      geom_histogram(bins = bins, fill = tf_pal[str_to_title(tf)]) +
      theme_classic() +
      ylab("Count of genes") +
      xlab(paste0("Count DE (FDR < ", fdr, ")")) +
      ggtitle(paste0(tf, " n=", n_tf)) +
      scale_x_continuous(breaks = pretty_breaks) +
      theme(axis.text = element_text(size = 35),
            axis.title = element_text(size = 30),
            plot.title = element_text(size = 35))

    ggsave(p, dpi = 300, device = "png", height = 6, width = 6,
           filename = paste0(plot_dir, species, "_", tf, "_decounts_FDR=", fdr, "_", date, ".png"))
    
  }
  
}


# Running saves individual hist plots for each TF to plot dir
plot_tf_hist(tf_list = tf_de$Human, meta = meta, species = "Human")
plot_tf_hist(tf_list = tf_de$Mouse, meta = meta, species = "Mouse")


# Scatter of Purity ~ Count DE 


purity_scatter <- function(tf_list,
                           species,
                           purity = "Purity",  # Purity | Signed_purity
                           na_filter = 2,
                           hi_purity = 0.8) {
  
  stopifnot(purity %in% c("Purity", "Signed_purity"))
  
  tfs <- names(tf_list)
  
  plist <- lapply(tfs, function(x) {
    
    plot_df <- tf_list[[x]] %>% 
      filter(Count_NA <= na_filter)
    
    ggplot(plot_df, aes(x = Count_DE, y = !!sym(purity))) +
      geom_rect(xmin = -Inf, 
                xmax = Inf, 
                ymin = -Inf, 
                ymax = hi_purity, 
                fill = "lightgrey",
                alpha = 0.3) +
      geom_rect(xmin = -Inf, 
                xmax = Inf, 
                ymin = hi_purity, 
                ymax = Inf, 
                fill = "forestgreen",
                alpha = 0.05) +
      geom_jitter(size = 1.8, shape = 21, colour = "black") +
      theme_classic() +
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


pur_scatter_hg <- purity_scatter(tf_list = tf_de$Human, purity = "Purity", species = "Human")
sign_pur_scatter_hg <- purity_scatter(tf_list = tf_de$Human, purity = "Signed_purity", species = "Human")

pur_scatter_mm <- purity_scatter(tf_list = tf_de$Mouse, purity = "Purity", species = "Mouse")
sign_pur_scatter_mm <- purity_scatter(tf_list = tf_de$Mouse, purity = "Signed_purity", species = "Mouse")


ggsave(plot_grid(plotlist = pur_scatter_hg, ncol = 4),
       dpi = 300, device = "png", height = 12, width = 24,
       filename = paste0(plot_dir, "Purity_vs_decount_scatter_human_all_FDR=", fdr, "_", date, ".png"))

ggsave(plot_grid(plotlist = sign_pur_scatter_hg, ncol = 4),
       dpi = 300, device = "png", height = 12, width = 24,
       filename = paste0(plot_dir, "Sign_purity_vs_decount_scatter_human_all_FDR=", fdr, "_", date, ".png"))


ggsave(plot_grid(plotlist = pur_scatter_mm, ncol = 4),
       dpi = 300, device = "png", height = 12, width = 24,
       filename = paste0(plot_dir, "Purity_vs_decount_scatter_mouse_all_FDR=", fdr, "_", date, ".png"))

ggsave(plot_grid(plotlist = sign_pur_scatter_mm, ncol = 4),
       dpi = 300, device = "png", height = 12, width = 24,
       filename = paste0(plot_dir, "Sign_purity_vs_decount_scatter_mouse_all_FDR=", fdr, "_", date, ".png"))


# Boxplot Purity ~ DE Counts, with rect colour indicating hi purity (>0.8)


purity_de_bplot <- function(tf_list, 
                            species,
                            na_filter = 2,
                            hi_purity = 0.8) {
  
  tfs <- names(tf_list)
  
  plist <- lapply(tfs, function(x) {
    
    plot_df <- tf_list[[x]] %>% 
      filter(Count_NA <= na_filter)
    
    ggplot(plot_df, aes(x = as.factor(Count_DE), y = Purity)) +
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

ggsave(plot_grid(plotlist = pur_bplot_hg, ncol = 4),
       dpi = 300, device = "png", height = 12, width = 24,
       filename = paste0(plot_dir, "Purity_vs_decount_bplot_human_all_FDR=", fdr, "_", date, ".png"))

ggsave(plot_grid(plotlist = pur_bplot_mm, ncol = 4),
       dpi = 300, device = "png", height = 12, width = 24,
       filename = paste0(plot_dir, "Purity_vs_decount_bplot_mouse_all_FDR=", fdr, "_", date, ".png"))

ggsave(plot_grid(plotlist = c(min_hg_plot, min_mm_plot), ncol = 2),
       dpi = 300, device = "png", height = 20, width = 18,
       filename = paste0(plot_dir, "Purity_vs_decount_bplot_mincount_FDR=", fdr, "_", date, ".png"))



# Histograms of purity


pur_distn <- function(tf_list,
                      species,
                      na_filter = 2,
                      hi_purity = 0.8) {
  
  tfs <- names(tf_list)
  
  plist <- lapply(tfs, function(x) {
    
    plot_df <- tf_list[[x]] %>% 
      filter(Count_NA <= na_filter)
    
    ggplot(plot_df, aes(x = Purity)) +
      geom_rect(xmin = 0, 
                xmax = hi_purity, 
                ymin = -Inf, 
                ymax = Inf, 
                fill = "lightgrey",
                alpha = 0.3) +
      geom_rect(xmin = hi_purity, 
                xmax = Inf, 
                ymin = -Inf, 
                ymax = Inf, 
                fill = "forestgreen",
                alpha = 0.05) +
      geom_histogram(bins = 20, colour = "black") +
      ggtitle(paste(species, x)) +
      theme_classic() +
      theme(axis.text = element_text(size = 15),
            axis.title = element_text(size = 25),
            plot.title = element_text(size = 25))
  })
  names(plist) <- tfs
  return(plist)
}


pur_distn_hg <- pur_distn(tf_list = tf_de$Human, species = "Human")
pur_distn_mm <- pur_distn(tf_list = tf_de$Mouse, species = "Mouse")


ggsave(plot_grid(plotlist = pur_distn_hg, ncol = 4),
       dpi = 300, device = "png", height = 12, width = 24,
       filename = paste0(plot_dir, "Purity_distn_human_all_FDR=", fdr, "_", date, ".png"))

ggsave(plot_grid(plotlist = pur_distn_mm, ncol = 4),
       dpi = 300, device = "png", height = 12, width = 24,
       filename = paste0(plot_dir, "Purity_distn_mouse_all_FDR=", fdr, "_", date, ".png"))


# Example heatmap of genes on opposite extremes of Purity


plot_df <- filter(tf_de$Human$ASCL1, Symbol %in% c("CABP7", "PDGFRA", "PCLAF"))

p5 <- fc_heatmap(tf = "ASCL1", tf_df = plot_df, fc_mat = mlist_hg$FC_mat, meta = meta)

png(width = 7, height = 6, units = "in", res = 300,
    filename = paste0(plot_dir, "Human_ASCL1_example_purity_heatmap_", date, ".png"))
p5
graphics.off()


# Example heatmaps of fold changes + DE prior (+/- binarized)


tf_hg <- "ASCL1"


# Sort by final perturbation rankings
df_hg <- tf_de$Human[[tf_hg]] %>% 
  arrange(desc(Count_DE), desc(Avg_abs_FC)) %>%
  slice_head(n = 15)


list_hg <- list(
  fc = fc_heatmap(tf = tf_hg, tf_df = df_hg, fc_mat = mlist_hg$FC_mat, meta = meta),
  depr = depr_heatmap(tf_df = df_hg, deprior = pdeg_hg),
  binary = depr_heatmap(tf_df = df_hg, deprior = pdeg_hg, deprior_binary = TRUE)
)


png(width = 6, height = 6, units = "in", res = 300,
    filename = paste0(plot_dir, "Human_", tf_hg, "_FC_heatmap_", date, ".png"))
list_hg[[1]]
graphics.off()

png(width = 6, height = 6, units = "in", res = 300,
    filename = paste0(plot_dir, "Human_", tf_hg, "_DEprior_heatmap_", date, ".png"))
list_hg[[2]]
graphics.off()

png(width = 6, height = 6, units = "in", res = 300,
    filename = paste0(plot_dir, "Human_", tf_hg, "_DEprior_binary_heatmap_", date, ".png"))
list_hg[[3]]
graphics.off()


tf_mm <- "Ascl1"


df_mm <- tf_de$Mouse[[tf_mm]] %>% 
  arrange(desc(Count_DE), desc(Avg_abs_FC)) %>% 
  slice_head(n = 15)

list_mm <- list(
  fc = fc_heatmap(tf = tf_mm, tf_df = df_mm, fc_mat = mlist_mm$FC_mat, meta = meta),
  depr = depr_heatmap(tf_df = df_mm, deprior = pdeg_mm),
  binary = depr_heatmap(tf_df = df_mm, deprior = pdeg_mm, deprior_binary = TRUE)
)


png(width = 7, height = 6, units = "in", res = 300,
    filename = paste0(plot_dir, "Mouse_", tf_mm, "_FC_heatmap_", date, ".png"))
list_mm[[1]]
graphics.off()

png(width = 7, height = 6, units = "in", res = 300,
    filename = paste0(plot_dir, "Mouse_", tf_mm, "_DEprior_heatmap_", date, ".png"))
list_mm[[2]]
graphics.off()

png(width = 7, height = 6, units = "in", res = 300,
    filename = paste0(plot_dir, "Mouse_", tf_mm, "_DEprior_binary_heatmap_", date, ".png"))
list_mm[[3]]
graphics.off()
