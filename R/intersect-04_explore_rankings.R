## This script is for the interactive exploration/plotting of the assembled
## gene rankings. It additionally saves out scatter plots of the gene targets
## -----------------------------------------------------------------------------

library(tidyverse)
library(WGCNA)
library(ggrepel)
source("~/regnetR/R/utils/plot_functions.R")

topn <- 500  # number of top genes to consider
fdr <- 0.1
plot_dir <- "~/Plots/Intersect/"
date <- "Apr2022"  # most recent data freeze

# List of all TR rankings and data matrices
rank_list <- readRDS(paste0("~/scratch/R_objects/", date, "_ranked_target_list.RDS"))
dat_list <- readRDS(paste0("~/scratch/R_objects/", date, "_all_data_list.RDS"))

pc_ortho <- read.delim("~/Data/Metadata/hg_mm_1to1_ortho_genes_DIOPT-v8.tsv", stringsAsFactors = FALSE)
tfs_hg <- names(rank_list$Human)
tfs_mm <- names(rank_list$Mouse)


# Correlate aggregate scores
# ------------------------------------------------------------------------------


get_cor <- function(mat, method = "spearman") {
  WGCNA::cor(select_if(mat, is.numeric), use = "pairwise.complete.obs", method = method)
}

keep_cols <- c("Mean_bind", "Count_DE", "Rank_binding", "Rank_perturbation", "Rank_integrated")

cor_list <- list(
  Human = lapply(rank_list$Human, function(x) get_cor(x[, keep_cols])),
  Mouse = lapply(rank_list$Mouse, function(x) get_cor(x[, keep_cols])),
  Ortho = lapply(rank_list$Ortho, function(x) get_cor(x[, keep_cols]))
)


# Genes that are top 500 in both species - the ortho rankings can prioritize 
# genes driven by one species. This specifically looks for topn ranked in both
# ------------------------------------------------------------------------------


top_ortho <- lapply(tfs_mm, function(x) {
  
  # Pulling species rank from top ortho rank. Note that some genes may be 
  # filtered from individual sets (eg, no binding) so must get consensus
  
  ortho_symbol <- filter(pc_ortho, ID %in% rank_list$Ortho[[x]]$Symbol[1:topn])
  
  human_rank <- rank_list$Human[[str_to_upper(x)]] %>% 
    filter(Symbol %in% ortho_symbol$Symbol_hg) %>% 
    arrange(match(Symbol, ortho_symbol$Symbol_hg))
  
  mouse_rank <- rank_list$Mouse[[x]] %>% 
    filter(Symbol %in% ortho_symbol$Symbol_mm) %>% 
    arrange(match(Symbol, ortho_symbol$Symbol_mm))
  
  ortho_rank <- data.frame(Symbol = ortho_symbol$ID,
                           Human = human_rank$Rank_integrated,
                           Mouse = mouse_rank$Rank_integrated
  ) %>%
    filter(Mouse < topn & Human < topn)
  
})
names(top_ortho) <- tfs_mm


# Hes1 and Mecp2 share 34 ortho genes in top 500, Runx1 shares 104
n_top_ortho <- unlist(lapply(top_ortho, nrow))
summary(n_top_ortho)


# Used for interactive viewing of a TR's rankings
# ------------------------------------------------------------------------------


tf <- "Neurod1"
cmeta <- filter(dat_list$Binding$Meta, str_to_title(Symbol) == tf)
pmeta <- filter(dat_list$Perturbation$Meta, str_to_title(Symbol) == tf)

# Subset to genes that are topn ranked in both data types
top_mm <- filter(rank_list$Mouse[[tf]], Rank_binding <= topn & Rank_perturbation <= topn)
top_hg <- filter(rank_list$Human[[str_to_upper(tf)]], Rank_binding <= topn & Rank_perturbation <= topn)
top_ortho_all <- filter(rank_list$Ortho[[str_to_title(tf)]], Rank_binding <= topn & Rank_perturbation <= topn)

# view(top_mm)
# view(top_hg)
# view(top_ortho[[tf]])
# view(rank_list$Mouse[[tf]])
# view(rank_list$Human[[str_to_upper(tf)]])
# view(rank_list$Ortho[[str_to_title(tf)]])


# Used for subsetting and plotting a specific gene
# ------------------------------------------------------------------------------

# Notable genes (highlighted in paper)
# ASCL1: DLL1, DLL3, DLL4, HES6, LFNG, CDC25B, KRTAP9-3, SHB, BMP7, JAG2, ID1, ID3, ZBTB18, CBFA2T3, KCNH2
# HES1: ATOH1, LTB, STARD7, BAHCC1, SOX12, E2F5, PFN1
# MECP2: PCDHGA7, IRAK1, ESRRG, SLC6A7, NRG2, SDK1, AUTS2
# MEF2C: ARID1A, HDAC5, HDAC9, MEF2D
# NEUROD1: PTPRK, UGCG, SFT2D2, SLC35F1, SOGA1, SORD, SFTD1, SRRM4, COG1, CXXC4, STXBP1
# PAX6: EPHA3, DNAJB6, MAB21L1, MAB21L2, ZNF608/Zfp608, ABHD4, BLVRA
# RUNX1: FDFT1, NFE2, PLEC, MYH9, LMNA, VIM, PXN, PHF19, SPI1, MOB3A
# TCF4: DEPP1, MFAP4, ZBTB18, TLCD1, GPSM2, MC4R



gene_list <- function(tf, gene, fdr = 0.1) {
  
  # Return a list of 2 dfs for chip and perturb data for requested gene
  
  pmeta <- dat_list$Perturbation$Meta %>% 
    filter(str_to_lower(Symbol) == str_to_lower(tf))
  
  phg <- intersect(colnames(dat_list$Perturbation$Human$FC_mat), pmeta$Experiment_ID)
  pmm <- intersect(colnames(dat_list$Perturbation$Mouse$FC_mat), pmeta$Experiment_ID)
  
  bmeta <- dat_list$Binding$Meta %>% 
    filter(str_to_lower(Symbol) == str_to_lower(tf))
  
  bhg <- intersect(colnames(dat_list$Binding$Human$QN_log), bmeta$Experiment_ID)
  bmm <- intersect(colnames(dat_list$Binding$Mouse$QN_log), bmeta$Experiment_ID)
  
  perturb <- data.frame(
    FC = c(dat_list$Perturbation$Human$FC_mat[str_to_upper(gene), phg],
           dat_list$Perturbation$Mouse$FC_mat[str_to_title(gene), pmm]),
    DE = c(dat_list$Perturbation$Human$FDR_mat[str_to_upper(gene), phg] < fdr,
           dat_list$Perturbation$Mouse$FDR_mat[str_to_title(gene), pmm] < fdr)
  ) %>% 
    rownames_to_column(var = "Experiment_ID") %>% 
    left_join(., pmeta, by = "Experiment_ID")
  
  
  bind <- data.frame(
    Bind = c(dat_list$Binding$Human$QN_log[str_to_upper(gene), bhg],
             dat_list$Binding$Mouse$QN_log[str_to_title(gene), bmm])
  ) %>%
    rownames_to_column(var = "Experiment_ID") %>% 
    left_join(., bmeta, by = "Experiment_ID")
  
  return(list(Perturb = perturb, Bind = bind))
  
}



tf <- "Ascl1"
gene <- "Dll1"
gene_stats <- gene_list(tf, gene)

# plotting binding data (all TRs)
bmeta_hg <- filter(dat_list$Binding$Meta, Species == "Human") %>% arrange(Symbol)
bgene_hg <- dat_list$Binding$Human$QN_log[str_to_upper(gene), bmeta_hg$Experiment_ID]
bmeta_mm <- filter(dat_list$Binding$Meta, Species == "Mouse") %>% arrange(Symbol)
bgene_mm <- dat_list$Binding$Mouse$QN_log[str_to_title(gene), bmeta_mm$Experiment_ID]

# plotting perturb data (all TRs)
pmeta_hg <- filter(dat_list$Perturbation$Meta, Species == "Human") %>% arrange(Symbol)

pgene_hg <- data.frame(
  FC = dat_list$Perturbation$Human$FC_mat[str_to_upper(gene), pmeta_hg$Experiment_ID],
  DE = dat_list$Perturbation$Human$FDR_mat[str_to_upper(gene), pmeta_hg$Experiment_ID] < fdr)

pmeta_mm <- filter(dat_list$Perturbation$Meta, Species == "Mouse") %>% arrange(Symbol)

pgene_mm <- data.frame(
  FC = dat_list$Perturbation$Mouse$FC_mat[str_to_title(gene), pmeta_mm$Experiment_ID],
  DE = dat_list$Perturbation$Mouse$FDR_mat[str_to_title(gene), pmeta_mm$Experiment_ID] < fdr)

# plotting within TF group
plot(gene_stats$Perturb$FC, col = ifelse(gene_stats$Perturb$DE, "red", "black"), main = paste(tf, gene, sep = "-"))
plot(gene_stats$Bind$Bind, col = ifelse(gene_stats$Bind$Species == "Mouse", "Goldenrod", "RoyalBlue"), main = paste(tf, gene, sep = "-"))

# Distn of human binding
boxplot(bgene_hg ~ bmeta_hg$Symbol)
plot(bgene_hg, col = ifelse(bmeta_hg$Symbol == str_to_upper(tf), "RoyalBlue", "black"))
plot(density(bgene_hg[names(bgene_hg) %in% cmeta$Experiment_ID]), col = "red")
lines(density(bgene_hg[!names(bgene_hg) %in% cmeta$Experiment_ID]))

# Distn of mouse binding
boxplot(bgene_mm ~ bmeta_mm$Symbol)
plot(bgene_mm, col = ifelse(bmeta_mm$Symbol == str_to_title(tf), "red", "black"))
plot(density(bgene_mm[names(bgene_mm) %in% cmeta$Experiment_ID]), col = "red")
lines(density(bgene_mm[!names(bgene_mm) %in% cmeta$Experiment_ID]))

# Human FC
boxplot(abs(pgene_hg$FC) ~ pmeta_hg$Symbol)
plot(pgene_hg$FC, col = ifelse(pmeta_hg$Symbol == str_to_upper(tf), "red", "black"))
plot(density(abs(pgene_hg[rownames(pgene_hg) %in% pmeta$Experiment_ID, "FC"]), na.rm = TRUE), col = "red")
lines(density(abs(pgene_hg[!rownames(pgene_hg) %in% pmeta$Experiment_ID, "FC"]), na.rm = TRUE))

# Mouse FC
boxplot(abs(pgene_mm$FC) ~ pmeta_mm$Symbol)
plot(pgene_mm$FC, col = ifelse(pmeta_mm$Symbol == str_to_title(tf), "red", "black"))
plot(density(abs(pgene_mm[rownames(pgene_mm) %in% pmeta$Experiment_ID, "FC"]), na.rm = TRUE), col = "red")
lines(density(abs(pgene_mm[!rownames(pgene_mm) %in% pmeta$Experiment_ID, "FC"]), na.rm = TRUE))


# Look at RUNX1 and gene scores in Macrophage cells
# Microglia development genes  # Csf1r Irf8 Spi1 Mafb
# https://www.annualreviews.org/doi/full/10.1146/annurev-immunol-032713-120240 
# https://www.biorxiv.org/content/10.1101/2021.05.30.446351v1.full.pdf 

tf <- "Runx1"
gene <- "Spi1"
gene_stats <- gene_list(tf, gene)

mpg <- c(
  "mESC derived macrophages",
  "Macrophages",
  "RAW264.7 cells",
  "Granulocyte-macrophage progenitors"
)

gene_stats$Perturb$Macrophage <- gene_stats$Perturb$Cell_Type %in% mpg
gene_stats$Bind$Macrophage <- gene_stats$Bind$Cell_Type %in% mpg

wilcox.test(gene_stats$Bind$Bind ~ gene_stats$Bind$Macrophage)
plot(gene_stats$Bind$Bind, col = ifelse(gene_stats$Bind$Macrophage, "red", "black"))
plot(gene_stats$Bind$Bind, col = ifelse(gene_stats$Bind$Species == "Mouse", "Goldenrod", "RoyalBlue"))
boxplot(gene_stats$Bind$Bind ~ gene_stats$Bind$Macrophage)



# Plots
# ------------------------------------------------------------------------------


# Scatter plots of mean bind score vs count DE, where curated targets are 
# coloured in blue, while curated targets that are in the topn of the aggregate
# ranking are in gold+labeled


plot_scatter <- function(df, tf) {
  
  # Create group var for presence in curation +/- topn
  
  df$Group <- vapply(1:nrow(df), function(x) {
    if (!df$Curated_target[x]) {
      return("Out")
    } else if (df$Curated_target[x] & df$Rank_integrated[x] < topn) {
      return ("Top")
    } else {
      return("Low-throughput")
    }
  }, FUN.VALUE = character(1))
  df$Group <- factor(df$Group, levels = c("Out", "Low-throughput", "Top"))
  
  
  ggplot() +
    geom_jitter(data = df[df$Group == "Out", ],
                aes(x = Mean_bind, y = Count_DE, fill = Group),
                alpha = 1, shape = 19, colour = "grey", size = 2, height = 0.3, width = 0.01) +
    geom_jitter(data = df[df$Group == "Low-throughput", ],
                aes(x = Mean_bind, y = Count_DE, fill = Group),
                alpha = 1, shape = 21, size = 3, col = "black", height = 0.3, width = 0.01) +
    geom_jitter(data = df[df$Group == "Top", ],
                aes(x = Mean_bind, y = Count_DE, fill = Group),
                alpha = 1, shape = 21, size = 3, col = "black", height = 0.3, width = 0.01) +
    geom_text_repel(data = df[df$Group == "Top", ],
                    aes(x = Mean_bind, y = Count_DE, label = Symbol),
                    force = 0.5, force_pull = 0.5, size = 5, max.overlaps = 20) +
    xlab("Mean bind score") +
    ylab("Count DE (FDR < 0.1)") +
    ggtitle(tf) +
    scale_y_continuous(breaks = pretty_breaks) +
    scale_fill_manual(name = "Low-throughput evidence",
                      values = c("Out" = "grey",
                                 "Low-throughput" = "royalblue",
                                 "Top" = "goldenrod"),
                      labels = c("None",
                                 "Curated",
                                 "Curated and top 500 integrated"),
                      breaks = levels(df$Group)) +
    theme_classic() +
    theme(
      axis.text = element_text(size = 25),
      axis.title = element_text(size = 25),
      plot.title = element_text(size = 25),
      legend.text = element_text(size = 20),
      legend.title = element_text(size = 20),
      plot.margin = margin(10, 10, 10, 10)
    )
}


plist_hg3 <- lapply(tfs_hg, function(x) plot_scatter(df = rank_list$Human[[x]], tf = x))
names(plist_hg3) <- tfs_hg


plist_mm3 <- lapply(tfs_mm, function(x) plot_scatter(df = rank_list$Mouse[[x]], tf = x))
names(plist_mm3) <- tfs_mm


for (tf in tfs_hg) {
  ggsave(plist_hg3[[tf]],
         dpi = 300, device = "png", height = 8, width = 12,
         filename = paste0(plot_dir, "Human_", tf, "_targets.png"))
}



# ASCL1 no legend for figure
ggsave(plist_hg3$ASCL1 + theme(legend.position = "none"),
       dpi = 300, device = "png", height = 8, width = 8,
       filename = paste0(plot_dir, "Human_ASCL1_targets_noleg.png"))



for (tf in tfs_mm) {
  ggsave(plist_mm3[[tf]],
         dpi = 300, device = "png", height = 8, width = 12,
         filename = paste0(plot_dir, "Mouse_", tf, "_targets.png"))
}

# Pax6 no legend for figure
ggsave(plist_mm3$Pax6 + theme(legend.position = "none"),
       dpi = 300, device = "png", height = 8, width = 8,
       filename = paste0(plot_dir, "Mouse_Pax6_targets_noleg.png"))



# Heatmap of the correlations between the aggregate scores


cor_df <- data.frame(
  Human = do.call(rbind, lapply(cor_list$Human, function(x) x["Rank_binding", "Rank_perturbation"])),
  Mouse = do.call(rbind, lapply(cor_list$Mouse, function(x) x["Rank_binding", "Rank_perturbation"])),
  Ortho = do.call(rbind, lapply(cor_list$Ortho, function(x) x["Rank_binding", "Rank_perturbation"]))
)
cor_df <- round(cor_df, 3)


pal_length <- 7
heatmap_pal <- c('#fff7fb','#ece2f0','#d0d1e6','#a6bddb','#67a9cf','#3690c0','#02818a')
color_min <- 0
color_max <- max(cor_df, na.rm = TRUE)
color_breaks <- seq(color_min, color_max, length.out = pal_length)

pheatmap(cor_df,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         color = heatmap_pal,
         display_numbers = TRUE,
         number_color = "black",
         fontsize = 20,
         border_color = "black",
         na_col = "white",
         cellwidth = 50,
         cellheight = 50,
         width = 5,
         height = 7,
         file = paste0(plot_dir, "countDE_meanbind_cor_heatmap_", date, ".png"))
