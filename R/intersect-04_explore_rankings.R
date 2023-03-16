## This script is for the interactive exploration/plotting of the assembled
## gene rankings. It additionally saves out scatter plots of the gene targets.
## -----------------------------------------------------------------------------

library(tidyverse)
library(WGCNA)
library(ggrepel)
source("R/setup-01_config.R")
source("R/utils/plot_functions.R")

topn <- 500  # number of top genes to consider
plot_dir <- file.path(iplot_dir, "Gene_rankings/")

# List of all TR rankings and data matrices
rank_list <- readRDS(rank_path)
dat_list <- readRDS(alldat_path)

# Curated targets
lt_all <- read.delim(curated_path_all, stringsAsFactors = FALSE)

# Mapping of orthologous genes
pc_ortho <- read.delim(ortho_path, stringsAsFactors = FALSE)
tfs_hg <- names(rank_list$Human)
tfs_mm <- names(rank_list$Mouse)

# Isolate meta for plotting
bmeta_hg <- filter(dat_list$Binding$Meta, Species == "Human") %>% arrange(Symbol)
bmeta_mm <- filter(dat_list$Binding$Meta, Species == "Mouse") %>% arrange(Symbol)
pmeta_hg <- filter(dat_list$Perturbation$Meta, Species == "Human") %>% arrange(Symbol)
pmeta_mm <- filter(dat_list$Perturbation$Meta, Species == "Mouse") %>% arrange(Symbol)


# Correlate aggregate scores
# ------------------------------------------------------------------------------


get_cor <- function(mat, method = "spearman") {
  WGCNA::cor(select_if(mat, is.numeric), 
             use = "pairwise.complete.obs", 
             method = method)
}

keep_cols <- c("Mean_bind", "Count_DE", "Rank_binding", "Rank_perturbation", "Rank_integrated")

cor_list <- list(
  Human = lapply(rank_list$Human, function(x) get_cor(x[, keep_cols])),
  Mouse = lapply(rank_list$Mouse, function(x) get_cor(x[, keep_cols])),
  Ortho = lapply(rank_list$Ortho, function(x) get_cor(x[, keep_cols]))
)


# Genes that are topn in both species - the ortho rankings can prioritize 
# genes driven by one species. This specifically looks for topn ranked in both.
# ------------------------------------------------------------------------------


top_ortho <- lapply(tfs_hg, function(x) {
  
  # Pulling species rank from top ortho rank. Note that some genes may be 
  # filtered from individual sets (eg, no binding) so must get consensus
  
  human_rank <- rank_list$Human[[x]][1:topn, ] %>% 
    dplyr::select(Symbol, Rank_integrated) %>% 
    dplyr::rename(Rank_hg = Rank_integrated)
  
  mouse_rank <- rank_list$Mouse[[str_to_title(x)]][1:topn, ] %>% 
    dplyr::select(Symbol, Rank_integrated) %>% 
    dplyr::rename(Rank_mm = Rank_integrated)
  
  ortho_rank <- rank_list$Ortho[[x]][1:topn, ] %>% 
    dplyr::select(Symbol, Rank_integrated) %>% 
    dplyr::rename(Rank_ortho = Rank_integrated) %>% 
    left_join(pc_ortho, by = c("Symbol" = "ID")) %>% 
    left_join(human_rank, by = c("Symbol_hg" = "Symbol")) %>% 
    left_join(mouse_rank, by = c("Symbol_mm" = "Symbol")) %>% 
    filter(Rank_hg <= topn, Rank_mm <= topn)
  
  return(ortho_rank)
  
})
names(top_ortho) <- tfs_hg


# TCF4 shares 29 ortho genes in top 500, RUNX1 shares 110
n_top_ortho <- unlist(lapply(top_ortho, nrow))


# The following is for interactive viewing of TR rankings
# ------------------------------------------------------------------------------


tf <- "NEUROD1"

# Subset to genes that are topn ranked in both data types, as rank product can
# be driven by one of the rankings
topn_mm <- filter(rank_list$Mouse[[str_to_title(tf)]], Rank_binding <= topn & Rank_perturbation <= topn)
topn_hg <- filter(rank_list$Human[[str_to_upper(tf)]], Rank_binding <= topn & Rank_perturbation <= topn)
topn_ortho <- filter(rank_list$Ortho[[str_to_upper(tf)]], Rank_binding <= topn & Rank_perturbation <= topn)

# view(top_mm)
# view(top_hg)
# view(top_ortho[[tf]])
# view(rank_list$Mouse[[tf]])
# view(rank_list$Human[[str_to_upper(tf)]])
# view(rank_list$Ortho[[str_to_upper(tf)]])


# Ranks of notable genes highlighted in paper
# ASCL1: c("DLL1", "DLL3", "DLL4", "HES6", "JAG2", "LFNG", "CDC25B", "ID1", "ID3", "ZBTB18", "CBFA2T3", "KCNH2", "KRTAP9-3", "SHB", "BMP7")
# HES1: c("LTB", "HES1", "ATOH1", "STARD7", "E2F5", "PFN1", "BAHCC1", "FBXO31")
# MECP2: c("PCDHGA7", "ESRRG", "SLC6A7", "NRG2", "SDK1", "AUTS2", "BDNF", "AFF1", "IRAK1", "CDS1", "SHANK2", "PLXNA2", "CNTNAP2", "TENM2", "TENM3")
# MEF2C: c("MEF2C", "NR4A1", "HDAC5", "HDAC9", "MEF2D", "KLF6", "KLF2", "KLF4", "ARID1A")
# NEUROD1: c("INSM1", "SRRM4", "COG1", "CXXC4", "STXBP1", "NOVA2", "PTPRK", "UGCG", "TRIM9", "SFT2D2", "SLC35F1", "ADGRL3", "SOGA1")
# PAX6: c("PAX6", "MAB21L1", "MAB21L2", "RLBP1", "EPHA3", "DNAJB6", "BLVRA", "MAP2", "ZNF608", "ABHD4")
# RUNX1: c("NFE2", "ITGB2", "MYH9", "PLEC", "LMNA", "VIM", "PXN", "SPI1", "RUNX1", "MOB3A", "PHF19")
# TCF4: c("TCF4", "CDKN1A", "DEPP1", "MFAP4", "ZBTB18", "GPSM2", "MC4R", "TLCD1")
# ------------------------------------------------------------------------------


# Isolate integrated rank for mouse in human of requested tf/genes into a df

get_ranks <- function(tf, gene_vec, rank_list) {
  
  gene_ranks <- lapply(gene_vec, function(x) {
    
    hg <- filter(rank_list$Human[[str_to_upper(tf)]], Symbol == str_to_upper(x))$Rank_integrated
    mm <- filter(rank_list$Mouse[[str_to_title(tf)]], Symbol == str_to_title(x))$Rank_integrated
    
    data.frame(
      Symbol = x,
      Human = ifelse(length(hg > 0), hg, NA),
      Mouse = ifelse(length(hg > 0), mm, NA)
    )
  })
  
  data.frame(do.call(rbind, gene_ranks))
}


tf <- "HES1"
gene_vec <- c("LTB", "HES1", "ATOH1", "STARD7", "E2F5", "PFN1", "BAHCC1", "FBXO31")
gene_ranks <- get_ranks(tf, gene_vec, rank_list)


# The following is for subsetting the effect sizes for a given gene for
# interactive plotting
# ------------------------------------------------------------------------------


# Return a list of the binding and perturb effect sizes for requested TR + gene

get_tf_gene <- function(tf, gene, dat_list, fdr = 0.1) {
  
  # isolate minimal metadata for tf
  
  pmeta <- dat_list$Perturbation$Meta %>% 
    filter(str_to_lower(Symbol) == str_to_lower(tf)) %>% 
    dplyr::select(Experiment_ID, Cell_Type, Species)
  
  bmeta <- dat_list$Binding$Meta %>% 
    filter(str_to_lower(Symbol) == str_to_lower(tf)) %>% 
    dplyr::select(Experiment_ID, Cell_Type, Species)
  
  # experiment IDs by species 
  
  phg <- intersect(colnames(dat_list$Perturbation$Human$FC_mat), pmeta$Experiment_ID)
  pmm <- intersect(colnames(dat_list$Perturbation$Mouse$FC_mat), pmeta$Experiment_ID)
  bhg <- intersect(colnames(dat_list$Binding$Human$QN_log), bmeta$Experiment_ID)
  bmm <- intersect(colnames(dat_list$Binding$Mouse$QN_log), bmeta$Experiment_ID)
  
  # perturbation data
  
  perturb <- data.frame(
    FC = c(dat_list$Perturbation$Human$FC_mat[str_to_upper(gene), phg],
           dat_list$Perturbation$Mouse$FC_mat[str_to_title(gene), pmm]),
    DE = c(dat_list$Perturbation$Human$FDR_mat[str_to_upper(gene), phg] < fdr,
           dat_list$Perturbation$Mouse$FDR_mat[str_to_title(gene), pmm] < fdr)
  ) %>% 
    rownames_to_column(var = "Experiment_ID") %>% 
    left_join(., pmeta, by = "Experiment_ID")
  
 # binding data
  
  bind <- data.frame(
    Binding_score = c(dat_list$Binding$Human$QN_log[str_to_upper(gene), bhg],
                      dat_list$Binding$Mouse$QN_log[str_to_title(gene), bmm])
  ) %>%
    rownames_to_column(var = "Experiment_ID") %>% 
    left_join(., bmeta, by = "Experiment_ID")
  
  return(list(Perturb = perturb, Bind = bind))
  
}


# Return a list of gene stats over all experiments, split by species

get_gene <- function(gene, dat_list, fdr = 0.1) {
  
  bind_hg <- data.frame(
    Bind = dat_list$Binding$Human$QN_log[str_to_upper(gene), ])
  
  bind_mm <- data.frame(
    Bind = dat_list$Binding$Mouse$QN_log[str_to_title(gene), ])
  
  perturb_hg <- data.frame(
    FC = dat_list$Perturbation$Human$FC_mat[str_to_upper(gene), ],
    DE = dat_list$Perturbation$Human$FDR_mat[str_to_upper(gene), ] < fdr)
  
  perturb_mm <- data.frame(
    FC = dat_list$Perturbation$Mouse$FC_mat[str_to_title(gene), ],
    DE = dat_list$Perturbation$Mouse$FDR_mat[str_to_title(gene), ] < fdr)
  
  return(list(Bind_hg = bind_hg,
              Bind_mm = bind_mm,
              Perturb_hg = perturb_hg, 
              Perturb_mm = perturb_mm))
  
}


tf <- "Ascl1"
gene <- "Dll1"
tf_gene_stats <- get_tf_gene(tf, gene, dat_list)
all_gene_stats <- get_gene(gene, dat_list)


# plotting within TF group

plot(tf_gene_stats$Perturb$FC, 
     col = ifelse(tf_gene_stats$Perturb$DE, "red", "black"), 
     main = paste(tf, gene, sep = "-"))

plot(tf_gene_stats$Perturb$FC, 
     col = ifelse(tf_gene_stats$Perturb$Species == "Mouse", "Goldenrod", "RoyalBlue"), 
     main = paste(tf, gene, sep = "-"))

plot(tf_gene_stats$Bind$Bind, 
     col = ifelse(tf_gene_stats$Bind$Species == "Mouse", "Goldenrod", "RoyalBlue"), 
     main = paste(tf, gene, sep = "-"))


# plotting across all experiments


# Distn of human binding

plot(all_gene_stats$Bind_hg$Bind, 
     col = ifelse(bmeta_hg$Symbol == str_to_upper(tf), "RoyalBlue", "black"))

boxplot(all_gene_stats$Bind_hg$Bind ~ bmeta_hg$Symbol)

boxplot(all_gene_stats$Bind_hg$Bind ~ bmeta_hg$Symbol == str_to_upper(tf))


# Distn of mouse binding

plot(all_gene_stats$Bind_mm$Bind, 
     col = ifelse(bmeta_mm$Symbol == str_to_title(tf), "goldenrod", "black"))

boxplot(all_gene_stats$Bind_mm$Bind ~ bmeta_mm$Symbol)

boxplot(all_gene_stats$Bind_mm$Bind ~ bmeta_mm$Symbol == str_to_title(tf))


# Human FC

plot(all_gene_stats$Perturb_hg$FC, 
     col = ifelse(pmeta_hg$Symbol == str_to_upper(tf), "RoyalBlue", "black"))

boxplot(abs(all_gene_stats$Perturb_hg$FC) ~ pmeta_hg$Symbol)

boxplot(abs(all_gene_stats$Perturb_hg$FC) ~ pmeta_hg$Symbol == str_to_upper(tf))


# Mouse FC

plot(all_gene_stats$Perturb_mm$FC, 
     col = ifelse(pmeta_mm$Symbol == str_to_title(tf), "goldenrod", "black"))

boxplot(abs(all_gene_stats$Perturb_mm$FC) ~ pmeta_mm$Symbol)

boxplot(abs(all_gene_stats$Perturb_mm$FC) ~ pmeta_mm$Symbol == str_to_title(tf))


# Look at RUNX1 and gene scores in Macrophage cells
# Microglia development genes  # Csf1r Irf8 Spi1 Mafb
# https://www.annualreviews.org/doi/full/10.1146/annurev-immunol-032713-120240 
# https://www.biorxiv.org/content/10.1101/2021.05.30.446351v1.full.pdf 
# ------------------------------------------------------------------------------


tf <- "Runx1"
gene <- "Spi1"
runx1_tf <- get_tf_gene(tf, gene, dat_list)
runx1_all <- get_gene(gene, dat_list)

mpg <- c(
  "mESC derived macrophages",
  "Macrophages",
  "RAW264.7 cells",
  "Granulocyte-macrophage progenitors"
)

runx1_tf$Perturb$Macrophage <- runx1_tf$Perturb$Cell_Type %in% mpg
runx1_tf$Bind$Macrophage <- runx1_tf$Bind$Cell_Type %in% mpg

wilcox.test(runx1_tf$Bind$Bind ~ runx1_tf$Bind$Macrophage)
plot(runx1_tf$Bind$Bind, col = ifelse(runx1_tf$Bind$Macrophage, "red", "black"))
plot(runx1_tf$Bind$Bind, col = ifelse(runx1_tf$Bind$Species == "Mouse", "Goldenrod", "RoyalBlue"))
boxplot(runx1_tf$Bind$Bind ~ runx1_tf$Bind$Macrophage)


# Scatter plots of mean bind score vs count DE, where curated targets are 
# coloured in blue, while curated targets that are in the topn of the aggregate
# ranking are in gold+labeled
# ------------------------------------------------------------------------------


plot_scatter <- function(df, tf) {
  
  # Create group var for presence in curation +/- topn
  
  df$Group <- vapply(1:nrow(df), function(x) {
    if (!df$Curated_target[x]) {
      return("Out")
    } else if (df$Curated_target[x] & df$Rank_integrated[x] < topn) {
      return("Top")
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
                    force = 1, force_pull = 2, size = 7, max.overlaps = 20) +
    xlab("Mean binding score") +
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
      axis.title = element_text(size = 30),
      plot.title = element_text(size = 30),
      legend.text = element_text(size = 25),
      legend.title = element_text(size = 25),
      plot.margin = margin(10, 10, 10, 10)
    )
}


plist_hg3 <- lapply(tfs_hg, function(x) plot_scatter(df = rank_list$Human[[x]], tf = x))
names(plist_hg3) <- tfs_hg


plist_mm3 <- lapply(tfs_mm, function(x) plot_scatter(df = rank_list$Mouse[[x]], tf = x))
names(plist_mm3) <- tfs_mm


plist_ortho3 <- lapply(tfs_hg, function(x) plot_scatter(df = rank_list$Ortho[[x]], tf = x))
names(plist_ortho3) <- tfs_hg


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


for (tf in tfs_hg) {
  ggsave(plist_ortho3[[tf]],
         dpi = 300, device = "png", height = 8, width = 12,
         filename = paste0(plot_dir, "Ortho_", tf, "_targets.png"))
}



for (tf in tfs_hg) {
  ggsave(plot_grid(plist_hg3[[str_to_upper(tf)]] + theme(legend.position = "none"), 
                   plist_mm3[[str_to_title(tf)]] + theme(legend.position = "none")),
         dpi = 300, device = "png", height = 8, width = 14,
         filename = paste0(plot_dir, "Combined_", tf, "_targets.png"))
}



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
