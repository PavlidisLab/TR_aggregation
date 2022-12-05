## This script explores the perturb metadata and gene x experiment matrices.
## Summarizes and plots the count of genes represented across experiments and 
## array vs RNA-seq experiments.
## -----------------------------------------------------------------------------

library(tidyverse)
library(pheatmap)
library(RColorBrewer)
source("R/setup-01_config.R")
source("R/utils/plot_functions.R")

plot_dir <- paste0(pplot_dir, "Meta_sample_matrix/")

# Load meta and lists of perturb effect size matrices
meta <- read.delim(file =  paste0(meta_dir, "batch1_tfperturb_meta_final_", date, ".tsv"), stringsAsFactors = FALSE)
mlist_hg <- readRDS(paste0(pmat_dir, "human_list_perturb_matrix_", date, ".RDS"))
mlist_mm <- readRDS(paste0(pmat_dir, "mouse_list_perturb_matrix_", date, ".RDS"))
mlist_ortho <- readRDS(paste0(pmat_dir, "ortho_list_perturb_matrix_", date, ".RDS"))

# Expression platform info
platform_meta <- read.delim(platform_path, stringsAsFactors = FALSE)


# Functions
#-------------------------------------------------------------------------------

# Helper to get the gene/row-wise sum of non-NAs for a perturb matrix

count_gene_measured <- function(mat) {
  apply(mat, 1, function(x) sum(!is.na(x)))
}


# Helper to get the col/experiment-wise sum of non-NAs for a perturb matrix

count_gene_coverage <- function(mat) {
  apply(mat, 2, function(x) sum(!is.na(x)))
}


# Describe meta
# ------------------------------------------------------------------------------


# distinct GEO IDs: multiple experiments can come from the same GEO accession,
# and must consider instances where Gemma splits an accession (which produces
# a digit suffix). NOTE - won't capture when same study is split over GSEs
n_gse <- n_distinct(str_replace(meta$GSE, "\\.[:digit:]", ""))

# Count by species
n_species <- table(meta$Species)

# Count TR by species
n_species_tr <- table(str_to_title(meta$Symbol), meta$Species)

# By perturbation type
n_pert <- table(meta$Perturbation)
n_species_pert <- table(meta$Perturbation, meta$Species)
prop_pert <- table(meta$Perturbation)/nrow(meta)

# Inspect perturbation notes
table(meta$Perturbation_Note, meta$Perturbation)

# Inpsect cell types
sort(table(meta$Cell_Type), decreasing = TRUE)


# Looking at gene coverage across experiments assembled in the perturb matrices
#-------------------------------------------------------------------------------


mlist_all  <- list(
  Human = mlist_hg$Tstat_mat,  # tstat/fc/prfc equivalent for NA detection
  Mouse = mlist_mm$Tstat_mat,
  Ortholog = mlist_ortho$Tstat_mat
)


# count of experiments for each perturb matrix
exp_counts <- unlist(lapply(mlist_all, ncol))


# how many times was each gene measured across the experiments
gene_counts <- lapply(mlist_all, count_gene_measured)


# count of total measured genes for each experiment - summarize by species
exp_gene_counts <- lapply(mlist_all, count_gene_coverage)
lapply(exp_gene_counts, summary)


# get the genes that were measured at least once
measured <- lapply(mlist_all, nrow)


# genes that were measured in every experiment
all_measured <- lapply(1:length(gene_counts), function(x) {
  gcounts <- gene_counts[[x]]
  ecounts <- exp_counts[[x]]
  length(gcounts[gcounts == ecounts])
})
names(all_measured) <- names(mlist_all)


# Platform curation - array vs sequencing using Nathaniel's platform metadata
#-------------------------------------------------------------------------------


platforms <- str_extract(meta$Platform, "^GPL[:digit:]+")
identical(length(platforms), nrow(meta))

tech <- unlist(lapply(platforms, function(x) {
  if (!(x %in% platform_meta$ad.Name)) {
    return(NA)
  } else {
    platform_meta[platform_meta$ad.Name == x, "ad.TechType"]
  }
}))

if (any(is.na(tech))) warning("Not all GPL identifiers were matched")


meta$Platform_type <- ifelse(tech == "SEQUENCING", "Sequencing", "Array")

# breakdown of experiments that were sequencing or array

tech_count <- table(meta$Platform_type)
tech_frac <- table(meta$Platform_type)/nrow(meta)


# use GPL identifiers to count unique platforms per species

unique_gpl_all <- meta %>%
  mutate(Platform = unlist(str_extract(Platform, "GPL[:digit:]+"))) %>% 
  group_by(Platform_type) %>% 
  summarise(n = n_distinct(Platform))


# example of sequencing experiment that does not have full coverage - pipeline
# removes lowly expressed genes

seq_exp <- mlist_mm$Tstat_mat[, "GSE122990_Neurod1_Mouse_Knockout"]
sum(is.na(seq_exp))


# microarray platforms with most gene coverage

exp_count_df <- data.frame(
  n_gene = c(exp_gene_counts$Human,
             exp_gene_counts$Mouse)) %>% 
  rownames_to_column("Experiment_ID")

micro_count <- meta %>% 
  filter(Platform_type == "Array") %>% 
  dplyr::select(Experiment_ID, Species, Platform) %>% 
  left_join(., y = exp_count_df, by = "Experiment_ID") 


# Plotting
#-------------------------------------------------------------------------------


# Count of experiments as stacked bar
# using same order as frequency of ChIP-seq experiments

tf_order <- rev(c("RUNX1", "MECP2", "ASCL1", "NEUROD1", "MEF2C", "PAX6", "TCF4", "HES1"))

p1a <- meta %>% 
  mutate(Symbol = str_to_upper(Symbol)) %>% 
  count(Symbol, Species) %>% 
  mutate(Symbol = factor(Symbol, levels = tf_order)) %>% 
  ggplot(., aes(x = Symbol, y = n, fill = Species)) +
  geom_bar(stat = "identity", colour = "black", width = 0.8) +
  ylim(0, 100) +
  ylab("Count of experiments") +
  theme_classic() +
  scale_fill_manual(values = c("royalblue", "goldenrod")) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 30),
        axis.text.y = element_text(size = 30),
        axis.text.x = element_text(size = 30, angle = 60, vjust = 1, hjust = 1),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15),
        legend.position = "bottom")
        
# Coord flip - ultimately used in paper

p1b <- p1a + 
  coord_flip() +
  ylab("Count of experiments") +
  theme(axis.title.y = element_blank(),
        axis.title.x = element_text(size = 30),
        legend.position = "none")
  
ggsave(p1b, dpi = 300, device = "png", height = 10, width = 12, 
       file = paste0(plot_dir, "batch1_tf_perturb_counts_", date, ".png"))


# bar chart of sequencing technology and perturbation type

p2a <- data.frame(tech_count) %>% 
  ggplot(., aes(x = Var1, y = Freq)) +
  geom_bar(stat = "identity", width = 0.8, col = "black", fill = "grey") +
  theme_classic() +
  ylab("Count of experiments") +
  ylim(c(0, 150)) +
  theme(axis.title.x = element_blank(),
        axis.text = element_text(size = 25),
        axis.title = element_text(size = 30))

ggsave(p2a, dpi = 300, device = "png", height = 8, width = 6, 
       file = paste0(plot_dir, "batch1_tf_perturb_tech_counts_", date, ".png"))


p2b <- meta %>% 
  count(Perturbation) %>% 
  ggplot(., aes(x = Perturbation, y = n)) +
  geom_bar(stat = "identity", width = 0.8, col = "black", fill = "grey") +
  theme_classic() +
  ylab("Count of experiments") +
  ylim(c(0, 150)) +
  theme(axis.title.x = element_blank(),
        axis.text = element_text(size = 25),
        axis.title = element_text(size = 30),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

ggsave(p2b, dpi = 300, device = "png", height = 8, width = 6, 
       file = paste0(plot_dir, "batch1_tf_perturb_type_counts_", date, ".png"))



# Histograms of experiment gene coverage


p3a <- 
  data.frame(Count_measured = exp_gene_counts$Human) %>% 
  ggplot(., aes(x = Count_measured)) + 
  geom_histogram(bins = 20) +
  ggtitle("Human experiment gene coverage (n=77)") +
  ylab("Count of experiments") +
  xlab("Count of genes") +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 5)) +
  theme_classic() +
  theme(axis.title = element_text(size = 25),
        axis.text = element_text(size = 25),
        title = element_text(size = 25))


p3b <- 
  data.frame(Count_measured = exp_gene_counts$Mouse) %>% 
  ggplot(., aes(x = Count_measured)) + 
  geom_histogram(bins = 30) +
  ggtitle("Mouse experiment gene coverage (n=165)") +
  ylab("Count of experiments") +
  xlab("Count of genes") +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 5)) +
  theme_classic() +
  theme(axis.title = element_text(size = 25),
        axis.text = element_text(size = 25),
        title = element_text(size = 25))


p3c <- 
  data.frame(Count_measured = exp_gene_counts$Ortholog) %>% 
  ggplot(., aes(x = Count_measured)) + 
  geom_histogram(bins = 30) +
  ggtitle("All experiment gene coverage (n=242)") +
  ylab("Count of experiments") +
  xlab("Count of genes") +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 5)) +
  theme_classic() +
  theme(axis.title = element_text(size = 25),
        axis.text = element_text(size = 25),
        title = element_text(size = 25))


ggsave(p3a, height = 9, width = 14, dpi = 300, device = "png",
       filename = paste0(plot_dir, "Gene_measure_count_human_", date, ".png"))

ggsave(p3b, height = 9, width = 14, dpi = 300, device = "png",
       filename = paste0(plot_dir, "Gene_measure_count_mouse_", date, ".png"))

ggsave(p3c, height = 9, width = 14, dpi = 300, device = "png",
       filename = paste0(plot_dir, "Gene_measure_count_ortho_", date, ".png"))


# heatmaps to show gene coverage across experiments

anno_col_hg <- list(Symbol = tf_pal_hg,
                    Platform_type = c(Array = "lightgrey", Sequencing = "grey30"))

anno_col_mm <- list(Symbol = tf_pal_mm,
                    Platform_type = c(Array = "lightgrey", Sequencing = "grey30"))


# mouse binary (gene present or not)

mat_mm_binary <- mlist_mm$Tstat_mat
mat_mm_binary[!is.na(mat_mm_binary)] <- 1
mat_mm_binary[is.na(mat_mm_binary)] <- 0


meta_anno_mm <- meta %>% 
  filter(Species == "Mouse") %>% 
  select(Symbol, Platform_type)

rownames(meta_anno_mm) <- colnames(mat_mm_binary)


# default cluster - demonstrate genes that are consistently measured or not


pheatmap(
  t(mat_mm_binary),
  col = c("black", "royalblue"),
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = FALSE,
  show_colnames = FALSE,
  annotation_row = meta_anno_mm,
  annotation_names_row = FALSE,
  annotation_colors = anno_col_mm,
  legend = FALSE,
  height = 6,
  width = 10,
  treeheight_row = 0,
  treeheight_col = 0,
  filename = paste0(plot_dir, "mouse_sample_matrix_clustered_", date, ".png")
)
  
  
# human binary
  
mat_hg_binary <- mlist_hg$Tstat_mat
mat_hg_binary[!is.na(mat_hg_binary)] <- 1
mat_hg_binary[is.na(mat_hg_binary)] <- 0
  
meta_anno_hg <- meta %>%
  filter(Species == "Human") %>%
  select(Symbol, Platform_type)

rownames(meta_anno_hg) <- colnames(mat_hg_binary)
  
# default cluster - demonstrate genes that are consistently measured or not
  
pheatmap(
  t(mat_hg_binary),
  col = c("black", "royalblue"),
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = FALSE,
  show_colnames = FALSE,
  annotation_row = meta_anno_hg,
  annotation_names_row = FALSE,
  annotation_colors = anno_col_hg,
  legend = FALSE,
  height = 6,
  width = 10,
  treeheight_row = 0,
  treeheight_col = 0,
  filename = paste0(plot_dir, "human_sample_matrix_clustered_", date, ".png")
)
  
  
# ortholog binary
  
mat_ortho_binary <- mlist_all$Ortholog
mat_ortho_binary[!is.na(mat_ortho_binary)] <- 1
mat_ortho_binary[is.na(mat_ortho_binary)] <- 0
  
meta_anno_all <- meta %>%
  select(Symbol, Platform_type) %>%
  mutate(Symbol = str_to_upper(Symbol))
  
rownames(meta_anno_all) <- colnames(mat_ortho_binary)
  
pheatmap(
  t(mat_ortho_binary),
  col = c("black", "royalblue"),
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = FALSE,
  show_colnames = FALSE,
  annotation_row = meta_anno_all,
  annotation_names_row = FALSE,
  annotation_colors = anno_col_hg,
  legend = FALSE,
  height = 7,
  width = 10,
  treeheight_row = 0,
  treeheight_col = 0,
  filename = paste0(plot_dir, "all_sample_matrix_clustered_", date, ".png")
)
