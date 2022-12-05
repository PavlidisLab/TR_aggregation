## This script tests each TR's gene rankings for status in the curated low-throughput 
## resource, and performs a precision recall analysis to determine how well the
## respective rankings recover curated targets.
## DREAM5 (PR analysis for ranking regulation): https://pubmed.ncbi.nlm.nih.gov/22796662/
## DREAM2: (Appendix details of PR): https://pubmed.ncbi.nlm.nih.gov/19348640/
## TODO: finalize plotting
## TODO: which targets are missing from each species
## TODO: currently just human/mouse, ortho needs symbol match in sample_target_auprc
## -----------------------------------------------------------------------------

library(tidyverse)
library(parallel)
library(pheatmap)
library(DescTools)
library(cowplot)
source("R/setup-01_config.R")
source("R/utils/plot_functions.R")
source("R/utils/ranking_functions.R")

plot_dir <- paste0(iplot_dir, "Gene_rankings/")

# List of all TR rankings and data matrices
rank_list <- readRDS(paste0(scratch_dir, date, "_ranked_target_list.RDS"))
dat_list <- readRDS(paste0(scratch_dir, date, "_all_data_list.RDS"))

# Curated targets
lt_all <- read.delim(paste0(meta_dir, "Curated_targets_all_July2022.tsv"), stringsAsFactors = FALSE)

pc_ortho <- read.delim(paste0(meta_dir, "hg_mm_1to1_ortho_genes_DIOPT-v8.tsv"), stringsAsFactors = FALSE)
tfs_hg <- names(rank_list$Human)
tfs_mm <- names(rank_list$Mouse)

ol_list <- readRDS(paste0(scratch_dir, date, "_intersect_similarity.RDS"))

# Describe count of targets and human/mouse presence
# TODO: human/mouse presence
# ------------------------------------------------------------------------------


lt_sub <- lt_all %>% 
  mutate(
    TF_Symbol = str_to_upper(TF_Symbol),
    Target_Symbol = str_to_upper(Target_Symbol),
    Experiment_Type = ifelse(is.na(Experiment_Type), "Unknown", Experiment_Type)
    ) %>% 
  filter(TF_Symbol %in% names(rank_list$Human)) 


# Count by low-throughpuut table
lt_sub %>% 
  group_by(TF_Symbol) %>% 
  summarise(n = n_distinct(Target_Symbol))


# Using human rankings to get the number of curated targets
n_target <- data.frame(
  Symbol = names(rank_list$Human),
  Count = unlist(lapply(rank_list$Human, function(x) sum(x$Curated_target)))
) %>% 
  arrange(desc(Count))

# How many unique targets across the TRs
n_unique <- n_distinct(unlist(
  lapply(rank_list$Human, function(x) filter(x, Curated_target)$Symbol)
))


# TODO: better spot for this note
# NOTE: TF-targets can have multiplicity in evidence for the same experiment
# type, given that they are done in unique cell types. Eg, ASCL1-DLL1 has 7 
# entries in human: 3x perturbation, 2x binding, 2x reporter assays, and 5
# in mouse: 2 binding, 2 perturbation, 1 reporter. However, also note that 9/12
# of these experiments come from the same study. RUNX1-SPI1 21 experiments


n_type <- lt_sub %>% 
  mutate(TF_Symbol = factor(TF_Symbol, levels = unique(n_target$Symbol))) %>% 
  count(Experiment_Type, TF_Symbol, name = "Count")


exp_per_target <- lapply(split(lt_sub, lt_sub$TF_Symbol), function(x) {
  
  a <- count(x, Target_Symbol, Experiment_Type)
  b <- matrix(0, nrow = n_distinct(a$Target_Symbol), ncol = 5)
  rownames(b) <- unique(a$Target_Symbol)
  colnames(b) <- c(unique(n_type$Experiment_Type), "Sum")
  
  for (i in 1:nrow(a)) {
    b[a$Target_Symbol[i], a$Experiment_Type[i]] <- a$n[i]
  }
  
  b[, "Sum"] <- rowSums(b)
  
  b <- b[order(b[, "Sum"], decreasing = TRUE), ]
  
  return(b)
})
names(exp_per_target) <- unique(lt_sub$TF_Symbol)



# Wilcoxon rank sum test for the three gene rankings by presence of curated
# evidence. # In mouse find evidence for all rankings save for Neurod1 
# perturbation and Tcf4 binding and integrated. In human it's a bit more spotty 
# (no MEF2C reaches sig) but the majority of comparisons still provide evidence
# ------------------------------------------------------------------------------


wilx_test <- function(df) {
  c(Perturbation = wilcox.test(df$Rank_perturbation ~ df$Curated_target)$p.value,
    Binding = wilcox.test(df$Rank_binding ~ df$Curated_target)$p.value,
    Integrated = wilcox.test(df$Rank_integrated ~ df$Curated_target)$p.value)
}


wilx_list <- 
  lapply(rank_list, function(x) signif(do.call(rbind, lapply(x, wilx_test)), 3))


wilx_sig <- lapply(wilx_list, function(x) x < 0.05)


# The following generates for each TR and each gene ranking, a precision recall 
# data frame where each entry has the calculated P+R for the top k ranked genes
# presence in the curated resource. These are then summarized with the area
# under the precision recall curve (AUPRC)
# ------------------------------------------------------------------------------


pr_list <- list(
  Human = lapply(rank_list$Human, get_all_pr),
  Mouse = lapply(rank_list$Mouse, get_all_pr),
  Ortho = lapply(rank_list$Ortho, get_all_pr)
)


auprc_list <- list(
  Human = do.call(rbind, lapply(rank_list$Human, get_all_auprc)),
  Mouse = do.call(rbind, lapply(rank_list$Mouse, get_all_auprc)),
  Ortho = do.call(rbind, lapply(rank_list$Ortho, get_all_auprc))
)


# Tally which rank was most performant

max_auprc <- lapply(auprc_list, function(x) {
  colnames(x)[apply(x, 1, which.max)]
})

tally_auprc <- lapply(max_auprc, table)


# Inspecting the top 5 curated targets for each rank


top5 <- function(df, rank) {
  df %>% 
    arrange(!!sym(rank)) %>% 
    filter(Curated_target) %>% 
    slice_head(n = 5) %>% 
    select(!!sym(rank), Symbol)
}


top5_all <- function(df) {
  list(
    Integrated = top5(df, "Rank_integrated"),
    Bind = top5(df, "Rank_binding"),
    Perturb = top5(df, "Rank_perturbation")
  )
}


top5_list <- list(
  Human = lapply(rank_list$Human, top5_all),
  Mouse = lapply(rank_list$Mouse, top5_all),
  Ortho = lapply(rank_list$Ortho, top5_all)
)



# Iteratively sample targets from curated resource and calculate AUPRC using
# a given TR's ranking (slow!). NOTE: Currently using the default sorting by 
# integrated ranking when calculating AUPRC for each sample
# ------------------------------------------------------------------------------


nreps = 1000
set.seed(20)

sample_target_list <- list(
  Human = lapply(rank_list$Human, sample_target_auprc, lt_df = lt_all, reps = nreps),
  Mouse = lapply(rank_list$Mouse, sample_target_auprc, lt_df = lt_all, reps = nreps, mouse = TRUE)
)

# Proportion of samples with greater AUPRC than observed 

prop_list <- list(
  Human = get_all_prop(auprc_list$Human, sample_target_list$Human),
  Mouse = get_all_prop(auprc_list$Mouse, sample_target_list$Mouse)
)


sample_target_list$Human$ASCL1 %>% head
prop_list$Human %>% head

# For each TR, get the AUPRC from the individual data sets composing the 
# aggregation, as well as individual rank product pairs. Focus on ordering
# perturbation by pvals
# ------------------------------------------------------------------------------



hg <- lapply(tfs_hg, function(x) {
  all_experiment_auprc(tf = x,
                       chip_meta = dat_list$Binding$Meta,
                       perturb_meta = dat_list$Perturbation$Meta,
                       bind_mat = dat_list$Binding$Human$QN_log,
                       perturb_mat = dat_list$Perturbation$Human$Pval_mat,
                       ranking_list = rank_list$Human)
})
names(hg) <- tfs_hg
  
  
  
mm <- lapply(tfs_mm, function(x) {
  all_experiment_auprc(tf = x,
                       chip_meta = dat_list$Binding$Meta,
                       perturb_meta = dat_list$Perturbation$Meta,
                       bind_mat = dat_list$Binding$Mouse$QN_log,
                       perturb_mat = dat_list$Perturbation$Mouse$Pval_mat,
                       ranking_list = rank_list$Mouse)
})
names(mm) <- tfs_mm
  
  
# Note for ortho just coerce casing of symbol to same to grab both
ortho <- lapply(tfs_hg, function(x) {
  all_experiment_auprc(
    tf = x,
    chip_meta = mutate(dat_list$Binding$Meta, Symbol = str_to_upper(Symbol)),
    perturb_meta = mutate(dat_list$Perturbation$Meta, Symbol = str_to_upper(Symbol)),
    bind_mat = dat_list$Binding$Ortho$QN_log,
    perturb_mat = dat_list$Perturbation$Ortho$Pval_mat,
    ranking_list = rank_list$Ortho)
})
names(ortho) <- tfs_hg


# Summarize the percentile AUPRC of the observed aggregated rankings relative
# to the distribution of all individual ChIP-seq and perturbation experiments,
# and their rank product pairings

perc_list <- list(
  Human = do.call(rbind, lapply(hg, `[[`, "AUPRC_percentile")),
  Mouse = do.call(rbind, lapply(mm, `[[`, "AUPRC_percentile")),
  Ortho = do.call(rbind, lapply(ortho, `[[`, "AUPRC_percentile"))
)


# Inspect example of individual experiments/rank products outperforming aggregate

filter_auprc <- function(auprc_list) {
  filter(auprc_list$AUPRC_df, AUPRC > auprc_list$AUPRC_agg$Integrated)
}


gt_list <- list(
  Human = lapply(hg, filter_auprc),
  Mouse = lapply(mm, filter_auprc),
  Ortho = lapply(ortho, filter_auprc)
)


# Most performant RUNX1 AUPRC rank product between AML ChIP-seq and K562 knockdown 
# CLO:0003679 (Human ErythroLeukemia) and CLO:0007050 (K562) most common cell types
# in low-throughput. So benchmark may be biased towards individual experiments
# of that cellular context

top_runx1 <- gt_list$Human$RUNX1 %>% arrange(desc(AUPRC))
sort(table(filter(lt_all, str_to_upper(TF_Symbol) == "RUNX1")$Cell_Type))

# Further, most performant mouse Neurod1 experiments tend to be pancreatic.
# Correspondingly, UBERON:0001264 (pancreas) is one of the top curated cell types
# (although max is cochlea... also lots of NAs)

top_neurod1 <- gt_list$Mouse$Neurod1 %>% arrange(desc(AUPRC))
sort(table(filter(lt_all, str_to_upper(TF_Symbol) == "NEUROD1")$Cell_Type))

# Topn 500 overlap of the highest performant Neurod1 pair
# ol_neurod1 <- ol_list$sim_list$Mouse$Pval_Genes[[top_neurod1$ID[1]]]
ol_neurod1 <- ol_list$sim_list$Mouse$Pval_Genes[["GSE30298_Neurod1_Mouse_Pancreatic-islets:GSE122992_Neurod1_Mouse_Overexpression"]]
targets_neurod1 <- filter(rank_list$Mouse$Neurod1, Curated_target)$Symbol
intersect(ol_neurod1, targets_neurod1)


# Rank product ordering of top pair
# TODO: shorthand for rank product

exp_chip <- data.frame(Bind = dat_list$Binding$Mouse$QN_log[, "GSE30298_Neurod1_Mouse_Pancreatic-islets"]) %>% rownames_to_column(var = "Symbol")
exp_perturb <- data.frame(Perturb = dat_list$Perturbation$Mouse$Pval_mat[, "GSE122992_Neurod1_Mouse_Overexpression"]) %>% rownames_to_column(var = "Symbol")

exp_df <- left_join(exp_chip, exp_perturb, by = "Symbol") %>% 
  mutate(Bind = rank(-Bind),
         Perturb = rank(Perturb),
         RP = rank(Bind/nrow(exp_chip) * Perturb/nrow(exp_perturb)))

filter(exp_df, Symbol %in% ol_neurod1)
filter(exp_df, Symbol %in% targets_neurod1)

# Mafa, G6pc2, Ins1, Ins2 example of curated targets with high rank product
# ranking between the top experiment pair but not in the aggregated rankings.
# Pancreatic genes

diff_df <- rank_list$Mouse$Neurod1 %>% 
  filter(Curated_target) %>% 
  left_join(exp_df, by = "Symbol")


# Plots
# ------------------------------------------------------------------------------


rank_cols <- c("Integrated" = "#1b9e77",
               "Perturbation" = "#d95f02",
               "Binding" = "#7570b3")


# Count of targets by TR

p1 <- n_target %>% 
  arrange(desc(Count)) %>% 
  mutate(Symbol = factor(Symbol, levels = unique(Symbol))) %>%
  ggplot(., aes(x = Symbol, y = Count)) +
  geom_bar(stat = "identity", fill = "slategrey", colour = "black") +
  ylab("Count of unique curated targets") +
  scale_y_continuous(breaks = seq(0, 160, 40)) +
  theme_classic() +
  theme(axis.title = element_text(size = 35),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 30),
        axis.text.x = element_text(size = 25, angle = 90, vjust = 0.5, hjust = 1))


p2 <- 
  ggplot(n_type, aes(fill = Experiment_Type, y = Count, x = TF_Symbol)) + 
  geom_bar(position = "stack", stat = "identity") +
  ylab("Count of experiments") +
  theme_classic() +
  scale_fill_manual(values = c('#66c2a5','#fc8d62','#8da0cb','#e78ac3')) +
  theme(axis.title = element_text(size = 35),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 30),
        axis.text.x = element_text(size = 25, angle = 90, vjust = 0.5, hjust = 1),
        legend.text = element_text(size = 25),
        legend.title = element_text(size = 25),
        legend.position = c(0.75, 0.85))



ggsave(plot_grid(p1, p2), dpi = 300, device = "png", height = 12, width = 24,
       filename = paste0(plot_dir, "count_curated_all.png"))


# Boxplot of ranks by presence in low-throughput


plot_box <- function(df, tf, wilx, tf_pal) {
  
  perturb <- ggplot(df, aes(x = Curated_target, y = Count_DE)) +
    geom_boxplot(width = 0.3, fill = "white", outlier.shape = NA) +
    geom_jitter(data = df[df$Curated_target,],
                aes(x = Curated_target, y = Count_DE),
                size = 3, height = 0, width = 0.1, shape = 21, alpha = 0.25,
                fill = tf_pal[tf]) +
    labs(title = tf,
         subtitle = paste0("P-value=", wilx[tf, "Perturbation"])) + 
    ylab("Count DE (FDR < 0.1)") +
    # ggtitle(tf) +
    theme_classic() +
    scale_y_continuous(breaks = pretty_breaks) +
    theme(axis.title = element_text(size = 25),
          # axis.title.x = element_blank(),
          axis.text.y = element_text(size = 25),
          axis.text.x = element_text(size = 25),
          plot.title = element_text(size = 20),
          plot.subtitle = element_text(size = 15))
  
  bind <- ggplot(df, aes(x = Curated_target, y = Mean_bind)) +
    geom_boxplot(width = 0.3, fill = "white", outlier.shape = NA) +
    geom_jitter(data = df[df$Curated_target,],
                aes(x = Curated_target, y = Mean_bind),
                size = 3, height = 0, width = 0.2, shape = 21, alpha = 0.25,
                fill = tf_pal[tf]) +
    labs(title = tf,
         subtitle = paste0("P-value=", wilx[tf, "Binding"])) + 
    ylab("Mean binding score") +
    # ggtitle(tf) +
    theme_classic() +
    theme(axis.title = element_text(size = 25),
          # axis.title.x = element_blank(),
          axis.text.y = element_text(size = 25),
          axis.text.x = element_text(size = 25),
          plot.title = element_text(size = 20),
          plot.subtitle = element_text(size = 15))
  
  cowplot::plot_grid(perturb, bind, nrow = 1)
  
}



plist_hg1 <- lapply(tfs_hg, function(x) {
  plot_box(rank_list$Human[[x]], tf = x, wilx = wilx_list$Human, tf_pal = tf_pal_hg)
}) 
names(plist_hg1) <- tfs_hg

plist_mm1 <- lapply(tfs_mm, function(x) {
  plot_box(rank_list$Mouse[[x]], tf = x, wilx = wilx_list$Mouse, tf_pal = tf_pal_mm)
}) 
names(plist_mm1) <- tfs_mm


ggsave(plot_grid(plotlist = plist_hg1, ncol = 2),
       dpi = 300, device = "png", height = 20, width = 20, bg = "white",
       filename = paste0(plot_dir, "aggregated_scores_by_curated_human_", date, ".png"))


ggsave(plot_grid(plotlist = plist_mm1, ncol = 2),
       dpi = 300, device = "png", height = 20, width = 20, bg = "white",
       filename = paste0(plot_dir, "aggregated_scores_by_curated_mouse_", date, ".png"))



# Precision vs recall curves


plot_pr <- function(pr_list, auc_df, tf, colours) {
  
  auc_labels <- c(
    Integrated = paste0("Integrated (AUC=", signif(auc_df[tf, "Integrated"], 3), ")"),
    Perturbation = paste0("Perturbation (AUC=", signif(auc_df[tf, "Perturbation"], 3), ")"),
    Binding = paste0("Binding (AUC=", signif(auc_df[tf, "Binding"], 3), ")"))
  
  ggplot() +
    geom_path(data = pr_list$Integrated, aes(x = Recall, y = Precision, col = "Integrated"), size = 1) +
    geom_path(data = pr_list$Perturbation, aes(x = Recall, y = Precision, col = "Perturbation"), size = 1) +
    geom_path(data = pr_list$Binding, aes(x = Recall, y = Precision, col = "Binding"), size = 1) +
    ggtitle(tf) +
    scale_color_manual(labels = auc_labels, values = colours) +
    theme_classic() +
    theme(axis.text = element_text(size = 25),
          axis.title = element_text(size = 25),
          plot.title = element_text(size = 25),
          legend.title = element_blank(),
          legend.text = element_text(size = 25),
          legend.position = c(0.65, 0.85))
  
}


plist_hg2 <- lapply(tfs_hg, function(x) {
  plot_pr(pr_list$Human[[x]], auc_df = auprc_list$Human, tf = x, colours = rank_cols)
})
names(plist_hg2) <- tfs_hg


plist_mm2 <- lapply(tfs_mm, function(x) {
  plot_pr(pr_list$Mouse[[x]], auc_df = auprc_list$Mouse, tf = x, colours = rank_cols)
})
names(plist_mm2) <- tfs_mm


ggsave(plist_hg2$ASCL1, dpi = 300, device = "png", height = 8, width = 8,
       filename = paste0(plot_dir, "Human_ASCL1_precisionrecall.png"))


ggsave(plist_mm2$Pax6, dpi = 300, device = "png", height = 8, width = 8,
       filename = paste0(plot_dir, "Mouse_pax6_precisionrecall.png"))


# Table of AUPRC values


auprc_table <- function(auprc_df, outfile) {
  
  pheatmap(auprc_df,
           cluster_cols = FALSE,
           cluster_rows = FALSE,
           color = "white",
           display_numbers = TRUE,
           number_format = "%.3f",
           number_color = "black",
           fontsize = 20,
           legend = FALSE,
           cellwidth = 50,
           cellheight = 50,
           angle_col = 90,
           labels_col = c("Integrated", "Binding", "Perturbation"),
           width = 4.5,
           height = 7.5,
           filename = outfile)
}


auprc_table(auprc_list$Human, paste0(plot_dir, "Human_AUPRC_table.png"))
auprc_table(auprc_list$Mouse, paste0(plot_dir, "Mouse_AUPRC_table.png"))
auprc_table(auprc_list$Ortho, paste0(plot_dir, "Ortho_AUPRC_table.png"))


# Distn of sampled AUPRCs overlaid with observed


plot_sample_auprc <- function(sample_list, auprc_df, tf) {
  
  data.frame(AUPRC = sample_list[[tf]]) %>% 
    ggplot(., aes(x = AUPRC)) +
    geom_histogram(bins = 100) +
    geom_vline(aes(xintercept = auprc_df[tf, "Integrated"], col = "Integrated"), size = 1.1, linetype = "solid") +
    geom_vline(aes(xintercept = auprc_df[tf, "Perturbation"], col = "Perturbation"), size = 1.1, linetype = "solid") +
    geom_vline(aes(xintercept = auprc_df[tf, "Binding"], col = "Binding"), size = 1.1, linetype = "solid") +
    xlim(c(NA, max(auprc_df[tf,])*1.2)) +  # pad xlim as observed typically greater
    ylab("Count") +
    ggtitle(paste0(tf, ": Sampled targets relative to curated")) +
    scale_colour_manual(values = rank_cols, labels = c("Integrated", "Binding", "Perturbation")) +
    theme_classic() +
    theme(axis.text = element_text(size = 25),
          axis.title = element_text(size = 25),
          plot.title = element_text(size = 25),
          legend.title = element_blank(),
          legend.text = element_text(size = 25),
          legend.position = c(0.90, 0.90),
          plot.margin = margin(10, 20, 10, 10))
}



plist_hg3 <- lapply(tfs_hg, function(x) {
  plot_sample_auprc(sample_target_list$Human, auprc_df = auprc_list$Human, x)
})
names(plist_hg3) <- tfs_hg


plist_mm3 <- lapply(tfs_mm, function(x) {
  plot_sample_auprc(sample_target_list$Mouse, auprc_df = auprc_list$Mouse, x)
})
names(plist_mm3) <- tfs_mm


ggsave(plist_hg3$ASCL1, dpi = 300, device = "png", height = 8, width = 12,
       filename = paste0(plot_dir, "Human_ASCL1_sample_AUPRC.png"))


# Table of proportion of sampled target AUPRCs that were greater than observed aggregated


proportion_table <- function(prop_df, outfile) {
  
  pheatmap(prop_df,
           cluster_cols = FALSE,
           cluster_rows = FALSE,
           color = "white",
           display_numbers = TRUE,
           number_format = "%.3f",
           number_color = "black",
           fontsize = 20,
           legend = FALSE,
           cellwidth = 50,
           cellheight = 50,
           angle_col = 90,
           labels_col = c("Integrated", "Binding", "Perturbation"),
           width = 4.5,
           height = 7.5,
           filename = outfile)
}


proportion_table(prop_list$Human, paste0(plot_dir, "Human_sampled_AUPRC_gt_observed_table.png"))
proportion_table(prop_list$Mouse, paste0(plot_dir, "Mouse_sampled_AUPRC_gt_observed_table.png"))


# Distribution of individual experiment AUPRCs overlaid with observed aggregated


plot_group_auprc <- function(auprc_list, tf) {
  
  ggplot(auprc_list$AUPRC_df, aes(x = Group, y = AUPRC)) +
    geom_violin(width = 0.4, fill = "lightslategrey") +
    geom_boxplot(width = 0.1, fill = "white") +
    geom_hline(aes(yintercept = auprc_list$AUPRC_agg$Integrated, col = "Integrated"), size = 1, linetype = "solid") +
    geom_hline(aes(yintercept = auprc_list$AUPRC_agg$Perturbation, col = "Perturbation"), size = 1, linetype = "solid") +
    geom_hline(aes(yintercept = auprc_list$AUPRC_agg$Binding, col = "Binding"), size = 1, linetype = "solid") +
    scale_color_manual("Integrated: ", values = rank_cols, labels = c("Integrated", "Binding", "Perturbation")) +
    scale_x_discrete("Individual experiments", labels = c("Binding", "Perturbation", "Rank product")) +
    ggtitle(paste0(tf, ": Aggregated relative to individual experiments")) +
    theme_classic() +
    theme(axis.text = element_text(size = 25),
          axis.title = element_text(size = 25),
          plot.title = element_text(size = 25),
          legend.title = element_text(size = 25),
          legend.text = element_text(size = 25),
          plot.margin = margin(10, 20, 10, 10))
}


plist_hg4 <- lapply(tfs_hg, function(x) plot_group_auprc(hg[[x]], x))
names(plist_hg4) <- tfs_hg

plist_mm4 <- lapply(tfs_mm, function(x) plot_group_auprc(mm[[x]], x))
names(plist_mm4) <- tfs_mm

plist_ortho4 <- lapply(tfs_hg, function(x) plot_group_auprc(ortho[[x]], x))
names(plist_ortho4) <- tfs_hg


# save out ASCL1 with no legend for publication
ggsave(plist_hg4$ASCL1 + theme(legend.position = "none"), 
       dpi = 300, device = "png", height = 8, width = 12,
       filename = paste0(plot_dir, "Human_ASCL1_single_experiment_AUPRC.png"))


# Table of percentile of aggregated AUPRCS relative to distn of individual experiments


percentile_table <- function(perc_df, outfile) {
  
  pheatmap(perc_df,
           cluster_cols = FALSE,
           cluster_rows = FALSE,
           color = "white",
           display_numbers = TRUE,
           number_format = "%.3f",
           number_color = "black",
           fontsize = 20,
           legend = FALSE,
           cellwidth = 50,
           cellheight = 50,
           angle_col = 90,
           labels_col = c("Integrated", "Binding", "Perturbation"),
           width = 4.5,
           height = 7.5,
           filename = outfile)
}


percentile_table(perc_list$Human, paste0(plot_dir, "Human_aggregate_AUPRC_percentile_table.png"))
percentile_table(perc_list$Mouse, paste0(plot_dir, "Mouse_aggregate_AUPRC_percentile_table.png"))
percentile_table(perc_list$Ortho, paste0(plot_dir, "Ortho_aggregate_AUPRC_percentile_table.png"))
