## Organizing plot functions used throughout the paper
## -----------------------------------------------------------------------------

library(tidyverse)
library(pheatmap)
library(cowplot)
library(RColorBrewer)
library(ggsci)
library(viridisLite)


# Colour definitions
# ------------------------------------------------------------------------------

# https://colorswall.com/palette/73/

tf_pal = c(
  Ascl1 = "#00188f",
  Hes1 = "#fff100",
  Mecp2 = "#009e49",
  Mef2c = "#00bcf2",
  Neurod1 = "#68217a",
  Pax6 = "#ffd8b1",
  Runx1 = "#e81123",
  Tcf4 = "#ff8c00"
)


pert_pal <- ggsci::pal_npg("nrc")(4)


pert_anno <- list(
  Perturbation = c(
    "Overexpression" = pert_pal[1],
    "Knockdown" = pert_pal[2],
    "Knockout" = pert_pal[3],
    "Mutant" = pert_pal[4]
  )
)


pal_length <- 100
bluered_pal <- colorRampPalette(c("#0571b0", "white", "#ca0020"))(pal_length)
depr_pal <- rev(magma(n = 11))



# General functions
# ------------------------------------------------------------------------------


# For better handling of integer axis breaks

pretty_breaks <- function(x) {
  # https://stackoverflow.com/questions/15622001/how-to-display-only-integer-values-on-an-axis-using-ggplot2
  unique(floor(pretty(seq(0, (max(x) + 1) * 1.1))))
}



# Experiment similarity functions
# ------------------------------------------------------------------------------


# Violin+boxplot of similarity stats by group

stat_vboxplot <- function(df, 
                          x_var = "Group", 
                          y_var, 
                          y_name, 
                          title,
                          ortho = FALSE) {
  
  stopifnot(identical(levels(df$Group), 
                      c("In_TR_in_species", "In_TR_out_species", "Out")))
  
  if (ortho) {
    levels(df$Group) <- c("TR+ S+", "TR+ S-", "Out")
  } else {
    levels(df$Group) <- c("Intra-TR", "TR+ S-", "Inter-TR")
  }
  
  ggplot(df, aes(x = !!sym(x_var), y = !!sym(y_var))) +
    geom_violin(width = 0.4, fill = "lightslategrey") +
    geom_boxplot(width = 0.1, fill = "white") +
    ylab(y_name) +
    ggtitle(title) +
    ylab(y_name) +
    theme_classic() +
    theme(axis.title = element_text(size = 25),
          axis.title.x = element_blank(),
          axis.text.y = element_text(size = 25),
          axis.text.x = element_text(size = 25),
          plot.title = element_text(size = 25, hjust = 0.5))
}


# stat_vboxplot for every unique TR


tf_vboxplot <- function(df, 
                        x_var = "Group", 
                        y_var, 
                        y_name, 
                        title_prefix,
                        ortho = FALSE) {
  
  tfs <- intersect(df$TF1, df$TF2)
  
  # list of plots by filtering df for tfs in plot_tfs
  
  tf_list <- lapply(tfs, function(x) {
    
    tf_df <- filter(df, TF1 == x | TF2 == x)
    
    stat_vboxplot(tf_df,
                  y_var = y_var,
                  y_name = y_name,
                  title = paste(title_prefix, x, sep = " "),
                  ortho = ortho)
  })
  names(tf_list) <- tfs
  return(tf_list)
}


# Density plots of stat by group


dplot <- function(df,
                  stat,
                  stat_name,
                  species) {
  
  if(species == "Human") {
    fill_col <- c("royalblue", "grey30")
    legend_text = c("Intra-TR", "Inter-TR")
  } else if (species == "Mouse") {
    fill_col <- c("goldenrod", "grey30")
    legend_text = c("Intra-TR", "Inter-TR")
  } else if (species == "Ortho") {
    fill_col <- c("forestgreen", "mediumvioletred", "lightgrey")
    legend_text = c("Intra-TR & within-species", "Intra-TR & cross-species", "Inter-TR (both species)")
  }
  
  
  ggplot(df, aes(x = !!sym(stat), fill = Group)) +
    geom_density(alpha = 0.6) +
    theme_classic() +
    ylab("Density") +
    xlab(stat_name) +
    ggtitle(species) +
    scale_fill_manual(values = fill_col, labels = legend_text) +
    theme(
      axis.text = element_text(size = 25),
      axis.title = element_text(size = 25),
      plot.title = element_text(hjust = 0.5, size = 30),
      legend.position = c(0.75, 0.85),
      legend.text = element_text(size = 20),
      legend.title = element_blank(),
      plot.margin = margin(10, 20, 10, 10)  # xaxis was getting clipped for pcor
    )
}


# ChIP-seq functions
# ------------------------------------------------------------------------------





# Perturbation functions
# ------------------------------------------------------------------------------


# Filters the FC matrix and meta for the relevant TF/symbols. Clips outlier FC
# values to min/max

fc_heatmap <- function(tf,
                       tf_df,
                       fc_mat,
                       meta,
                       FC_min = -2.5,
                       FC_max = 2.5,
                       legend_arg = TRUE,
                       anno_legend = TRUE) {
  
  
  symbols <- tf_df$Symbol
  tf_meta <- filter(meta, Symbol == tf) %>% arrange(Perturbation)
  fc_mat <- fc_mat[symbols, tf_meta$Experiment_ID]
  
  anno <- tf_meta[, "Perturbation", drop = FALSE]
  anno$Perturbation <- factor(anno$Perturbation, levels = unique(anno$Perturbation))
  droplevels(anno)
  rownames(anno) <- colnames(fc_mat)
  

  color_breaks <- seq(FC_min, FC_max, length.out = pal_length)
  
  pheatmap(
    fc_mat,
    color = bluered_pal,
    breaks = color_breaks,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    show_colnames = FALSE,
    border_col = "black",
    annotation_col = anno,
    annotation_colors = pert_anno,
    annotation_names_col = FALSE,
    annotation_legend = anno_legend,
    cellheight = 20,
    cellwidth = 20,
    fontsize_row = 20,
    legend = legend_arg,
    height = 6,
    width = 9
  )
}


# Heatmap for DE prior - if binary is TRUE, will use cutoff to create black/white
# cells instead of continuous 


depr_heatmap <- function(tf_df,
                         deprior,
                         deprior_binary = FALSE,
                         binary_cutoff = 0.9,
                         legend_arg = TRUE) {

  symbols <- tf_df$Symbol

  deprior <- deprior %>%
    filter(Symbol %in% symbols) %>%
    arrange(match(Symbol, symbols)) %>%
    dplyr::select(DE_Prior_Rank)
  
  if (deprior_binary) {
    
    deprior$DE_Prior_Rank <- as.numeric(deprior$DE_Prior_Rank > binary_cutoff)
    pal <- c("white", "black")
    
    # incredibly hacky fix to deal with the fact that pheatmap can't handle
    # when all values are the same (why...)
    if (all(deprior$DE_Prior_Rank == 0)) {
      pal <- "white"
      deprior$DE_Prior_Rank[1] <- 1
    } else if (all(deprior$DE_Prior_Rank == 1)) {
      pal <- "black"
      deprior$DE_Prior_Rank[1] <- 0
    }
    
    legend_arg <- FALSE
    
  } else {
    pal <- depr_pal
  }
  
  pheatmap(
      deprior[, "DE_Prior_Rank", drop = FALSE],
      color = pal,
      cluster_rows = FALSE,
      cluster_cols = FALSE,
      show_colnames = FALSE,
      show_rownames = FALSE,
      border_col = "black",
      cellheight = 20,
      cellwidth = 20,
      legend = legend_arg,
      height = 6,
      width = 9
    )
}

