## Functions to create and work with matrices and dataframes containing 
## similarity metrics
## -----------------------------------------------------------------------------

library(parallel)
library(Matrix)
library(GeneOverlap)
library(plyr)
library(dplyr)
library(WGCNA)


# Assumes m x n mat with named rows and columns. Returns a dataframe of all
# the row-col elements - if symmetric is TRUE, then only keep unique pairs
# https://stackoverflow.com/questions/28035001/

mat_to_df <- function(mat, symmetric = TRUE) {
  
  if (symmetric) {
    df <- data.frame(
      Row = rownames(mat)[row(mat)[lower.tri(mat)]],
      Col = colnames(mat)[col(mat)[lower.tri(mat)]],
      Value = mat[lower.tri(mat)],
      stringsAsFactors = FALSE
    )
  } else {
    df <- data.frame(
      Row = rownames(mat)[row(mat)],
      Col = colnames(mat)[col(mat)],
      Value = c(mat),
      stringsAsFactors = FALSE
    )
  }
  return(df)
}


# Remove columns that are NA in every entry
# https://stackoverflow.com/questions/2643939/

filter_na_cols <- function(mat) {

  as.matrix(Filter(function(x) !all(is.na(x)), as.data.frame(mat)))
}


# Create a data frame of elements in a named mat and add columns for the 
# compared TFs/species and whether they are the same. If perturb is
# TRUE, will also add a column for whether the perturbations are the same.
# Symmetric matrix only unique elements are kept.
# Assumes mat has rows and cols named as: [GSE]_[TF]_[Species]_[Context|Perturbation]

format_pair_df <- function(mat, symmetric = TRUE, perturb = FALSE) {
  
  df <- mat_to_df(mat, symmetric = symmetric)
  
  row_split <- str_split(df$Row, "_", simplify = TRUE)
  col_split <- str_split(df$Col, "_", simplify = TRUE)
  
  geo1 <- row_split[, 1]
  geo2 <- col_split[, 1]
  
  tf1 <- str_to_title(row_split[, 2])
  tf2 <- str_to_title(col_split[, 2])
  
  species1 <- row_split[, 3]
  species2 <- col_split[, 3]
  
  context1 <- row_split[, 4]
  context2 <- col_split[, 4]
  
  group <- vapply(1:nrow(row_split), function(i) {
    if (tf1[i] == tf2[i] & species1[i] == species2[i]) {
      return("In_TR_in_species")
    } else if (tf1[i] == tf2[i] & species1[i] != species2[i]) {
      return("In_TR_out_species")
    } else {
      return("Out")
    }
  }, FUN.VALUE = character(1))
  
  group <- factor(group,
                  levels = c("In_TR_in_species", "In_TR_out_species", "Out"))
  
  pair_df <- data.frame(df, 
                        GEO1 = geo1,
                        GEO2 = geo2,
                        TF1 = tf1, 
                        TF2 = tf2, 
                        Group = group)
  
  if (perturb) pair_df$Same_perturb <- context1 == context2
  
  return(pair_df)
}


# Split a pair into a list for the intra-TR experiments for each TR

split_pair_df <- function(pair_df) {
  
  tfs <- unique(c(as.character(pair_df$TF1), as.character(pair_df$TF2)))
  
  tf_list <- lapply(tfs, function(x) {
    filter(pair_df, as.character(TF1) == x & as.character(TF2) == x)
  })
  names(tf_list) <- tfs
  
  return(tf_list)
}


# Helper to extract top n genes from a named vector of values

get_topn_genes <- function(vec, n, decrease = TRUE) {
  
  stopifnot(is.numeric(vec) | length(names(vec)) == 0)
  names(head(sort(vec, decreasing = decrease), n))
}


# Assumes vec1 and vec2 are named - return intersect names

get_common <- function(vec1, vec2) {
  
  intersect(names(vec1[!is.na(vec1)]), names(vec2[!is.na(vec2)]))
}


# Helper that takes rows of a paired df and uses mat to find out how many genes
# were measured in common between the two experiments

count_common <- function(df, mat, ncores) {
  
  l <- mclapply(1:nrow(df), function(x) {
    exp1 <- df$Row[x]
    exp2 <- df$Col[x]
    length(intersect(names(na.omit(mat[, exp1])), names(na.omit(mat[, exp2]))))
  }, mc.cores = ncores)
  
  return(unlist(l))
}


# Performs gene overlap for every pairwise column of a matrix, keeping only
# the top n gene symbols for each experiment. Used on ChIP-seq on Perturbation
# effect sizes matrices. Fills out 4 matrices of overlap metrics: jaccard 
# coefficients, FET pvals, odds ratio, and intersect length.

get_overlap_matrices <- function(mat,  # gene x exp matrix
                                 topn,  # how many top symbols to keep
                                 decrease_arg = TRUE, # decrease sort?
                                 common_only = FALSE) { # only consider shared symbols
  
  # init output matrices
  jac_mat <- matrix(NA, nrow = ncol(mat), ncol = ncol(mat))
  rownames(jac_mat) <- colnames(jac_mat) <- colnames(mat)
  fetpval_mat <- or_mat <- intersect_mat <- jac_mat
  
  for (i in 1:ncol(mat)) {
    for (j in 1:ncol(mat)) {
      if (i == j) {
        next
      }
      
      exp1 <- mat[, i]
      exp2 <- mat[, j]
      
      if (common_only) {
        exp1 <- exp1[names(exp1) %in% get_common(exp1, exp2)]
        exp2 <- exp2[names(exp2) %in% get_common(exp1, exp2)]
      }
      
      a <- get_topn_genes(exp1, topn, decrease = decrease_arg)
      b <- get_topn_genes(exp2, topn, decrease = decrease_arg)
      
      gene_overlap <- testGeneOverlap(newGeneOverlap(  # overlap object
        listA = a,
        listB = b,
        genome.size = nrow(mat)))
      
      # fill matrices with overlap stats
      jac_mat[i, j] <- gene_overlap@Jaccard
      fetpval_mat[i, j] <- gene_overlap@pval
      intersect_mat[i, j] <- length(gene_overlap@intersection)
      or_mat[i, j] <- log(gene_overlap@odds.ratio)
    }
  }
  
  return(list(Intersect = intersect_mat, 
              Jaccard = jac_mat, 
              FET_pval = fetpval_mat,
              Log_OR = or_mat))
}


# Performs gene overlap for every pairwise column of a matrix, keeping only
# the top n gene symbols for each experiment. Used *BETWEEN* columns of ChIP-seq
# and Perturbation effect sizes matrices. Fills out 4 matrices of overlap 
# metrics: jaccard coefficients, FET pvals, odds ratio, intersect length. Along
# with these matrices, also export a character vector of the overlapping genes
# between each comparison.

overlap_chip_perturb <- function(chip_mat,  # gene x exp matrix,
                                 perturb_mat,  # gene x exp matrix
                                 topn,  # how many top symbols to keep
                                 decrease_perturb = TRUE, # decrease sort for perturb?
                                 common_only = FALSE) { # only consider shared symbols
  
  # init output matrices
  jac_mat <- matrix(NA, nrow = ncol(chip_mat), ncol = ncol(perturb_mat))
  rownames(jac_mat) <- colnames(chip_mat)
  colnames(jac_mat) <- colnames(perturb_mat)
  fetpval_mat <- or_mat <- intersect_mat <- jac_mat
  
  # init list of overlapping genes
  overlap_genes <- 
    vector(mode = "list", length = ncol(chip_mat) * ncol(perturb_mat))
  
  i_vec <- rep(seq_len(ncol(chip_mat)), each = ncol(perturb_mat))
  j_vec <- rep(seq_len(ncol(perturb_mat)), times = ncol(chip_mat))
  
  names(overlap_genes) <- 
    paste0(colnames(chip_mat)[i_vec], ":", colnames(perturb_mat)[j_vec])
  
  
  for (i in 1:ncol(chip_mat)) {
    for (j in 1:ncol(perturb_mat)) {
      
      chip_exp <- chip_mat[, i]
      perturb_exp <- perturb_mat[, j]
      
      if (common_only) {
        common <- get_common(chip_exp, perturb_exp)
        chip_exp <- chip_exp[common]
        perturb_exp <- perturb_exp[common]
      }
      
      top_chip <- get_topn_genes(chip_exp, topn, decrease = TRUE)
      top_perturb <- get_topn_genes(perturb_exp, topn, decrease = decrease_perturb)
      
      gene_overlap <- testGeneOverlap(newGeneOverlap(  # overlap object
        listA = top_chip,
        listB = top_perturb,
        genome.size = nrow(chip_mat)))
      
      # fill matrices with overlap stats
      jac_mat[i, j] <- gene_overlap@Jaccard
      fetpval_mat[i, j] <- gene_overlap@pval
      intersect_mat[i, j] <- length(gene_overlap@intersection)
      or_mat[i, j] <- log(gene_overlap@odds.ratio)
      
      # fill list of genes
      overlap_name <- paste0(colnames(chip_mat)[i], ":", colnames(perturb_mat)[j])
      overlap_genes[[overlap_name]] <- gene_overlap@intersection
    }
  }
  
  return(list(Intersect = intersect_mat, 
              Jaccard = jac_mat, 
              FET_pval = fetpval_mat,
              Log_OR = or_mat,
              Genes = overlap_genes))
}


# Mat as a gene x exp matrix - returns a list of 5 matrices of similarity
# metrics between the cols/exps of mat: Pcor, length of intersect of topn, 
# Jaccard of topn, FET pval of topn, and log odds ratio of topn

chip_sim_list <- function(mat, topn) {

  overlap <- get_overlap_matrices(mat, topn)
  pcor <- WGCNA::cor(mat)
  sim_list <- append(list(Pcor = pcor), overlap)
  
  return(sim_list)
}


# Generates list of similarity matrices between columns of perturbation mat.
# Pearson cor of raw and absolute fold changes. And for each of sorting by
# pvals, positive FCs and negative FCs: length of intersect of topn, Jaccard of 
# topn, FET pval of topn, and log odds ratio of topn

perturb_sim_list <- function(fc_mat, 
                             pval_mat, 
                             topn, 
                             common_arg = TRUE) {
  
  pcor <- WGCNA::cor(fc_mat, use = "pairwise.complete.obs")
  
  pcor_abs <- WGCNA::cor(abs(fc_mat), use = "pairwise.complete.obs")
  
  # Pval - sort min to max
  overlap_pval <- get_overlap_matrices(
    mat = pval_mat, 
    topn = topn, 
    decrease_arg = FALSE,
    common_only = common_arg)
  
  # Upreg - sort max to min
  overlap_up <- get_overlap_matrices(
    mat = fc_mat, 
    topn = topn, 
    decrease_arg = TRUE,
    common_only = common_arg)
  
  # Downreg - sort min to max
  overlap_down <- get_overlap_matrices(
    mat = fc_mat, 
    topn = topn, 
    decrease_arg = FALSE,
    common_only = common_arg)
 
  names(overlap_pval) <- paste("Pval", names(overlap_pval), sep = "_")
  names(overlap_up) <- paste("Upreg", names(overlap_up), sep = "_")
  names(overlap_down) <- paste("Downreg", names(overlap_down), sep = "_")
  
  sim_list <- c(list(Pcor = pcor), 
                list(Pcor_abs = pcor_abs), 
                overlap_pval, 
                overlap_up, 
                overlap_down)
  
  return(sim_list)
}


# Perform overlap_chip_perturb() for each of sorting perturbation data by
# pvals, positive FC, and negative FC

intersect_sim_list <- function(chip_mat, 
                               fc_mat, 
                               pval_mat, 
                               topn, 
                               common_arg = TRUE) {
  
  # Pval - sort min to max
  overlap_pval <- overlap_chip_perturb(
    chip_mat = chip_mat,
    perturb_mat = pval_mat,
    topn = topn,
    decrease_perturb = FALSE,
    common_only = common_arg)
  
  # Upreg - sort max to min
  overlap_up <- overlap_chip_perturb(
    chip_mat = chip_mat,
    perturb_mat = fc_mat,
    topn = topn,
    decrease_perturb = TRUE,
    common_only = common_arg)
  
  # Downreg - sort min to max
  overlap_down <- overlap_chip_perturb(
    chip_mat = chip_mat,
    perturb_mat = fc_mat,
    topn = topn,
    decrease_perturb = FALSE,
    common_only = common_arg)
  
  names(overlap_pval) <- paste("Pval", names(overlap_pval), sep = "_")
  names(overlap_up) <- paste("Upreg", names(overlap_up), sep = "_")
  names(overlap_down) <- paste("Downreg", names(overlap_down), sep = "_")
  
  sim_list <- c(overlap_pval, overlap_up, overlap_down)
  
  return(sim_list)
}


# Takes in a list of named similarity matrices, formats each into a df of the
# unique pairs, rename generic stat column to the current stat, then merge all
# into a single df. If add_common is TRUE, add column of how many genes are
# commonly measured between the paired experiments. Symmetric refers to the
# original gene x experiment matrices used to generate similarity between 
# exps/cols: perturb/chip produce symmetric matrices (so only want unique
# elements) while intersect produces non-symmetric (want everything)

format_and_merge <- function(mat_list, 
                             gene_mat = NULL, 
                             add_common = FALSE,
                             symmetric_arg = TRUE) {
  
  stats <- names(mat_list)
  
  # Hacky removal of intersect object's list of genes (which is not sim mat)
  stats <- stats[!str_detect(stats, "^.*_Genes")]
  
  # List of formatted dfs for each similarity statistic 
  sim_list <- lapply(stats, function(x) {
    
    df <- format_pair_df(mat_list[[x]], symmetric = symmetric_arg) %>%
      mutate(df, Key = paste(Row, Col, sep = ":")) %>% 
      relocate(Value, .after = last_col()) %>% 
      dplyr::rename({{x}} := Value)
    
    return(df)
  })
  
  # Get all stat names and then join all dfs into one
  join_names <- Reduce(intersect, lapply(sim_list, colnames))
  
  df <- plyr::join_all(sim_list, by = join_names) %>% select(-Key)
  
  if (add_common) {
    df$Common_genes <- count_common(df, gene_mat, ncores)
  }

  return(df)
}


# Return a dataframe of the min/max/mean/median of numeric cols by Group

get_summary <- function(df) {
  
  df %>% 
    group_by(Group) %>% 
    dplyr::summarise(across(where(is.numeric), list(
      Min = ~min(.x, na.rm = TRUE),
      Max = ~max(.x, na.rm = TRUE),
      Mean = ~mean(.x, na.rm = TRUE),
      Median = ~median(.x, na.rm = TRUE)
    )))
}


# get_summary for every unique TF in df

tf_summary <- function(df) {
  
  tfs <- unique(c(as.character(df$TF1), as.character(df$TF2)))
  
  tf_list <- lapply(tfs, function(x) {
    tf_df <- filter(df, as.character(TF1) == x | as.character(TF2) == x)
    get_summary(tf_df)
  })
  names(tf_list) <- tfs
  
  return(tf_list)
}


