## Functions to create and work with matrices and dataframes containing 
## similarity metrics
## TODO: re-consider overlap matrices: parallel impl as well as whether all 4
## matrices are required
## -----------------------------------------------------------------------------


library(parallel)
library(Matrix)
library(GeneOverlap)
library(plyr)
library(dplyr)
library(WGCNA)
library(ggplot2)
library(cowplot)


mat_to_df <- function(mat, symmetric = TRUE) {
  
  # Assumes m x n mat with named rows and columns. Returns a dataframe of all
  # the row-col elements - if symmetric is TRUE, then only keep unique pairs
  # https://stackoverflow.com/questions/28035001/
  
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



filter_na_cols <- function(mat) {
  # Remove columns that are NA in every entry
  # https://stackoverflow.com/questions/2643939/
  as.matrix(Filter(function(x) !all(is.na(x)), as.data.frame(mat)))
}


jaccard <- function(matrix) {
  # Takes in a named matrix and returns the jaccard coeff for each pairwise 
  # comparison as a matrix
  # https://stats.stackexchange.com/a/89947
  
  # common values:
  A <- Matrix::tcrossprod(matrix)
  # indexes for non-zero common values
  im <-  which(A > 0, arr.ind=TRUE)
  # counts for each row
  b <-  rowSums(matrix)
  # only non-zero values of common
  Aim <-  A[im]
  # Jacard formula: #common / (#i + #j - #common)
  jmat <-  as.matrix(Matrix::sparseMatrix(
    i = im[,1],
    j = im[,2],
    x = Aim / (b[im[,1]] + b[im[,2]] - Aim),
    dims = dim(A)
  ))
  
  rownames(jmat) <- colnames(jmat) <- rownames(matrix)
  
  return(jmat)
}



format_pair_df <- function(mat, symmetric = TRUE, perturb = FALSE) {
  
  # create data frame of unique comparisons in a named mat and then add columns
  # for the compared TFs/species and whether they are the same. If perturb is
  # TRUE, will also add a column for whether the perurbations match
  # Assumes mat has rows and cols named as: [GSE]_[TF]_[Species]_[Context|Perturbation]
  
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
      return ("In_TR_in_species")
    } else if (tf1[i] == tf2[i] & species1[i] != species2[i]) {
      return ("In_TR_out_species")
    } else {
      return ("Out")
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


get_topn_genes <- function(vec, n, decrease = TRUE) {
  # Helper to extract top n genes from a named vector of values
  stopifnot(is.numeric(vec) | length(names(vec)) == 0)
  names(head(sort(vec, decreasing = decrease), n))
}


get_common <- function(vec1, vec2) {
  # Assumes vec1 and vec2 are named - return intersect names
  common <- intersect(names(vec1[!is.na(vec1)]), names(vec2[!is.na(vec2)]))
}


get_overlap_matrices <- function(mat,  # gene x exp matrix
                                 topn,  # how many top symbols to keep
                                 decrease_arg = TRUE, # decrease sort?
                                 common_only = FALSE) { # only consider shared symbols
  
  
  # performs gene overlap for every pairwise column of a matrix, keeping only
  # the top n gene symbols for each experiment. fills out 4 matrices of these
  # comparisons: jaccard coefficients, FET pvals, odds ratio, and the length
  # of the intersect
  
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


chip_sim_list <- function(mat, topn) {
  
  # Mat as a gene x exp matrix - returns a list of 5 matrices of similarity
  # metrics between the cols/exps of mat: Pcor, length of intersect of topn, 
  # Jaccard of topn, FET pval of topn, and log odds ratio of topn
  
  overlap <- get_overlap_matrices(mat, topn)
  pcor <- WGCNA::cor(mat)
  sim_list <- append(list(Pcor = pcor), overlap)
  return(sim_list)
}


perturb_sim_list <- function(fc_mat, 
                             pval_mat, 
                             topn, 
                             common_arg = TRUE) {
  
  # mat_list must contain gene x exp pval and FC perturb matrices - generates 
  # list of matrices containing exp by exp similarity metrics: Pcor +/- abs FC, 
  # and for each of topn pval/upreg/downreg: length of intersect, Jaccard, FET 
  # pval, and log odds ratio of topn
  
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


count_common <- function(df, mat, ncores = 8) {
  # Helper that takes rows of a paired df and uses mat to find out
  # how many genes were measured in common between the correlated experiments
  
  unlist(mclapply(1:nrow(df), function(x) {
    exp1 <- df$Row[x]
    exp2 <- df$Col[x]
    length(intersect(names(na.omit(mat[, exp1])), names(na.omit(mat[, exp2]))))
  }, mc.cores = ncores))
}


format_and_merge <- function(mat_list, 
                             gene_mat = NULL, 
                             add_common = FALSE,
                             symmetric_arg = TRUE) {
  
  # Take in a list of named similarity matrices, format each into a df of the
  # unique pairs, rename generic stat column to the current stat, then merge all
  # into a single df. If add_common is TRUE, add column of how many genes are
  # commonly measured between the paired experiments. Symmetric refers to the
  # original gene x experiment matrices used to generate similarity between 
  # exps/cols: perturb/chip produce symmetric matrices (so only want unique
  # elements) while intersect produces non-symmetric (want everything)
  
  stats <- names(mat_list)
  stats <- stats[!str_detect(stats, "^.*_Genes")]  # Hacky removal of intersect list of genes (not sim mat)
  
  
  sim_list <- lapply(stats, function(x) {
    
    df <- format_pair_df(mat_list[[x]], symmetric = symmetric_arg)
    
    df <- mutate(df, Key = paste(Row, Col, sep = ":")) %>% 
      relocate(Value, .after = last_col())
    
    colnames(df)[colnames(df) == "Value"] <- x
    
    return(df)
  })
  
  join_names <- Reduce(intersect, lapply(sim_list, colnames))
  df <- plyr::join_all(sim_list, by = join_names) %>% select(-Key)
  
  if (add_common) {
    df$Common_genes <- count_common(df, gene_mat)
  }

  return(df)
}


overlap_chip_perturb <- function(chip_mat,  # gene x exp matrix,
                                 perturb_mat,  # gene x exp matrix
                                 topn,  # how many top symbols to keep
                                 decrease_perturb = TRUE, # decrease sort for perturb?
                                 common_only = FALSE) { # only consider shared symbols
  
  
  # performs gene overlap for every pairwise columns/experiments of chip binding
  # matrix and a peturb effect size matrix, keeping only the top n gene symbols 
  # for each experiment. fills out 4 matrices of these comparisons: jaccard, 
  # FET pvals, odds ratio, and the length of the intersect
  
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


intersect_sim_list <- function(chip_mat, 
                               fc_mat, 
                               pval_mat, 
                               topn, 
                               common_arg = TRUE) {
  
  # Perform overlap_chip_perturb() for each of sorting perturbation data by
  # pvals, positive FC, and negative FC
  
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


get_summary <- function(df) {
  
  # Return a dataframe of the min/max/mean/median of numeric cols by Group
  
  df %>% 
    group_by(Group) %>% 
    select_if(is.numeric) %>%
    dplyr::summarise(across(where(is.numeric), list(
      Min = ~min(.x, na.rm = TRUE),
      Max = ~max(.x, na.rm = TRUE),
      Mean = ~mean(.x, na.rm = TRUE),
      Median = ~median(.x, na.rm = TRUE)
    )))
}


tf_summary <- function(df) {
  
  # get_summary for every unique TF in df
  
  tfs <- unique(c(as.character(df$TF1), as.character(df$TF2)))
  tf_list <- lapply(tfs, function(x) {
    tf_df <- filter(df, as.character(TF1) == x | as.character(TF2) == x)
    get_summary(tf_df)
  })
  names(tf_list) <- tfs
  return(tf_list)
}



# The following functions were originally used for overlapping chip and perturb
# matrix. They separated the retrieval of the geneoverlap object and the 
# construction of the stat matrices. Use of mapply however should be considered
# over current usage of for loops

# top_overlap <- function(chip_mat, 
#                         perturb_mat, 
#                         topn,
#                         decrease_sort = FALSE,
#                         common_only = TRUE,
#                         cores = 8) {
#   
#   # chip_mat and perturb_mat are gene x experiment matrices of equal nrows (genes).
#   # uses GeneOverlap to find the overlap of the 'topn' scoring genes between
#   # each column of chip and perturb mat. returns a list of size 
#   # ncol(chip_mat) * ncol(row_mat) containing the size of the intersect, the
#   # jaccard, the Fisher's Exact Test pval, the log odds ratio, and the vector
#   # of overlapping genes
#   
#   stopifnot(identical(nrow(chip_mat), nrow(perturb_mat)))
#   all_common <- intersect(rownames(chip_mat), rownames(perturb_mat))
#   
#   # indices take form (chip_mat[, i], perturb_mat[, j]): (1,1), (1,2), (1,3) ...
#   i_vec <- rep(seq_len(ncol(chip_mat)), each = ncol(perturb_mat))
#   j_vec <- rep(seq_len(ncol(perturb_mat)), times = ncol(chip_mat))
#   intersect_names <- paste0(colnames(chip_mat)[i_vec], ":", colnames(perturb_mat)[j_vec])
#   
#   overlap_list <- mcmapply(i_vec, j_vec, FUN = function(i, j) {
#     
#     exp1 <- chip_mat[, i]
#     exp2 <- perturb_mat[, j]
#     exp2 <- exp2[!is.na(exp2)]
#     common <- intersect(names(exp1), names(exp2))
# 
#     top1 <- get_topn_genes(exp1, topn)
#     top2 <- get_topn_genes(exp2, topn, decrease = decrease_sort)
#     
#     # generate GeneOverlap object
#     gene_overlap <-
#       testGeneOverlap(newGeneOverlap(
#         listA = top1,
#         listB = top2,
#         genome.size = length(all_common)
#       ))
#     
#     # extract relevant statistics/gene identities into a list
#     list(Intersect = length(gene_overlap@intersection),
#          Jaccard = gene_overlap@Jaccard,
#          FET_pval = gene_overlap@pval,
#          Log_OR = log(gene_overlap@odds.ratio),
#          Total_common_genes = length(common),
#          Genes = gene_overlap@intersection)
#     
#   }, mc.cores = cores, SIMPLIFY = FALSE)
#   
#   names(overlap_list) <- intersect_names
#   return(overlap_list)
# }
# 
# 
# overlap_matrix <- function(overlap_list, chip_mat, perturb_mat, overlap_type) {
#   # Generate a matrix of dim ncol(chip_mat) x ncol(perturb_mat) that contains
#   # the overlap statistics for each experiment of chip_mat and perturb_mat
#   
#   stopifnot(overlap_type %in% c("Intersect", "Jaccard", "FET_pval", "Log_OR"))
#   
#   mat <- matrix(unlist(lapply(overlap_list, `[`, overlap_type)),
#                 nrow = ncol(chip_mat),
#                 ncol = ncol(perturb_mat),
#                 byrow = TRUE)  # list output is indexed (1,1), (1,2), (1,3) ...
#   
#   rownames(mat) <- colnames(chip_mat)
#   colnames(mat) <- colnames(perturb_mat)
#   return(mat)
#   
# }



# The following functions were an attempt to implement a version of generating
# a cor table by first generating all choose(n, 2) pairs of column wise vectors
# to correlate, and then performing these correlations via parallel. However,
# it ended up being slower than the normal stats::cor method (and much slower
# than WGCNA::cor). 


# get_colwise_pairs <- function(mat) {
#   # Helper to get a list of all the pairwise combination of column in mat
#   combn(1:ncol(mat), 2, simplify = FALSE)
# }


# cor_colwise_pairs <- function(pair, mat, ...) {
#   # Helper to get the correlation of a single pair from mat. Assumes pair
#   # contains the column indices to correlate
#   stopifnot(length(pair) == 2 & is.integer(pair))
#   cor(mat[, pair[1]], mat[, pair[2]], ...)
# }


# mc_cor_colwise <- function(mat_arg, 
#                            cores = 8, 
#                            use_arg = "pairwise.complete.obs",
#                            method_arg = "pearson") {
#   
#   # https://davetang.org/muse/2012/01/31/creating-a-correlation-matrix-with-r/
#   # Given an n x m gene matrix, return a choose(m, 2) x 3 dataframe of all unique
#   # pairwise combinations of correlations in the column of the matrix
#   
#   pairs <- get_colwise_pairs(mat)
#   
#   cor_list <- mclapply(pairs, function(x) {
#     cor_colwise_pairs(
#       pair = x,
#       mat = mat_arg,
#       use = use_arg,
#       method = method_arg
#     )
#   }, mc.cores = cores)
#   
#   # Format into dataframe, matching gene names from original matrix
#   cor_df <- as.data.frame(cbind(do.call(rbind, pairs), unlist(cor_list)))
#   names(cor_df) <- c("Gene1", "Gene2", "Cor")
#   cor_df$Gene1 <- colnames(mat)[cor_df$Gene1]
#   cor_df$Gene2 <- colnames(mat)[cor_df$Gene2]
#   cor_df$Cor <- round(cor_df$Cor, 5)
#   
#   return (cor_df)
#   
# }
