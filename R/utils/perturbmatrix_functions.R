## Functions for working with perturbation data formatted as a gene x experiment
## matrix of effect sizes (FC (fold changes), FDR (false discovery rate))


count_table <- function(fdr_mat, fdr) {
  
  # Given an gene x experiment FDR mat, return a dataframe summarizing the 
  # count of times a gene was under the FDR threshold as well as the count NA
  
  deg_mat <- fdr_mat < fdr
  count_de <- rowSums(deg_mat, na.rm = TRUE)
  count_na <- apply(deg_mat, 1, function(x) sum(is.na(x)))
  count_exps <- ncol(deg_mat)
  frac_measured <- round(count_de / (count_exps - count_na), 3)
  frac_all <- round(count_de / count_exps, 3)
  
  df <- data.frame(
    Symbol = names(count_de),
    Count_DE = count_de,
    Fraction_DE_measured = frac_measured,
    Fraction_DE_all = frac_all,
    Count_NA = count_na,
    stringsAsFactors = FALSE
  )
  return (df[order(df$Count_DE, decreasing = TRUE), ])
}


tally_direction <- function(fc_mat, meta) {
  
  # Return a data frame that has for each gene the count of times it was down/up
  # in GoF (overexpression) and LoF (KO + KD) experiments
  
  gof_exps <- filter(meta, Perturbation == "Overexpression")$Experiment_ID
  lof_exps <- filter(meta, Perturbation %in% c("Knockdown", "Knockout"))$Experiment_ID
  gof_mat <- fc_mat[, gof_exps, drop = FALSE]
  lof_mat <- fc_mat[, lof_exps, drop = FALSE]
  
  gof_down <- apply(gof_mat, 1, function(x) sum(x < 0, na.rm = TRUE))
  gof_up <- apply(gof_mat, 1, function(x) sum(x > 0, na.rm = TRUE))
  lof_down <- apply(lof_mat, 1, function(x) sum(x < 0, na.rm = TRUE))
  lof_up <- apply(lof_mat, 1, function(x) sum(x > 0, na.rm = TRUE))
  
  return(data.frame(
    Symbol = rownames(fc_mat),
    GoF_down = gof_down,
    GoF_up = gof_up,
    LoF_down = lof_down,
    LoF_up = lof_up
  ))
}


# is_consistent_across <- function(count_table) {
#   
#   # Does a gene have nonzero counts as upregulated or downregulated in both
#   # LoF and GoF experiments?
#   
#   count_table$Consistent_across <- ifelse(
#     (count_table$GoF_up > 0 & count_table$LoF_up > 0) |
#       (count_table$GoF_down > 0 & count_table$LoF_down > 0),
#     FALSE,
#     TRUE
#   )
#   return(count_table)
# }


is_consistent_across <- function(count_table) {
  
  # Does a gene have nonzero counts as upregulated or downregulated in both
  # LoF and GoF experiments?
  
  gof_class <- ifelse(count_table$GoF_down > count_table$GoF_up, "GoF_down", "GoF_up")
  lof_class <- ifelse(count_table$LoF_down > count_table$LoF_up, "LoF_down", "LoF_up")
  
  count_table$Consistent_across <- ifelse(
    (gof_class == "GoF_down" & lof_class == "LoF_down") |
      (gof_class == "GoF_up" & lof_class == "LoF_up"),
    FALSE,
    TRUE
  )
  return(count_table)
}



# TODO: this doesn't address ties, which are just auto-assigned to down-reg...

add_purity <- function(count_table, use_total = "all") {
  
  # Treating GoF and LoF experiments as two classes, gives a measure [0,1] for
  # each gene of whether most experiments had the same up/down direction within
  # each class. Doesn't consider consistency (eg, purity can be 1 when all GoF
  # experiments are down and all LoF experiments are also down)
  # https://stats.stackexchange.com/questions/95731/how-to-calculate-purity
  # https://nlp.stanford.edu/IR-book/html/htmledition/evaluation-of-clustering-1.html
  
  stopifnot(use_total %in% c("all", "gene"))
  
  # are genes primarily up or down regulated in each of GoF/LoF experiments
  gof_class <- ifelse(count_table$GoF_down > count_table$GoF_up, "GoF_down", "GoF_up")
  lof_class <- ifelse(count_table$LoF_down > count_table$LoF_up, "LoF_down", "LoF_up")
  
  if (use_total == "all") {  # infer total fixed # of experiments
    
    all_msrd <- which(count_table$Count_NA == 0)[1]
    total <- sum(count_table[all_msrd, c("LoF_up", "LoF_down", "GoF_up", "GoF_down")])
    
    count_table$Purity <- vapply(1:nrow(count_table), function(x) {
      (count_table[x, gof_class[x]] +  count_table[x, lof_class[x]]) / total
    }, FUN.VALUE = numeric(1))
    
  } else {  # each gene has its own total of non-NAs
    
    total <- rowSums(count_table[, c("LoF_up", "LoF_down", "GoF_up", "GoF_down")])
    
    count_table$Purity <- vapply(1:nrow(count_table), function(x) {
      (count_table[x, gof_class[x]] +  count_table[x, lof_class[x]]) / total[x]
    }, FUN.VALUE = numeric(1))
    
  }
  return(count_table)
}


add_signed_purity <- function(count_table) {
  # If Consistent_across is FALSE (eg, dominant FC direction is same for LoF
  # and GoF), make Purity negative
  
  stopifnot(c("Consistent_across", "Purity") %in% colnames(count_table))
  
  count_table$Signed_purity <- 
    ifelse(count_table$Consistent_across, 
           count_table$Purity, 
           -(count_table$Purity))
  
  return(count_table)
}



avg_abs_fc <- function(fc_mat) {
  
  # Get the absolute average fold change for each gene in fc_mat as a df
  
  data.frame(Avg_abs_FC = rowMeans(abs(fc_mat), na.rm = TRUE)) %>% 
    rownames_to_column(var = "Symbol")
}


process_all <- function(fdr_mat, 
                        fdr, 
                        fc_mat, 
                        meta, 
                        purity_arg = "gene",
                        de_prior = NULL) {
  
  # Return a dataframe that contains for each symbol:
  # 1) DE counts/fraction measured/NA counts; 
  # 2) The count of times a gene was up/down in each of GoF and LoF experiments; 
  # 3) Whether the counts were consistent in direction between LoF and GoFs;
  # 4) The GoF/LoF class purity; 
  # 5) The average absolute FC
  
  
  de_count <- count_table(fdr_mat, fdr)
  
  direction <- tally_direction(fc_mat, meta) %>% 
    is_consistent_across()
  
  fc <- avg_abs_fc(fc_mat)
  
  df <- left_join(de_count, direction, by = "Symbol") %>% 
    left_join(fc, by = "Symbol") %>% 
    add_purity(use_total = purity_arg) %>% 
    add_signed_purity()
  
  if (!is.null(de_prior)) {
    df <- left_join(df, de_prior[, c("Symbol", "DE_Prior_Rank")], by = "Symbol")
  }
  
  return(df)
}


process_tf <- function(fdr_mat, 
                       fdr, 
                       fc_mat, 
                       meta, 
                       use_total = "all",
                       de_prior = NULL,
                       ortho = FALSE,
                       cores = 8) {
  
  # Performs process_all() for each unique symbol's experiments in meta/FC/FDR
  # and returns a list of a df for each symbol
  
  if (ortho) meta$Symbol <- str_to_title(meta$Symbol)
  
  tfs <- unique(meta$Symbol)
  
  tf_list <- mclapply(tfs, function(x) {
    
    tf_meta <- filter(meta, Symbol == x)
    
    process_all(fdr_mat = fdr_mat[, tf_meta$Experiment_ID], 
                fdr = fdr, 
                fc_mat = fc_mat[, tf_meta$Experiment_ID], 
                meta = tf_meta,
                de_prior = de_prior)
  }, mc.cores = cores)
  names(tf_list) <- tfs
  
  return(tf_list)
}


merge_ortho_counts <- function(count_ortho, count_hg, count_mm, pc_ortho) {
  
  # Given a list of DE counts for ortho and each species, return the ortho list
  # with two extra columns for the counts deriving from each species
  
  tfs <- str_to_title(c(names(count_ortho), names(count_hg), names(count_mm)))
  tfs <- unique(tfs)
  
  tf_list <- lapply(tfs, function(x) {
    
    mm <- count_mm[[x]][, c("Symbol", "Count_DE")] %>% 
      filter(Symbol %in% pc_ortho$Symbol_mm) %>% 
      dplyr::rename(Symbol_mm = Symbol, Count_DE_mm = Count_DE) %>% 
      left_join(pc_ortho, by = "Symbol_mm") %>% 
      dplyr::rename(Symbol = ID)
    
    hg <- count_hg[[str_to_upper(x)]][, c("Symbol", "Count_DE")] %>% 
      filter(Symbol %in% pc_ortho$Symbol_hg) %>% 
      dplyr::rename(Symbol_hg = Symbol, Count_DE_hg = Count_DE) %>% 
      left_join(pc_ortho, by = "Symbol_hg") %>% 
      dplyr::rename(Symbol = ID)
    
    left_join(count_ortho[[x]],
              mm[, c("Count_DE_mm", "Symbol")],
              by = "Symbol") %>%
      left_join(hg[, c("Count_DE_hg", "Symbol")],
                by = "Symbol")
  })
  names(tf_list) <- tfs
  return(tf_list)
}
