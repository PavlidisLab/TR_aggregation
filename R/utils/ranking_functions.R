## Code for precision-recall analysis of gene ranking's ability to recover
## curated targets.
## -----------------------------------------------------------------------------

library(tidyverse)
library(ROCR)
library(DescTools)
library(parallel)
library(data.table)  # frank() for rank using two variables


# Appends to df a column of the integer rank of perturbation evidence using the
# count of differentially expressed genes and the average absolute fold change. 

rank_perturb <- function(df) {
  
  stopifnot(c("Count_DE", "Avg_abs_FC") %in% colnames(df))
  
  df$Rank_perturbation <- data.table::frank(
    list(-df$Count_DE, -df$Avg_abs_FC), ties.method = "min")
  
  return(df)
}


# Appends to df a column of the integer rank of the binding evidence using the 
# mean binding score 

rank_binding <- function(df) {
  
  stopifnot("Mean_bind" %in% colnames(df))
  
  df$Rank_binding <- rank(-df$Mean_bind, ties.method = "min")
  
  return(df) 
}


# Add the rank product of the perturbation and binding evidence, as in the BETA 
# algo. Because "top ranked" genes from the RP have arbitrary small units, take
# the rank of the rank product to yield a final integer rank where rank=1 has 
# the strongest evidence.
# na_perturb arg controls behavior for genes with all NAs for perturbation
# data. If TRUE, these genes are treated as 'NAs' for rank product calculation 
# and thus will be assigned low importance. If FALSE, the original perturb 
# ranking is used, and the gene can still have elevated importance if the 
# binding ranking is high.
# BETA: https://pubmed.ncbi.nlm.nih.gov/24263090/
# R application: https://pubmed.ncbi.nlm.nih.gov/32894066/


rank_product <- function(df, na_perturb = FALSE) {
  
  stopifnot(
    c("Rank_binding", "Rank_perturbation", "Avg_abs_FC") %in% colnames(df))
  
  if (na_perturb) {
    
    n_not_na <- sum(!is.na(df$Avg_abs_FC))
    
    df <- df %>%
      mutate(
        Temp_rank = ifelse(is.na(Avg_abs_FC), NA, Rank_perturbation),
        Rank_prod = (Rank_binding / nrow(df)) * (Temp_rank / n_not_na),
        Rank_integrated = rank(Rank_prod, ties.method = "min")
      ) %>%
      dplyr::select(-c(Temp_rank, Rank_prod))
  
  } else {
      
    rank_prod <- (df$Rank_binding/nrow(df)) * (df$Rank_perturbation/nrow(df))
    df$Rank_integrated <- rank(rank_prod, ties.method = "min")
  }
  
  return(df)
}


# Append the perturbation, binding, and rank product ("Rank_integrated") to df 
# and order by rank integrated.

add_ranks <- function(df, ...) {
  
  df %>% 
    rank_binding() %>% 
    rank_perturb() %>% 
    rank_product(...) %>% 
    arrange(Rank_integrated)
  
}


# Uses the ROCR package to return a precision recall data frame where each entry
# has the calculated P+R for the top k ranked genes presence in the curated
# resource. Assumes that rows of rank_df are ordered by importance, such that
# the first row/gene has the highest importance.

get_pr <- function(rank_df) {
  
  stopifnot("Curated_target" %in% colnames(rank_df))
  
  # positives/has curated evidence=1 negatives/no known evidence=0 
  labels <- as.factor(as.numeric(rank_df$Curated_target))
  
  # convert ranks to monotonically decreasing scores representing rank importance
  scores <- 1/(1:length(labels))

  pred <- ROCR::prediction(predictions = scores, labels = labels)
  perf <- ROCR::performance(pred, measure = "prec", x.measure = "rec")
  
  pr_df <- data.frame(Precision = unlist(perf@y.values),
                      Recall = unlist(perf@x.values))
  
  return(pr_df)
}


# Uses the ROCR package to return return the area under the precision recall 
# curve for the ordered ranking of genes in rank_df

get_auprc <- function(rank_df) {
  
  stopifnot("Curated_target" %in% colnames(rank_df))
  
  # positives/has curated evidence=1 negatives/no known evidence=0 
  labels <- as.factor(as.numeric(rank_df$Curated_target))
  
  # convert ranks to monotonically decreasing scores representing rank importance
  scores <- 1/(1:length(labels))

  pred <- ROCR::prediction(predictions = scores, labels = labels)
  aupr <- ROCR::performance(pred, measure = "aucpr")@y.values[[1]]
  
  return(aupr)
}


# Get a list of PR data frames for the integrated, binding, and perturbation
# rankings in rank_df

get_all_pr <- function(rank_df) {
  
  stopifnot(c("Rank_integrated", "Rank_perturbation", "Rank_binding") %in% colnames(rank_df))
  
  # Order by respective rankings
  df_int <- arrange(rank_df, Rank_integrated)
  df_bind <- arrange(rank_df, Rank_binding)
  df_perturb <- arrange(rank_df, Rank_perturbation)
  
  list(
    Integrated = get_pr(df_int),
    Binding = get_pr(df_bind),
    Perturbation = get_pr(df_perturb)
  )
}


# Get a dataframe of AUPRCs for the integrated, binding, and perturbation
# rankings in rank_df

get_all_auprc <- function(rank_df) {
  
  stopifnot(c("Rank_integrated", "Rank_perturbation", "Rank_binding") %in% colnames(rank_df))
  
  # Order by respective rankings
  df_int <- arrange(rank_df, Rank_integrated)
  df_bind <- arrange(rank_df, Rank_binding)
  df_perturb <- arrange(rank_df, Rank_perturbation)
  
  data.frame(
    Integrated = get_auprc(df_int),
    Binding = get_auprc(df_bind),
    Perturbation = get_auprc(df_perturb)
  )
}


# Generates a vector of sampled AUPRC values, using the original gene rankings
# in rank_df and sampling size-matched targets from the data frame of curated 
# targets (lt_df)

sample_target_auprc <- function(rank_df, 
                                lt_df, 
                                reps = 1000, 
                                mouse = FALSE,
                                ncores = 1) {
  
  n_targets <- length(filter(rank_df, Curated_target)$Symbol)
  all_targets <- unique(str_to_upper(lt_df$Target_Symbol))
  
  stopifnot(n_targets > 0, length(all_targets) > 0)

  # Hacky case matching for mouse  
  if (mouse) {
    all_targets <- str_to_title(all_targets)
    all_targets <- all_targets[all_targets %in% rank_df$Symbol]
  }
  
  aupr <- mclapply(1:reps, function(x) {
    
    # Replace curated target (labels) with sampled targets
    sample_targets <- sample(all_targets, size = n_targets, replace = FALSE)
    rank_df$Curated_target <- rank_df$Symbol %in% sample_targets
    get_auprc(rank_df)
    
  }, mc.cores = ncores)
  
  return(unlist(aupr))
}


# Returns a table of the proportion of TR-specific sampled AUPRCs in 
# sampled_list that are greater than the observed AUPRCs in auprc_df

get_all_prop <- function(auprc_df, sampled_list) {
  
  # Inner function returns a single proportion
  get_prop <- function(obs, vec, n) { 
    sum(obs < vec) / n 
  }

  tfs <- rownames(auprc_df)
  stopifnot(identical(rownames(auprc_df), names(sampled_list)))
  
  # Loop through each TR, comparing the observed AUPRC with sampled
  prop_list <- lapply(tfs, function(x) {
    sampled <- unlist(sampled_list[[x]])
    n <- length(sampled)
    auprc <- auprc_df[x, ]
    apply(auprc, 2, function(y) get_prop(obs = y, vec = sampled, n))
  })
  names(prop_list) <- tfs

  return(do.call(rbind, prop_list))
}


# The following generates the AUPRC values for individual experiments and for 
# the rank product between each ChIP-seq and perturbation experiment.
# Returns a list that contains (1) these values in a df; (2) The observed
# AUPRCs for the aggregated rankings, and (3) The percentile of the aggregate 
# AUPRCs relative to the emprical distribution of individual experiment AUPRCs


# AUPRCs for each experiment in the ranked matrix

experiment_auprc <- function(ids, rank_mat, targets, ncores = 1) {
  
  mclapply(ids, function(x) {
    
    ranking <- names(sort(rank_mat[, x], na.last = TRUE))
    rank_df <- data.frame(Symbol = ranking, Curated_target = ranking %in% targets)
    get_auprc(rank_df)
  
  }, mc.cores = ncores)
}


# AUPRCs of the rank product of each paired ChIP-seq+perturb experiment

rp_auprc <- function(bind_ids, perturb_ids, bmat, pmat, targets, ncores) {
  
  mclapply(bind_ids, function(x) {
    
    lapply(perturb_ids, function(y) {
      
      rank_df <-
        data.frame(Symbol = rownames(bmat),
                   Bind = bmat[, x],
                   Perturb = pmat[, y]) %>%
        mutate(
          Rank_product = (Bind / nrow(bmat)) * (Perturb / nrow(pmat)),
          Curated_target = Symbol %in% targets
        ) %>%
        arrange(Rank_product)
      
      get_auprc(rank_df)
      
    })
  }, mc.cores = ncores)
}


# Main function: This is called within a lapply call that loops over tfs.
# Ranking_list is the list of TF-specific aggregated data gene rankings.

all_experiment_auprc <- function(tf,
                                 chip_meta,
                                 perturb_meta,
                                 bind_mat,
                                 perturb_mat,
                                 ranking_list,
                                 ncores = 1) {
  
  stopifnot(identical(rownames(bind_mat), rownames(perturb_mat)))
  
  # Isolate relevant TF-specific data
  bind_ids <- filter(chip_meta, Symbol == tf)$Experiment_ID
  perturb_ids <- filter(perturb_meta, Symbol == tf)$Experiment_ID
  targets <- filter(ranking_list[[tf]], Curated_target)$Symbol
  
  # Convert effect size matrices to ranks
  perturb_rank <- apply(perturb_mat[, perturb_ids, drop = FALSE], 2, rank, ties = "min")
  bind_rank <- apply(-bind_mat[, bind_ids, drop = FALSE], 2, rank, ties = "min")
  
  # AUPRCs for combined/perturb/bind aggregated rankings for the given TF
  all_auprcs <- get_all_auprc(ranking_list[[tf]])
  
  # AUPRCs for each specific ChIP-seq experiment
  bind_auprcs <- experiment_auprc(ids = bind_ids, 
                                  rank_mat = bind_rank,
                                  targets = targets, 
                                  ncores = ncores)
 
  # AUPRCs for each specific perturbation experiment
  perturb_auprcs <- experiment_auprc(ids = perturb_ids,
                                     rank_mat = perturb_rank,
                                     targets = targets,
                                     ncores = ncores)
  
  # AUPRCs of the rank product of each paired ChIP-seq+perturb experiment
  rp_auprcs <- 
    rp_auprc(bind_ids, perturb_ids, bind_rank, perturb_rank, targets, ncores)
  
  # Convert list of rank products into long df with paired experiment IDs
  rp_df <- reshape2::melt(rp_auprcs)
  rp_df$L1 <- colnames(bind_rank)[rp_df$L1]
  rp_df$L2 <- colnames(perturb_rank)[rp_df$L2]
  rp_df$ID <- paste(rp_df$L1, rp_df$L2, sep = ":")
  
  # Summary df of all individual AUPRCs
  auprc_df <- data.frame(
    Group = c(rep("Bind", length(bind_auprcs)),
              rep("Perturb", length(perturb_auprcs)),
              rep("Rank_product", nrow(rp_df))),
    ID = c(bind_ids, perturb_ids, rp_df$ID),
    AUPRC = c(unlist(bind_auprcs), unlist(perturb_auprcs), rp_df$value)
  )
  
  # Empirical distn of individual AUPRCs 
  auprc_emp <- ecdf(auprc_df$AUPRC)
  
  # Organize the percentile of observed aggregated AUPRC to the emp. distn
  auprc_perc <- data.frame(
    Integrated = auprc_emp(all_auprcs$Integrated),
    Binding = auprc_emp(all_auprcs$Binding),
    Perturbation = auprc_emp(all_auprcs$Perturbation))
  
  return(list(
    AUPRC_df = auprc_df,
    AUPRC_aggregated = all_auprcs,
    AUPRC_percentile = auprc_perc
  ))
}
