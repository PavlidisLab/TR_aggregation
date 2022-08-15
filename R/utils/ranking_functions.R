## Code for precision-recall analysis
## TODO: consider breaking up all exp AUPRC function
## TODO: collapse get_pr and get_vec_pr
## TODO: param doc
## -----------------------------------------------------------------------------

library(tidyverse)
library(ROCR)
library(DescTools)
library(parallel)


# The following uses the ROCR package to extract the area under the precision
# recall curve (AUPRC). The expected input is a data frame of TR-targets, ordered
# by gene importance. Also include functionality to generate sampled AUPRCs
# from an external df of curated targets.
# ------------------------------------------------------------------------------



get_pr <- function(df) {
  
  # Uses the ROCR package to return a data frame of precision and recall, where
  # each is calculated for every element/gene of df. 
  
  stopifnot("Curated_target" %in% colnames(df))
  
  # positives/has curated evidence=1 negatives/no known evidence=0 
  labels <- as.factor(as.numeric(df$Curated_target))
  
  # convert ranks to monotonically decreasing scores representing rank importance
  scores <- 1/(1:length(labels))
  
  pred <- ROCR::prediction(predictions = scores, labels = labels)
  perf <- ROCR::performance(pred, measure = "prec", x.measure = "rec")
  pr_df <- data.frame(Precision = unlist(perf@y.values), Recall = unlist(perf@x.values))
  return(pr_df)
}




get_auprc <- function(df) {
  
  # Uses the ROCR package to return the area under the precision recall curve
  # for the ordered ranking of gene in df as curated targets
  
  stopifnot("Curated_target" %in% colnames(df))
  
  # positives/has curated evidence=1 negatives/no known evidence=0 
  labels <- as.factor(as.numeric(df$Curated_target))
  
  # convert ranks to monotonically decreasing scores representing rank importance
  scores <- 1/(1:length(labels))
  
  pred <- ROCR::prediction(predictions = scores, labels = labels)
  aupr <- performance(pred, measure = "aucpr")@y.values[[1]]
  return(aupr)
}



get_all_auprc <- function(df) {
  
  # Get a dataframe of AUPRCs for the aggregated, binding, and perturbation
  # rankings in df
  
  stopifnot(c("Rank_integrated", "Rank_perturbation", "Rank_binding") %in% colnames(df))
  
  # Order by respective rankings
  df_int <- arrange(df, Rank_integrated)
  df_bind <- arrange(df, Rank_binding)
  df_perturb <- arrange(df, Rank_perturbation)
  
  data.frame(
    Integrated = get_auprc(df_int),
    Binding = get_auprc(df_bind),
    Perturbation = get_auprc(df_perturb)
  )
}


get_all_pr <- function(df) {
  
  # Get a list of PR data frames for the aggregated, binding, and perturbation
  # rankings in df
  
  stopifnot(c("Rank_integrated", "Rank_perturbation", "Rank_binding") %in% colnames(df))
  
  # Order by respective rankings
  df_int <- arrange(df, Rank_integrated)
  df_bind <- arrange(df, Rank_binding)
  df_perturb <- arrange(df, Rank_perturbation)
  
  list(
    Integrated = get_pr(df_int),
    Binding = get_pr(df_bind),
    Perturbation = get_pr(df_perturb)
  )
}



# The following generates sampled AUPRC values by using the original gene
# rankings (integrated by default) and sampling targets from the data frame
# of curated targets (lt_df)
# ------------------------------------------------------------------------------


sample_target_auprc <- function(df, 
                                lt_df, 
                                reps = 1000, 
                                mouse = FALSE,
                                ncore = 8) {
  
  # Get the AUPR for a random set of curated targets, equal in number to curated
  # targets in df, sampling from all available targets in lt_df
  
  n_targets <- length(filter(df, Curated_target)$Symbol)
  all_targets <- unique(str_to_upper(lt_df$Target_Symbol))
  
  stopifnot(n_targets > 0 || length(all_targets) > 0)
  
  # Hacky case matching for mouse  
  if (mouse) {
    all_targets <- str_to_title(all_targets)
    all_targets <- all_targets[all_targets %in% df$Symbol]
  }
  
  aupr <- mclapply(1:reps, function(x) {
    # Replace curated target (labels) with sampled targets
    sample_targets <- sample(all_targets, size = n_targets, replace = FALSE)
    df$Curated_target <- df$Symbol %in% sample_targets
    get_auprc(df)
  }, mc.cores = ncore)
  
  return(unlist(aupr))
}



# Proportion of sampled AUPRCS with greater AUPRC than observed 


get_all_prop <- function(auprc_df, sampled_list) {
  
  
  get_prop <- function(obs, vec, n) sum(obs < vec) / n
  
  tfs <- rownames(auprc_df)
  stopifnot(identical(rownames(auprc_df), names(sampled_list)))
  
  prop_list <- lapply(tfs, function(x) {
    sampled <- unlist(sampled_list[[x]])
    n <- length(sampled)
    auprc <- auprc_df[x, ]
    apply(auprc, 2, function(y) get_prop(obs = y, vec = sampled, n))
  })
  names(prop_list) <- tfs
  
  return(do.call(rbind, prop_list))
}



# The following generates the AUPRC values for individual experiments and
# their rank product pairings, returning a list that contains the values, as 
# well as the percentile of the aggregate AUPRCs relative to the distribution
# of individual experiment AUPRCs
# ------------------------------------------------------------------------------



get_vec_auprc <- function(vec, targets) {
  
  # Uses the ROCR package to return the area under the precision recall curve
  # for the ordered ranking of targets in vec
  
  # positives/has curated evidence=1 negatives/no known evidence=0 
  labels <- as.factor(as.numeric(vec %in% targets))
  
  # convert ranks to monotonically decreasing scores representing rank importance
  scores <- 1/(1:length(labels))
  
  pred <- ROCR::prediction(predictions = scores, labels = labels)
  aupr <- performance(pred, measure = "aucpr")@y.values[[1]]
  return(aupr)
}



# This function calculates the AUPRC of each individual experiment for the 
# requested TF, as well as the AUPRC of each individual pairing of ChIP-seq
# and perturbation experiments. 

all_experiment_auprc <- function(tf,
                                 chip_meta,
                                 perturb_meta,
                                 bind_mat,
                                 perturb_mat,
                                 ranking_list) {
  
  
  # isolate relevant TF-specific data
  bind_ids <- filter(chip_meta, Symbol == tf)$Experiment_ID
  perturb_ids <- filter(perturb_meta, Symbol == tf)$Experiment_ID
  targets <- filter(ranking_list[[tf]], Curated_target)$Symbol
  perturb_mat <- perturb_mat[, perturb_ids, drop = FALSE]
  bind_mat <- bind_mat[, bind_ids, drop = FALSE]
  
  
  # AUPRCs for combined/perturb/bind aggregated rankings
  all_auprcs <- get_all_auprc(ranking_list[[tf]])
  
  
  # AUPRCs for each specific ChIP-seq experiment
  bind_auprcs <- lapply(bind_ids, function(x) {
    ranking <- names(sort(bind_mat[, x], decreasing = TRUE))
    get_vec_auprc(ranking, targets)
  })
  
  
  # AUPRCs for each specific perturbation experiment
  perturb_auprcs <- lapply(perturb_ids, function(x) {
    ranking <- names(sort(perturb_mat[, x]))
    get_vec_auprc(ranking, targets)
  })
  
  
  # AUPRCs of the rank product of each paired experiment
  rp_auprcs <- mclapply(perturb_ids, function(x) {
    
    l <- lapply(bind_ids, function(y) {
      
      bind_ranking <- rank(-bind_mat[, y, drop = FALSE])
      perturb_ranking <- rank(perturb_mat[, x, drop = FALSE])
      
      df <- data.frame(Symbol = rownames(bind_mat),
                       Bind = bind_ranking, 
                       Perturb = perturb_ranking)
      
      df$Rank_product <- bind_ranking/nrow(df) * perturb_ranking/nrow(df)
      
      ranking <- arrange(df, Rank_product)$Symbol
      
      get_vec_auprc(ranking, targets)
      
    })
  }, mc.cores = 8)
  
  
  # Make rank product into long df with experiment IDs
  rp_df <- reshape2::melt(rp_auprcs)
  rp_df$L2 <- colnames(bind_mat)[rp_df$L2]
  rp_df$L1 <- colnames(perturb_mat)[rp_df$L1]
  rp_df$ID <- paste(rp_df$L1, rp_df$L2, sep = "_")
  
  # Summary df
  auprc_df <- data.frame(
    Group = c(rep("Bind", length(bind_auprcs)), 
              rep("Perturb", length(perturb_auprcs)),
              rep("Rank_product", nrow(rp_df))),
    ID = c(bind_ids, perturb_ids, rp_df$ID),
    AUPRC = c(unlist(bind_auprcs), unlist(perturb_auprcs), rp_df$value)
  )
  
  # empirical distn of AUPRC from all experiments and get percentiles of observed aggregate
  auprc_emp <- ecdf(auprc_df$AUPRC)
  
  auprc_perc <- data.frame(
    Aggregate = auprc_emp(all_auprcs$Integrated),
    Bind = auprc_emp(all_auprcs$Binding),
    Perturb = auprc_emp(all_auprcs$Perturbation))
  
  return(list(
    AUPRC_df = auprc_df,
    AUPRC_agg = all_auprcs,
    AUPRC_percentile = auprc_perc
  ))
}



# I originally wrote these to calculate AUPR down a list, but the ROCR version
# was faster (needed for the sampling procedure). Keeping for posterity
# ------------------------------------------------------------------------------


# get_pr <- function(ranking, targets) {
#   
#   # Ranking as a vector of ordered symbols, targets as positive labels. Returns
#   # a dataframe of precision and recall for each element in ranking (k)
#   
#   pos <- as.numeric(ranking %in% targets)
#   
#   pr <- mclapply(1:length(pos), function(k) {
#     
#     tp <- sum(pos[1:k])
#     fp <- k - tp
#     prec <- tp / (tp + fp) # TP/(TP+FP)
#     rec <- tp / length(targets) # TP/P
#     
#     if (is.na(prec)) prec <- 0
#     if (is.na(rec)) rec <- 0
#     
#     data.frame(Precision = prec, Recall = rec)
#     
#   }, mc.cores = 8)
#   
#   pr <- do.call(rbind, pr)
#   return(pr)
# }



# get_aucpr <- function(pr_df, method = "spline") {
#   
#   # Get the area under the precision-recall curve
#   # https://search.r-project.org/CRAN/refmans/DescTools/html/AUC.html
#   
#   suppressWarnings(DescTools::AUC(
#     x = pr_df$Recall,
#     y = pr_df$Precision,
#     method = method
#   ))
# }