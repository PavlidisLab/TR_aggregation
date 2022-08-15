## Code for importing and working with Gemma data. 
## TODO: Most of this was ported in from a quick first pass in 19_. Needs to be
## cleaned and tested

library(tidyverse)
library(gemmaAPI)


list_result_sets <- function(GSE,
                             results_dir = "~/Data/Expression_files/Gemma/Resultsets/") {
  # list the resultset .txt files for the input GSE
  # GSE : chr of form "GSEXXX"
  
  list.files(path = paste0(results_dir, GSE),
             pattern = "resultset")
}


load_result_set <- function(GSE, 
                            file,
                            results_dir = "~/Data/Expression_files/Gemma/Resultsets/") {
  # Read in the resultset table for the given GSE accession and file name.
  # Change gene symbol column name for consistency with peak tables
  # GSE : chr of form "GSEXXX"
  # file : chr, example "resultset_ID480379.data.txt"
  
  table <- read.delim(
    paste0(results_dir, GSE, "/", file),
    comment.char = "#",
    stringsAsFactors = FALSE,
    na.strings = "NA"
  )
  colnames(table)[colnames(table) == "Gene_Symbol"] <- "Symbol"
  warning_blank_table(table)
  return (table)
}


warning_blank_table <- function (rs) {
  # Generate a warning if the supplied resultset table has no rows or no
  # mapped symbols
  # rs : resultset dataframe downloaded from gemma
  
  if (nrow(rs) == 0) {
    warning("No rows in resultset table")
  } else if (all(rs$Symbol == "")) {
    warning("All mapped symbols are blank")
  }
}


match_colnames <- function(rs, symbol) {
  
  # Return unique colname patterns that match symbol. Helper to curate desired
  # columns of interest when multiple exist
  
  cols <- NA
  match <- str_detect(str_to_lower(colnames(rs)), str_to_lower(symbol))
  cols <- colnames(rs)[match]
  cols <- unique(str_replace_all(cols, "FoldChange_|Tstat_|PValue_", ""))
  return(cols)
}


which_symbol <- function(rs, symbol, exact = TRUE) {
  # retrieve the indices corresponding to the input symbol
  # rs : resultset dataframe downloaded from gemma
  # symbol : chr, example "TCF4"
  # exact : logical, if TRUE then return the exact pattern match to symbol,
  #         if FALSE allow for fuzzier matching
  
  stopifnot("Symbol" %in% colnames(rs))
  if (!exact) {
    which(str_detect(rs$Symbol, symbol))
  } else {
    which(str_detect(rs$Symbol, paste0("^", symbol, "$")))
  }
}


filter_symbol <- function(rs, symbol, exact = TRUE) {
  # extract the rows corresponding to the input symbol
  # rs : resultset dataframe downloaded from gemma
  # symbol : chr, example "TCF4"
  # exact : logical, if TRUE then return the exact pattern match to symbol,
  #         if FALSE allow for fuzzier matching
  
  
  stopifnot("Symbol" %in% colnames(rs))
  if (!exact) {
    rs[which_symbol(rs, symbol, exact = FALSE), ]
  } else {
    rs[which_symbol(rs, symbol), ]
  }
}


which_tstat_column <- function(rs) {
  # provide the index of the Tstat columns.
  # rs : resultset dataframe downloaded from gemma
  
  which(str_detect(
    string = colnames(rs), 
    pattern = "^Tstat"))
}


which_fc_column <- function(rs) {
  # provide the index of the fold change columns.
  # rs : resultset dataframe downloaded from gemma
  
  which(str_detect(
    string = colnames(rs), 
    pattern = "^FoldChange"))
}


which_pval_column <- function(rs) {
  # provide the index of the unadjusted pvalue columns.
  # rs : resultset dataframe downloaded from gemma
  
  which(str_detect(
    string = colnames(rs), 
    pattern = "^PValue"))
}


keep_match_cols <- function(rs, string) {
  # Return rs with only the descriptive columns and stat columns whose suffix  
  # match the provided string. 
  
  cols <- str_replace_all(colnames(rs), "Tstat_|PValue_|FoldChange_", "")
  match_cols <- which(cols == string)
  return(rs[, c(1:4, match_cols)])
}


strip_colnames <- function(rs) {
  # Remove the experiment info suffix from Pval/Tstat/FC
  
  stopifnot(ncol(rs) == 7)
  colnames(rs)[5:7] <- str_replace_all(colnames(rs)[5:7], "_.*", "")
  return(rs)
}


which_max_tstat <- function(rs) {
  # provide the index of the maximum absolute Tstat in the resultset table
  # rs : resultset dataframe downloaded from gemma
  
  which.max(abs(rs[, which_tstat_column(rs)]))
}


get_max_tstat <- function(rs) {
  # extract the row corresponding to the max T stat
  # rs : resultset dataframe downloaded from gemma
  
  rs[which_max_tstat(rs), ]
}


get_max_tstat_by_symbol <- function(rs, symbol, exact = TRUE) {
  # extract the row corresponding to the max T stat for the supplied symbol.
  # note that some symbols may have multiple mapped probes, and so this only
  # returns the absolute maximum tstat
  # rs : resultset dataframe downloaded from gemma
  # symbol : chr, example "TCF4"
  # exact : logical, if TRUE then return the exact pattern match to symbol,
  #         if FALSE allow for fuzzier matching
  
  if (!exact) {
    filter_table <- filter_symbol(rs, symbol, exact = FALSE)
  } else {
    filter_table <- filter_symbol(rs, symbol)
  }
  max_ix <- which_max_tstat(filter_table)
  max_probe <- filter_table[max_ix, "Element_Name"]
  rs[which(rs$Element_Name == max_probe), ]
}


which_max_fc <- function(rs) {
  # provide the index of the maximum absolute fold change in the resultset table
  # rs : resultset dataframe downloaded from gemma
  
  which.max(abs(rs[, which_fc_column(rs)]))
}


get_max_fc <- function(rs) {
  # extract the row corresponding to the max T stat
  # rs : resultset dataframe downloaded from gemma
  
  rs[which_max_fc(rs), ]
}


get_max_fc_by_symbol <- function(rs, symbol, exact = TRUE) {
  # extract the row corresponding to the max fold change for the supplied symbol.
  # note that some symbols may have multiple mapped probes, and so this only
  # returns the absolute maximum fold change
  # rs : resultset dataframe downloaded from gemma
  # symbol : chr, example "TCF4"
  # exact : logical, if TRUE then return the exact pattern match to symbol,
  #         if FALSE allow for fuzzier matching
  
  if (!exact) {
    filter_table <- filter_symbol(rs, symbol, exact = FALSE)
  } else {
    filter_table <- filter_symbol(rs, symbol)
  }
  max_ix <- which_max_fc(filter_table)
  max_probe <- filter_table[max_ix, "Element_Name"]
  rs[which(rs$Element_Name == max_probe), ]
}


get_rank_fc <- function(rs) {
  # Return a vector containing the ranks of the absolute fold change of the
  # input gemma resultset table. Higher rank = more differentially expressed
  # rs : resultset dataframe downloaded from gemma

  rank(abs(rs[, which_fc_column(rs)]))
}


get_perc_rank_fc <- function(rs) {
  # Return a vector containing the percentile ranks of the absolute fold change of the
  # input gemma resultset table. 1 = the max differentially expressed, 0 the min
  # rs : resultset dataframe downloaded from gemma
  
  fc_rank <- get_rank_fc(rs)
  round(fc_rank/length(fc_rank), 5)
}


keep_single_symbols <- function(rs) {
  # Remove probes/elements lacking a symbol or mapped to multiple symbols
  # rs : resultset dataframe downloaded from gemma
 
  cleaned_table <- filter(rs,
                          (!str_detect(Symbol, "\\|") &
                             str_detect(Symbol, boundary("character"))))
  warning_blank_table(cleaned_table)
  return(cleaned_table)
}


filter_by_max_tstat <- function(rs) {
  # Given a resultset table, find any symbols that have multiple probes,
  # and only keep the probe with the highest tstat. Take first row of highest 
  # fold change if there is a tie
  
  rs %>% 
    group_by(Symbol) %>% 
    filter(abs(Tstat) == max(abs(Tstat))) %>% 
    filter(abs(FoldChange) == max(abs(FoldChange))) %>% 
    dplyr::slice(1)
}






## OLD

# get_relevant_cols <- function(resultset_table, symbol) {
#   # given a gemma resultset table and a symbol, return only the columns of the 
#   # resultset table that are matched to the symbol
#   # resultset_table : resultset table downloaded from gemma
#   # symbol : chr, example "TCF4"
#   
#   result_names <- colnames(resultset_table)
#   gene_names_ix <- which(c("Element_Name", "Symbol", "Gene_Name", "NCBI_ID") %in% result_names)
#   
#   target_fc <- which(str_detect(
#     string = result_names,
#     pattern =  paste0("^FoldChange_", ".*", symbol, "\\.")
#   ))
#   
#   target_tstat <- which(str_detect(
#     string = result_names,
#     pattern =  paste0("^Tstat_", ".*", symbol, "\\.")
#   ))
#   
#   target_pval <- which(str_detect(
#     string = result_names,
#     pattern =  paste0("^PValue_", ".*", symbol, "\\.")
#   ))
#   
#   return(resultset_table[, c(gene_names_ix, target_fc, target_pval, target_tstat)])
# }

