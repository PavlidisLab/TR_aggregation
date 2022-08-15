## Ready table for curation into google sheets. Meta is constructed as one
## experiment per row, but the corresponding resultsets may have multiple 
## relevant or irrelevant studies, thus requiring curation/matching. 
## -----------------------------------------------------------------------------

library(googlesheets4)
library(tidyverse)
source("~/regnetR/R/utils/gemma_functions.R")

gsheets_id <- "1oXo1jfoPYcX94bq2F6Bqm1Es1glc8g9mnJvYAO37Vw4"
date <- "Apr2022"
exprs_dir <- "~/Data/Expression_files/Gemma/Resultsets/"
meta_loaded_out <- paste0("~/Data/Metadata/Gemma/batch1_tfperturb_meta_final_", date, ".tsv")

# master metadata
meta <- read_sheet(ss = gsheets_id, 
                   sheet = "Master_batch1", 
                   trim_ws = TRUE, 
                   col_types = "c")

meta_loaded <- meta[meta$Results_Loaded == 1, ]

# ids of form GSE_TF_Species_Perturbation - if duplicates, add integer suffix

ids <- paste(str_replace(meta_loaded$GSE, "\\.[:digit:]", ""),
             meta_loaded$Symbol,
             meta_loaded$Species,
             meta_loaded$Perturbation,
             sep = "_")

ids <- make.unique(ids, sep = "-")


meta_loaded <- meta_loaded %>% 
  dplyr::select(-c(Gemma_Link, GEO_Link, Results_Loaded)) %>% 
  mutate(Experiment_ID = ids) %>% 
  arrange(Symbol)


rs_files <- lapply(meta_loaded$GSE, function(x) {
  list.files(path = paste0(exprs_dir, x), pattern = "resultset")
})


# Link IDs to result sets


check_rs <- function(rs_df) {
  
  for (i in 1:nrow(rs_df)) {
    # Loop through rs_df to open associated rs files and check colnames for
    # match to expected symbol, or multi matches which are flagged for curation.
      
    if (str_detect(rs_df$Resultset_ID[i], ";")) {  # multi rs in a GEO dir
      
      rs_ids <- unlist(str_split(rs_df$Resultset_ID[i], ";"))
      final_ids <- rs_ids
      final_cols <- list()
      
      # If rs is not relevant (no symbol match), remove the ID. Otherwise,
      # include all relevant cols for curation
      
      for (id in rs_ids) {
        rs <- try(load_result_set(GSE = rs_df$GSE[i], file = id))
        cols <- match_colnames(rs, rs_df$Symbol[i])
        if (length(cols) == 0) {
          final_ids <- setdiff(final_ids, id)
        } else {
          final_cols <- c(final_cols, cols)
        }
      }
      
      rs_df$Resultset_ID[i] <- paste(final_ids, collapse = ";")
      rs_df$Match_Col[i] <- paste(unique(unlist(final_cols)), collapse = ";")
    
    } else { # single resultset in a GEO dir
        
      rs <- load_result_set(GSE = rs_df$GSE[i], file = rs_df$Resultset_ID[i])
      cols <- paste(match_colnames(rs, rs_df$Symbol[i]), collapse = ";")
      rs_df$Match_Col[i] <- cols
    }
    
    if (any(str_detect(c(rs_df$Resultset_ID[i], rs_df$Match_Col[i]), ";"))) {
      rs_df$Curate[i] <- TRUE
    }
  }
  return(rs_df)
}


rs_df <- data.frame(
  GSE = meta_loaded$GSE,
  Symbol = meta_loaded$Symbol,
  Curate = FALSE,
  Resultset_ID = unlist(lapply(rs_files, paste, collapse = ";")), # ';' as delim
  Experiment_ID = meta_loaded$Experiment_ID,
  Match_Col = NA,
  Cell_Type = meta_loaded$Cell_Type,
  Notes = meta_loaded$Note,
  stringsAsFactors = FALSE
)


rs_df <- check_rs(rs_df)


# Order so those that need curation are at top and GSEs are clumped
rs_df <- arrange(rs_df, !Curate, GSE)

# save out sheet for curation and all meta of loaded

write_sheet(data = rs_df, 
            ss = gsheets_id, 
            sheet = paste0("Loaded_resultset_IDs_", date))

write.table(meta_loaded, 
            sep = "\t",
            quote = FALSE,
            row.names = FALSE,
            file = meta_loaded_out)
