## Ready table for curation into google sheets. Meta is constructed as one
## experiment per row, but the corresponding resultsets may have multiple 
## relevant or irrelevant studies, thus requiring curation/matching. 
## -----------------------------------------------------------------------------

library(googlesheets4)
library(tidyverse)
source("R/setup-01_config.R")
source("R/utils/gemma_functions.R")


# Load curated perturbation metadata
meta <- read_sheet(ss = gsheets_perturb, 
                   sheet = "Master_batch1", 
                   trim_ws = TRUE, 
                   col_types = "c")

# Loaded means the experiment is in Gemma and resultsets have been downloaded
meta_loaded <- meta[meta$Results_Loaded == 1, ]


# Make experiment ids of form GSE_TF_Species_Perturbation. If duplicates, add 
# integer suffix

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


# Get resultset files for each GEO identifier
rs_files <- lapply(meta_loaded$GSE, list_result_sets, rs_dir)


# Link experiment IDs to result sets. A GEO directory may contain multiple 
# result set files beyond the one of interest, and a single file may also contain
# contrasts beyond the experiment of interest, this information must be matched 
# and curated.


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


# Loop through rs_df to open associated rs files and check colnames for
# match to expected symbol, or multi matches which are flagged for curation.

check_rs <- function(rs_df, rs_dir) {
  
  for (i in 1:nrow(rs_df)) {
    
    # Multiple resultset files found in the supplied GEO dir
    if (str_detect(rs_df$Resultset_ID[i], ";")) {  
      
      rs_ids <- unlist(str_split(rs_df$Resultset_ID[i], ";"))
      final_ids <- rs_ids
      final_cols <- list()
      
      # If rs is not relevant (no symbol match), remove the ID. Otherwise,
      # include all relevant cols for curation
      
      for (id in rs_ids) {
        
        rs <- try(load_result_set(GSE = rs_df$GSE[i],
                                  file = id,
                                  results_dir = rs_dir))
        
        cols <- match_colnames(rs, rs_df$Symbol[i])
        
        if (length(cols) == 0) {
          final_ids <- setdiff(final_ids, id)
        } else {
          final_cols <- c(final_cols, cols)
        }
        
      }

      rs_df$Resultset_ID[i] <- paste(final_ids, collapse = ";")
      rs_df$Match_Col[i] <- paste(unique(unlist(final_cols)), collapse = ";")
    
    # A single resultset in the GEO dir
    } else { 
        
      rs <- try(load_result_set(GSE = rs_df$GSE[i],
                                file = rs_df$Resultset_ID[i],
                                results_dir = rs_dir))
      
      cols <- paste(match_colnames(rs, rs_df$Symbol[i]), collapse = ";")
      
      rs_df$Match_Col[i] <- cols
      
    }
    
    if (any(str_detect(c(rs_df$Resultset_ID[i], rs_df$Match_Col[i]), ";"))) {
      rs_df$Curate[i] <- TRUE
    }
  }
  
  return(rs_df)
}


rs_df <- check_rs(rs_df, rs_dir)


# Order so those that need curation are at top and GSEs are clumped
rs_df <- arrange(rs_df, !Curate, GSE)


# save out sheet for curation and all meta of loaded

write_sheet(data = rs_df, 
            ss = gsheets_perturb, 
            sheet = paste0("Loaded_resultset_IDs_", date))


write.table(meta_loaded, 
            sep = "\t",
            quote = FALSE,
            row.names = FALSE,
            file = perturb_meta_path)
