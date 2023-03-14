## Generating an experiment x protein coding gene matrix, where elements are
## binding scores. After building for mouse and human, create mat of all
## experiments and orthologous genes.
## -----------------------------------------------------------------------------

library(tidyverse)
library(parallel)
source("R/setup-01_config.R")
source("R/utils/range_table_functions.R")

# Load metadata and link IDs to their corresponding directory
run_ids <- read.delim(chip_run_path, stringsAsFactors = FALSE)
meta <- read.delim(chip_meta_path, stringsAsFactors = FALSE)

stopifnot(all(meta$Experiment_ID %in% run_ids$Experiment_ID))
stopifnot(!any(run_ids$Dir == "" | is.na(run_ids$Dir)))

run_ids <- filter(run_ids, Experiment_ID %in% meta$Experiment_ID)
ids_hg <- run_ids[run_ids$Species == "Human", ]
ids_mm <- run_ids[run_ids$Species == "Mouse", ]

# Load protein coding tables for gene annotation
pc_hg <- read.delim(ref_path_hg, stringsAsFactors = FALSE)
pc_mm <- read.delim(ref_path_mm, stringsAsFactors = FALSE)
pc_ortho <- read.delim(ortho_path, stringsAsFactors = FALSE)

# Remove duplicated pseudoautosomal genes in human (keep X copy)
dupl <- pc_hg$Symbol[duplicated(pc_hg$Symbol)]
pc_hg <- filter(pc_hg, !(Symbol %in% dupl & Chromosome == "Y"))  

# Protein coding tables to GR objects
pc_hg <- pc_to_gr(pc_hg)
pc_mm <- pc_to_gr(pc_mm)

# Load ENCODE blacklists and convert to GRanges
bl_hg <- bl_to_gr(read.delim(bl_path_hg, stringsAsFactors = FALSE))
bl_mm <- bl_to_gr(read.delim(bl_path_mm, stringsAsFactors = FALSE))


# Into list:
# 1) Load individual peak table
# 2) Get peak summits
# 3) remove blacklisted
# 4) Get binding scores
# 5) Bind list into matrix
#-------------------------------------------------------------------------------



load_and_score <- function(input_df, bl_gr, pc_gr, cores) {
  
  scores_l <- lapply(1:nrow(input_df), function(x) {
    
    # MECP2 overlap peaks (histone pipeline), all others IDR peaks
    pipeline_type <- 
      ifelse(str_detect(input_df$Experiment_ID[x], "HISTONE"), "overlap", "idr")
    
    peak_gr <- 
      read_encpeak(input_df$Dir[x], peakset = pipeline_type) %>% 
      peak_to_gr() %>%
      filter_blacklist(bl_gr)
    
    scores <- binding_scores(pc_gr, 
                             peak_gr, 
                             method = "Ouyang",
                             decay_constant = 5e3,
                             ncore = cores)
    
    message(input_df$Experiment_ID[x], " complete ", Sys.time())
    return(scores)
  })
  
  mat <- do.call(cbind, scores_l)
  colnames(mat) <- input_df$Experiment_ID
  return(mat)
}


# Human
mat_hg <- load_and_score(ids_hg, bl_hg, pc_hg, cores)

# Mouse
mat_mm <- load_and_score(ids_mm, bl_mm, pc_mm, cores)


# Build ortho mat, using all experiments and 1:1 orthologs
#-------------------------------------------------------------------------------


symbol_ortho <- pc_ortho %>% 
  filter(Symbol_hg %in% rownames(mat_hg) & 
           Symbol_mm %in% rownames(mat_mm))
  
mat_ortho <- cbind(mat_hg[symbol_ortho$Symbol_hg, ], 
                   mat_mm[symbol_ortho$Symbol_mm, ])

mat_ortho <- mat_ortho[, meta$Experiment_ID]

rownames(mat_ortho) <- symbol_ortho$ID


stopifnot(identical(
  mat_hg["ACTB", "GSE76147_TCF4_Human_GEN2-2"],
  mat_ortho["ACTB_Actb", "GSE76147_TCF4_Human_GEN2-2"]
))


stopifnot(identical(
  mat_mm["Crybb1", "GSE132069_Pax6_Mouse_Control-KD"],
  mat_ortho["CRYBB1_Crybb1", "GSE132069_Pax6_Mouse_Control-KD"]
))


# Save 
#-------------------------------------------------------------------------------


saveRDS(
  object = mat_hg,
  file = bsmat_hg_path
)


saveRDS(
  object = mat_mm,
  file = bsmat_mm_path
)


saveRDS(
  object = mat_ortho,
  file = bsmat_ortho_path
)
