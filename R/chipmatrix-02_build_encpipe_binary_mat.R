## Generating a binary experiment x protein coding gene matrix, where elements 
## indicate whether or not a peak was within input distance to the TSS
## -----------------------------------------------------------------------------

library(tidyverse)
library(parallel)
source("~/regnetR/R/utils/range_table_functions.R")

date <- "Apr2022"
dist <- 25e3
unique_symbols_flag <- TRUE  # every row of pc anno table or just unique symbols
peakset <- "idr"  # idr for TF (but overlap still for mecp2) or all overlap

pipeline_dir <- "/cosmos/data/pipeline-output/chipseq-encode-pipeline/chip"
out_dir <- "~/Data/Annotated_objects/Bind_matrices/Encpipe/"

pc_hg <- read.delim("~/Data/Metadata/refseq_select_hg38.tsv", stringsAsFactors = FALSE)
pc_mm <- read.delim("~/Data/Metadata/refseq_select_mm10.tsv", stringsAsFactors = FALSE)
pc_ortho <- read.delim("~/Data/Metadata/hg_mm_1to1_ortho_genes_DIOPT-v8.tsv", stringsAsFactors = FALSE)

bl_hg <- read.delim("~/Data/Chromosome_info/blacklist_hg38.tsv", stringsAsFactors = FALSE)
bl_mm <- read.delim("~/Data/Chromosome_info/blacklist_mm10.tsv", stringsAsFactors = FALSE)

# Need to link meta/experiment IDs with corresponding directory
run_ids <- read.delim(paste0("~/Data/Metadata/Chipseq/batch1_run_dirs_", date, ".tsv"), stringsAsFactors = FALSE)
meta <- read.delim(paste0("~/Data/Metadata/Chipseq/batch1_chip_meta_final_", date, ".tsv"), stringsAsFactors = FALSE)
stopifnot(all(meta$Experiment_ID %in% run_ids$Experiment_ID))
run_ids <- filter(run_ids, Experiment_ID %in% meta$Experiment_ID)

stopifnot(!any(run_ids$Dir == "" | is.na(run_ids$Dir)))

ids_hg <- run_ids[run_ids$Species == "Human", ]
ids_mm <- run_ids[run_ids$Species == "Mouse", ]

# Remove duplicated pseudoautosomal genes (keep X copy)
dupl <- pc_hg$Symbol[duplicated(pc_hg$Symbol)]
pc_hg <- filter(pc_hg, !(Symbol %in% dupl & Chromosome == "Y"))

# TODO: hacky fix from ensembl->refseq, need to re-write anno function
pc_hg <- pc_hg %>% 
  mutate(Gene_ID = str_replace(Refseq_ID, "_", ""), Transcript_ID = NA)

pc_mm <- pc_mm %>% 
  mutate(Gene_ID = str_replace(Refseq_ID, "_", ""), Transcript_ID = NA)


# Call distance anno for every run, which produces a binary vector of gene IDs
# or unique symbol, filled if a peak was found within the input distance. 
# Bind these vectors into a matrix
#-------------------------------------------------------------------------------


# Human

bind_list_hg <- lapply(1:nrow(ids_hg), function(i) {
  
  # hacky: check if mecp2 histone run
  
  pipeline_type <- 
    ifelse(str_detect(ids_hg$Experiment_ID[i], "HISTONE"), "overlap", peakset)
  
  peak_table <- read_encpeak(ids_hg$Dir[i], peakset = pipeline_type)
  
  bind_vec <- distance_anno(
    peak_table = peak_table,
    gene_table = pc_hg,
    bl_table = bl_hg,
    distance = dist,
    unique_symbols = unique_symbols_flag
  )

  message(ids_hg$Experiment_ID[i], " complete ", Sys.time())
  
  return(bind_vec)

})

mat_hg <- do.call(cbind, bind_list_hg)
colnames(mat_hg) <- ids_hg$Experiment_ID


# Mouse

bind_list_mm <- lapply(1:nrow(ids_mm), function(i) {
  
  # hacky: check if mecp2 histone run
  
  pipeline_type <- 
    ifelse(str_detect(ids_mm$Experiment_ID[i], "HISTONE"), "overlap", peakset)
  
  peak_table <- read_encpeak(ids_mm$Dir[i], peakset = pipeline_type)
  
  bind_vec <- distance_anno(
    peak_table = peak_table,
    gene_table = pc_mm,
    bl_table = bl_mm,
    distance = dist,
    unique_symbols = unique_symbols_flag
  )
  
  message(ids_mm$Experiment_ID[i], " complete ", Sys.time())
  
  return(bind_vec)
  
})

mat_mm <- do.call(cbind, bind_list_mm)
colnames(mat_mm) <- ids_mm$Experiment_ID



# Build ortho mat, using all experiments and 1:1 orthologs
#-------------------------------------------------------------------------------


subset_to_ortho <- function(mat, pc_vec) {
  # subset and order rows of mat to vec of ortho protein coding genes
  mat <- mat[rownames(mat) %in% pc_vec, ]
  mat <- mat[order(match(rownames(mat), pc_vec)), ]
  return(mat)
}


mat_ortho <- cbind(subset_to_ortho(mat_hg, pc_ortho$Symbol_hg),
                   subset_to_ortho(mat_mm, pc_ortho$Symbol_mm))

mat_ortho <- mat_ortho[, meta$Experiment_ID]
rownames(mat_ortho) <- pc_ortho$ID


stopifnot(identical(
  mat_hg["ACTB", "GSE76147_TCF4_Human_GEN2-2"],
  mat_ortho["ACTB_Actb", "GSE76147_TCF4_Human_GEN2-2"]
))


stopifnot(identical(
  mat_mm["Crybb1", "GSE132069_Pax6_Mouse_Control-KD"],
  mat_ortho["CRYBB1_Crybb1", "GSE132069_Pax6_Mouse_Control-KD"]
))



# Save 
# ------------------------------------------------------------------------------


saveRDS(
  object = mat_hg,
  file = paste0(out_dir, "binary_refseq_human_batch1_", date, "_", "distance=", dist/1e3, "kb_uniquesymbol=", unique_symbols_flag, "_peakset=", peakset, ".RDS")
)


saveRDS(
  object = mat_mm,
  file = paste0(out_dir, "binary_refseq_mouse_batch1_", date, "_", "distance=", dist/1e3, "kb_uniquesymbol=", unique_symbols_flag, "_peakset=", peakset, ".RDS")
)


saveRDS(
  object = mat_ortho,
  file = paste0(out_dir, "binary_refseq_ortho_batch1_", date, "_", "distance=", dist/1e3, "kb_uniquesymbol=", unique_symbols_flag, "_peakset=", peakset, ".RDS")
)
