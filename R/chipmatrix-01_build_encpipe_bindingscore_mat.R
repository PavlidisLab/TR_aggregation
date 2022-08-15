## Generating an experiment x protein coding gene matrix, where elements are
## binding scores. After building for mouse and human, create mat of all
## experiments and orthologous genes.
## TODO: bindscore functionality re-write from multi TSS ensembl to refseq
## -----------------------------------------------------------------------------

library(tidyverse)
library(parallel)
source("~/regnetR/R/utils/range_table_functions.R")

date <- "Apr2022"
decay_constant <- 5000
intensity_flag <- FALSE  # should peak intensity scores be included?
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


# build the matrices. 
# all == all transcripts, max only take max score for a symbol
# for all, rownames are of form: Symbol_Chr_Start_Stop_TranscriptID

bscore_mat_all_hg <- matrix(0, nrow = nrow(pc_hg), ncol = nrow(ids_hg))
bscore_mat_max_hg <- matrix(0, nrow = n_distinct(pc_hg$Symbol), ncol = nrow(ids_hg))

bscore_mat_all_mm <- matrix(0, nrow = nrow(pc_mm), ncol = nrow(ids_mm))
bscore_mat_max_mm <- matrix(0, nrow = n_distinct(pc_mm$Symbol), ncol = nrow(ids_mm))

colnames(bscore_mat_all_hg) <- colnames(bscore_mat_max_hg) <- ids_hg$Experiment_ID
rownames(bscore_mat_all_hg) <- get_gene_names(pc_hg)
rownames(bscore_mat_max_hg) <- unique(pc_hg$Symbol)

colnames(bscore_mat_all_mm) <- colnames(bscore_mat_max_mm) <- ids_mm$Experiment_ID
rownames(bscore_mat_all_mm) <- get_gene_names(pc_mm)
rownames(bscore_mat_max_mm) <- unique(pc_mm$Symbol)

chr_order <- unique(pc_hg$Chromosome)


# 1) Load individual peak table. 
# 2) Get peak summits
# 3) remove blacklisted
# 4) Get binding scores
# 5) Fill matrices row-wise
#-------------------------------------------------------------------------------


# Human


bl_gr_hg <- bl_to_gr(bl_hg)


for (i in 1:nrow(ids_hg)) {

  # load mecp2 overlap peakset, IDR or overlap (peakset arg) for TFs
  
  pipeline_type <- 
    ifelse(str_detect(ids_hg$Experiment_ID[i], "HISTONE"), "overlap", peakset)
  
  peak_table <- read_encpeak(ids_hg$Dir[i], peakset = pipeline_type)
  
  # match format to pcoding, remove blacklisted regions, get peak summit position
  # hacky - currently blacklist filter requires GR object but ouyang assignment
  # works with a df. convert to GR for filtering and revert back
  
  peak_table <- peak_table %>% 
    peak_to_gr() %>% 
    filter_blacklist(bl_gr = bl_gr_hg) %>% 
    gr_to_df()
  
  
  # call score for every transcript, and take max score for each symbol
  bscore_scores_all_hg <- all_ouyang_scores(peak_table,
                                            pc_hg,
                                            decay_constant,
                                            cores = 8,
                                            use_intensity = intensity_flag)
  
  bscore_scores_max_hg <- max_score_per_symbol(bscore_scores_all_hg, pc_hg)
  
  stopifnot(identical(names(bscore_scores_all_hg), rownames(bscore_mat_all_hg)))
  stopifnot(identical(names(bscore_scores_max_hg), rownames(bscore_mat_max_hg)))
  
  bscore_mat_all_hg[ , i] <- bscore_scores_all_hg
  bscore_mat_max_hg[ , i] <- bscore_scores_max_hg
  
  message(ids_hg$Experiment_ID[i], " complete ", Sys.time())
}


# Mouse


bl_gr_mm <- bl_to_gr(bl_mm)


for (i in 1:nrow(ids_mm)) {
  
  # load mecp2 overlap peakset, IDR for TFs
  
  pipeline_type <- 
    ifelse(str_detect(ids_mm$Experiment_ID[i], "HISTONE"), "overlap", peakset)
  
  peak_table <- read_encpeak(ids_mm$Dir[i], peakset = pipeline_type)
  
  # match format to pcoding, remove blacklisted regions, get peak summit position
  # hacky - currently blacklist filter requires GR object but ouyang assignment
  # works with a df. convert to GR for filtering and revert back
  
  peak_table <- peak_table %>% 
    peak_to_gr() %>% 
    filter_blacklist(bl_gr = bl_gr_mm) %>% 
    gr_to_df()
  
  # call score for every transcript, and take max score for each symbol
  bscore_scores_all_mm <- all_ouyang_scores(peak_table, 
                                            pc_mm, 
                                            decay_constant, 
                                            cores = 8,
                                            use_intensity = intensity_flag)
  
  bscore_scores_max_mm <- max_score_per_symbol(bscore_scores_all_mm, pc_mm)
  
  stopifnot(identical(names(bscore_scores_all_mm), rownames(bscore_mat_all_mm)))
  stopifnot(identical(names(bscore_scores_max_mm), rownames(bscore_mat_max_mm)))
  
  bscore_mat_all_mm[ , i] <- bscore_scores_all_mm
  bscore_mat_max_mm[ , i] <- bscore_scores_max_mm
  
  message(ids_mm$Experiment_ID[i], " complete ", Sys.time())
}



# Build ortho mat, using all experiments and 1:1 orthologs
#-------------------------------------------------------------------------------


subset_to_ortho <- function(mat, pc_vec) {
  # subset and order rows of mat to vec of ortho protein coding genes
  mat <- mat[rownames(mat) %in% pc_vec, ]
  mat <- mat[order(match(rownames(mat), pc_vec)), ]
  return(mat)
}


mat_ortho <- cbind(subset_to_ortho(bscore_mat_max_hg, pc_ortho$Symbol_hg),
                   subset_to_ortho(bscore_mat_max_mm, pc_ortho$Symbol_mm))

mat_ortho <- mat_ortho[, meta$Experiment_ID]
rownames(mat_ortho) <- pc_ortho$ID


stopifnot(identical(
  bscore_mat_max_hg["ACTB", "GSE76147_TCF4_Human_GEN2-2"],
  mat_ortho["ACTB_Actb", "GSE76147_TCF4_Human_GEN2-2"]
))


stopifnot(identical(
  bscore_mat_max_mm["Crybb1", "GSE132069_Pax6_Mouse_Control-KD"],
  mat_ortho["CRYBB1_Crybb1", "GSE132069_Pax6_Mouse_Control-KD"]
))


# Save 
#-------------------------------------------------------------------------------


saveRDS(
  object = bscore_mat_max_hg,
  file = paste0(out_dir, "ouyang_refseq_human_batch1_", date, "_", "dc=", decay_constant, "_intensity=", intensity_flag, "_peakset=", peakset, ".RDS")
)


saveRDS(
  object = bscore_mat_max_mm,
  file = paste0(out_dir, "ouyang_refseq_mouse_batch1_", date, "_", "dc=", decay_constant, "_intensity=", intensity_flag, "_peakset=", peakset, ".RDS")
)


saveRDS(
  object = mat_ortho,
  file = paste0(out_dir, "ouyang_refseq_ortho_batch1_", date, "_", "dc=", decay_constant, "_intensity=", intensity_flag, "_peakset=", peakset, ".RDS")
)