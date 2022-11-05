## Saves out two sources of input information for the ENCODE ChIP-seq pipeline.
## The first is the list of SRX IDs for each TF, to be used to download the 
## corresponding fastq files. The second is a tsv containing which samples
## are associated with a Experiment_ID, which is the identifier used for the set
## of samples belonging to an experiment to be entered into the pipeline.
## -----------------------------------------------------------------------------

library(tidyverse)
library(googlesheets4)
source("R/setup-01_config.R")

# This is where text files organizing experiments for encode pipeline are saved
exp_out <- paste0(encode_dir, "inputs/")

# This is where text files organizing samples->fastqs are saved
fastq_out <- paste0(fastq_dir, "inputs/")

# Load the up to date curated metadata sheet

pipeline_meta <- read_sheet(
  ss = gsheets_chip,
  sheet = "Master_batch1",
  trim_ws = TRUE,
  col_types = "c"
)

pipeline_meta$Species <- str_replace(pipeline_meta$Species, "Mouse", "Mm")
pipeline_meta$Species <- str_replace(pipeline_meta$Species, "Human", "Hg")

# Assume that entries with a non-NA Experiment_ID have been curated

pipeline_meta <- pipeline_meta[!is.na(pipeline_meta$Experiment_ID), ]


# Save out input tsvs for fasdtq download
# ------------------------------------------------------------------------------


# removing problematic experiments - these would have already gone through the
# pipeline and were curated as unsuccessful. Coerce NAs (waiting for pipeline)
# to 0s to prevent NA in str detect

pipeline_meta$Complete[is.na(pipeline_meta$Complete)] <- 0
skip_exp <- str_detect(pipeline_meta$Complete, "Fastq failure")
pipeline_meta <- pipeline_meta[!skip_exp, ]

tfs <- unique(str_to_lower(pipeline_meta$Symbol))

# for every TF, split human and mouse meta, and create a table of SRX IDs
# and whether or not the sample is an input control. save out

# TODO: This was a hacky start that needs to be fleshed out. By TF/species is
# agitating for newer, smaller inputs, and the overall workflow should be 
# better attached to the pipeline

for (tf in tfs) {
  
  meta_hg <- filter(pipeline_meta, Symbol == str_to_upper(tf))
  
  if (nrow(meta_hg) > 0) {
    
    output_hg <- data.frame(
      SRX_ID = unique(c(meta_hg$ID, unlist(str_split(meta_hg$Input_ID, ", ")))),
      Species = meta_hg$Species[1],
      Symbol = meta_hg$Symbol[1])
    
    output_hg <- output_hg[output_hg$SRX_ID != "NA", ]
    output_hg$Control <- !(output_hg$SRX_ID %in% meta_hg$ID)
    output_hg$Species <- str_replace(output_hg$Species, "Human", "Hg")
    
    write.table(
      x = output_hg,
      quote = FALSE,
      row.names = FALSE,
      col.names = FALSE,
      sep = "\t",
      file = paste0(fastq_dir, tf, "_hg_fastq_dl_input.tsv")
    )
  }
  
  meta_mm <- filter(pipeline_meta, Symbol == str_to_title(tf))
  
  if (nrow(meta_mm) > 0) {
    
    output_mm <- data.frame(
      SRX_ID = unique(c(meta_mm$ID, unlist(str_split(meta_mm$Input_ID, ", ")))),
      Species = meta_mm$Species[1],
      Symbol = meta_mm$Symbol[1]
    )
    
    output_mm <- output_mm[output_mm$SRX_ID != "NA", ]
    output_mm$Control <- !(output_mm$SRX_ID %in% meta_mm$ID)
    output_mm$Species <- str_replace(output_mm$Species, "Mouse", "Mm")
    
    write.table(
      x = output_mm,
      quote = FALSE,
      row.names = FALSE,
      col.names = FALSE,
      sep = "\t",
      file = paste0(fastq_dir, tf, "_mm_fastq_dl_input.tsv")
    )
  }
}


# Save out run inputs for ENCODE pipeline
# ------------------------------------------------------------------------------


pipeline_runs <- unique(pipeline_meta$Experiment_ID)


for (run in pipeline_runs) {
  
  out_file = paste0(exp_out, run, ".tsv")
  
  if (!file.exists(out_file)) {
    
    run_table <- filter(pipeline_meta, Experiment_ID == run)
    
    # str split for multiple inputs
    exp_ids <- run_table$ID
    input_ids <- unique(unlist(str_split(run_table$Input_ID, pattern = ", ")))
    
    run_table <- data.frame(
      SRX_ID = c(exp_ids, input_ids),
      Species = run_table$Species[1],
      Symbol = run_table$Symbol[1]
    )
    
    run_table <- run_table[run_table$SRX_ID != "NA", ]
    
    run_table$Control <- run_table$SRX_ID %in% input_ids
    
    write.table(
      x = run_table,
      quote = FALSE,
      row.names = FALSE,
      col.names = FALSE,
      sep = "\t",
      file = out_file
    )
  }
}
