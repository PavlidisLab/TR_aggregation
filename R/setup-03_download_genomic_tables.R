## Download and formats gene/range annotation tables:
## 1) Protein coding genes (both ensembl and refseq)
## 2) List of TFs from Animal TF DB 
## 3) DIOPT 1:1 orthologous genes between mouse and human [NOTE: Pavlab server only]
## 4) ENCODE blacklisted regions (for filtering ChIP-seq peaks)
## 5) ENCODE candidate cis-reg elements (cCREs)
## 6) DE prior ranking [NOTE: Pavlab server only]
## 7) Literature curated targets from Chu 2021
## -----------------------------------------------------------------------------

library(tidyverse)
library(biomaRt)
source("R/setup-01_config.R")


# 1) Protein coding tables
# Ensembl has every TSS for a gene. Refseq select is only using one TSS per gene.
# All analysis was done on hg38/mm10. Fixing ensembl at V98.
# Using ensembl range formatting: no 'chr' prefix and strand as 1/-1
# https://www.ncbi.nlm.nih.gov/refseq/refseq_select/
# ------------------------------------------------------------------------------


download_refseq <- function(outfile, 
                            species) {  # Human|Mouse
  
  if (file.exists(outfile)) return(message(outfile, " already exists!"))
  
  if (species == "Mouse") {
    link <- "http://hgdownload.soe.ucsc.edu/goldenPath/mm10/database/ncbiRefSeqSelect.txt.gz"
    chr <- c(1:19, "MT", "X", "Y")
  } else if (species == "Human") {
    link <- "http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/ncbiRefSeqSelect.txt.gz"
    chr <- c(1:22, "MT", "X", "Y")
  }
  
  download.file(link, outfile)
  
  refseq <- read.delim(outfile, stringsAsFactors = FALSE, header = FALSE)
  
  refseq <- dplyr::select(refseq, c(V3, V5, V6, V4, V2, V13))
  
  colnames(refseq) <- c("Chromosome",
                        "Start",
                        "End",
                        "Strand",
                        "Refseq_ID",
                        "Symbol")
  
  refseq <- refseq %>% 
    mutate(
      Strand = ifelse(Strand == "+", 1, -1),
      Transcription_start_site = ifelse(Strand == 1, Start, End),
      Chromosome = str_replace(Chromosome, "chr", "")) %>% 
    filter(Chromosome %in% chr) %>% 
    arrange(match(Chromosome, chr), Transcription_start_site) %>% 
    dplyr::relocate(Transcription_start_site, .after = Chromosome)
  
  
  write.table(refseq,
              quote = FALSE,
              row.names = FALSE,
              sep = "\t",
              file = outfile)
  
}


download_ensembl_pcoding <- function(outfile,
                                     species,  # Human|Mouse
                                     version = "98") {
  
  if (file.exists(outfile)) return(message(outfile, " already exists!"))
  
  if (species == "Mouse") {
    symbol = "mgi_symbol"
    species_data = "mmusculus_gene_ensembl"
    chr_filter <- c(1:19, "MT", "X", "Y")
  } else if (species == "Human") {
    symbol = "hgnc_symbol"
    species_data = "hsapiens_gene_ensembl"
    chr_filter <- c(1:22, "MT", "X", "Y")
  }
  
  attributes <- c(
    "chromosome_name",
    "transcription_start_site",
    "transcript_start",
    "transcript_end",
    "strand",
    "ensembl_gene_id",
    symbol,
    "ucsc",
    "gene_biotype"
  )
  
  ens_mart <- useEnsembl(biomart = "ensembl",
                         dataset = species_data,
                         version = version)
  
  anno_table <- getBM(
    attributes = attributes,
    filters = "chromosome_name",
    values = chr_filter,
    mart = ens_mart,
    useCache = FALSE
  )
  
  # only protein coding gene type and order the table by chromosome then by TSS
  anno_table <- anno_table %>% 
    filter(gene_biotype == "protein_coding") %>%
    arrange(match(chromosome_name, chr_filter), transcription_start_site)
  
  anno_table$gene_biotype <- NULL
  
  colnames(anno_table) <- c(
    "Chromosome",
    "Transcription_start_site",
    "Start",
    "End",
    "Strand",
    "Gene_ID",
    "Symbol",
    "Transcript_ID"
  )
  
  write.table(anno_table,
              quote = FALSE,
              row.names = FALSE,
              sep = "\t",
              file = outfile)
  
}


download_refseq(outfile = ref_path_hg, species = "Human")
download_refseq(outfile = ref_path_mm, species = "Mouse")

# Getting time out issues - not using in final analysis so commented out for now
# download_ensembl_pcoding(outfile = ens_path_hg, species = "Human")
# download_ensembl_pcoding(outfile = ens_path_mm, species = "Mouse")


# 2) List of TFs http://bioinfo.life.hust.edu.cn/AnimalTFDB/
# ------------------------------------------------------------------------------


tf_url_hg <- "http://bioinfo.life.hust.edu.cn/static/AnimalTFDB3/download/Homo_sapiens_TF"
tf_url_mm <- "http://bioinfo.life.hust.edu.cn/static/AnimalTFDB3/download/Mus_musculus_TF"


if (!file.exists(tf_path_hg)) {
  download.file(url = tf_url_hg, destfile = tf_path_hg)
}


if (!file.exists(tf_path_mm)) {
  download.file(url = tf_url_mm, destfile = tf_path_mm)
}


# 3) Filter/organize high-confidence 1:1 orthologs between mouse and human.
# https://www.flyrnai.org/diopt
# Note that this was a provided data dump that lives on Pavlab servers.
# Recommended heuristic (passed on by Sanja): keep scores >= 5, require mutual 
# bestscore, and remove symbols with more than one match
# ------------------------------------------------------------------------------


diopt <- read.delim(diopt_path, stringsAsFactors = FALSE)

# Protein coding genes
pc_hg <- read.delim(ref_path_hg, stringsAsFactors = FALSE)
pc_mm <- read.delim(ref_path_mm, stringsAsFactors = FALSE)

symbol_hg <- unique(pc_hg$Symbol)
symbol_mm <- unique(pc_mm$Symbol)

diopt_filt <- filter(diopt,
                     score >= 5 &
                     human_symbol %in% pc_hg$Symbol &
                     symbol2 %in% pc_mm$Symbol &
                     ID_type2 == "MGI" &
                     best_score == "Yes" &
                     best_score_rev == "Yes")

split_hg <- split(diopt_filt, diopt_filt$human_symbol)
split_mm <- split(diopt_filt, diopt_filt$symbol2)

which_gt1_hg <- which(unlist(lapply(split_hg, nrow)) > 1)
which_gt1_mm <- which(unlist(lapply(split_mm, nrow)) > 1)

diopt_filt <- filter(diopt_filt, 
                        !(human_symbol) %in% names(which_gt1_hg) &
                          !(symbol2) %in% names(which_gt1_mm))

stopifnot(all(diopt_filt$human_symbol %in% pc_hg$Symbol))
stopifnot(all(diopt_filt$symbol2 %in% pc_mm$Symbol))

stopifnot(identical(n_distinct(diopt_filt$symbol2), 
                    n_distinct(diopt_filt$human_symbol)))

# create a df with a key as not all symbols have an exact 1:1 naming match

symbols <- data.frame(
  Symbol_hg = diopt_filt$human_symbol,
  Symbol_mm = diopt_filt$symbol2,
  ID = paste(diopt_filt$human_symbol, diopt_filt$symbol2, sep = "_")
)

stopifnot(identical(n_distinct(symbols$ID), nrow(symbols)))


if (!file.exists(ortho_path)) {
  write.table(symbols,
              sep = "\t",
              quote = FALSE,
              file = ortho_path)
}


# 4) ENCODE blacklisted regions
# ------------------------------------------------------------------------------


chr_hg <- c(1:22, "MT", "X", "Y")
chr_mm <- c(1:19, "MT", "X", "Y")

bl_url_hg <- "https://www.encodeproject.org/files/ENCFF356LFX/@@download/ENCFF356LFX.bed.gz"
bl_url_mm <- "https://www.encodeproject.org/files/ENCFF547MET/@@download/ENCFF547MET.bed.gz"


download_blacklist <- function(outfile, url, chr) {
  
  if (!file.exists(outfile)) {
    
    download.file(url, destfile = outfile)
    bl <- read.delim(outfile, header = FALSE)
    colnames(bl) <- c("Chromosome", "Start", "End")
    bl$Chromosome <- gsub("chr", "", bl$Chromosome)
    bl <- bl[order(match(bl$Chromosome, chr)),]
    
    write.table(bl,
                sep = "\t",
                quote = FALSE,
                row.names = FALSE,
                file = outfile)
  }
  
}


download_blacklist(bl_path_hg, bl_url_hg, chr_hg)
download_blacklist(bl_path_mm, bl_url_mm, chr_mm)



# 5) ENCODE cCREs
# https://screen.encodeproject.org/
# ------------------------------------------------------------------------------


ccre_url_hg <- "https://api.wenglab.org/screen_v13/fdownloads/V3/GRCh38-cCREs.bed"
ccre_url_mm <- "https://api.wenglab.org/screen_v13/fdownloads/V3/mm10-cCREs.bed"



download_ccre <- function(outfile, url, chr) {
  
  if (!file.exists(outfile)) {
    
    download.file(url, destfile = outfile)
    
    ccre <- read.delim(outfile, header = FALSE)
    colnames(ccre) <- c("Chromosome", "Start", "End", "Drop", "ID", "Group")
    ccre$Drop <- NULL # remove unknown ID column
    ccre$Chromosome <- gsub("chr", "", ccre$Chromosome)
    ccre <- ccre[order(match(ccre$Chromosome, chr)),]
    
    write.table(ccre,
                sep = "\t",
                quote = FALSE,
                row.names = FALSE,
                file = outfile)
  }
  
}


download_ccre(ccre_path_hg, ccre_url_hg, chr_hg)
download_ccre(ccre_path_mm, ccre_url_mm, chr_mm)


# 6) DE prior ranking (likelihood of being measured as DE across studies)
# Original PNAS version human only: https://doi.org/10.1073/pnas.1802973116
# Nathaniel updated/expanded version: unpublished at time of writing, but using 
# same algo. as in PNAS paper but more and diverse platforms
# ------------------------------------------------------------------------------


format_depr <- function(infile, outfile) {
  
  if (!file.exists(outfile)) {
    
    depr <- readRDS(infile)
    
    depr <- depr %>%
      dplyr::rename(Symbol = gene.Name,
                    DE_prior_rank = deStrict.mfxRank) %>%
      dplyr::select(Symbol, DE_prior_rank)
    
    write.table(depr,
                sep = "\t",
                quote = FALSE,
                row.names = FALSE,
                outfile)
  }
}


format_depr(depr_in_hg, depr_path_hg)
format_depr(depr_in_mm, depr_path_mm)


# Curated TR-tarters from Chu 2021. Records refers to experiments curated in
# the paper and includes additional information relative to all, which aggregates
# data from other curated resources like TRRUST.
# https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1009484
# ------------------------------------------------------------------------------

# S3: Experiments curated in study
chu2021_url_records <- "https://doi.org/10.1371/journal.pcbi.1009484.s025"

# S8: All records including external
chu2021_url_all <-  "https://doi.org/10.1371/journal.pcbi.1009484.s030"


if (!file.exists(chu2021_path_records)) {
  download.file(chu2021_url_records, destfile = chu2021_path_records)
}


if (!file.exists(chu2021_path_all)) {
  download.file(chu2021_url_all, destfile = chu2021_path_all)
}
