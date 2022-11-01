## Download and formats gene/range annotation tables:
## 1) Protein coding genes (both ensembl and refseq)
## 2) List of TFs from Animal TF DB 
## 3) DIOPT 1:1 orthologous genes between mouse and human [NOTE: Pavlab server only]
## 4) ENCODE blacklisted regions (for filtering ChIP-seq peaks)
## 5) ENCODE candidate cis-reg element (cCREs)
## 6) DE prior ranking [NOTE: Pavlab server only]
## -----------------------------------------------------------------------------

library(tidyverse)
library(biomaRt)

# Out files

ref_hg <- "~/Data/Metadata/refseq_select_hg38.tsv"
ens_hg <- "~/Data/Metadata/ensembl_human_protein_coding_V98.tsv"
ref_mm <- "~/Data/Metadata/refseq_select_mm10.tsv"
ens_mm <- "~/Data/Metadata/ensembl_mouse_protein_coding_V98.tsv"

tf_hg <- "~/Data/Metadata/human_tfs.tsv"
tf_mm <- "~/Data/Metadata/mouse_tfs.tsv"

ortho <- "~/Data/Metadata/hg_mm_1to1_ortho_genes_DIOPT-v8.tsv"

bl_hg <- "~/Data/Chromosome_info/blacklist_hg38.tsv"
bl_mm <- "~/Data/Chromosome_info/blacklist_mm10.tsv"

ccre_hg <- "~/Data/Chromosome_info/cCREs_V3_hg38.bed"
ccre_mm <- "~/Data/Chromosome_info/cCREs_V3_mm10.bed"

depr_out_hg <- "~/Data/Metadata/DE_prior_hg.tsv"
depr_out_mm <- "~/Data/Metadata/DE_prior_mm.tsv"
depr_out_pnas <- "~/Data/Metadata/DE_prior_PNAS.tsv"


# 1) Protein coding tables
# Ensembl has every TSS for a gene. Refseq select is only using one TSS per gene.
# All analysis was done on hg38/mm10. Fixing ensembl at V98.
# Using ensembl range formatting: no 'chr' prefix and strand as 1/-1
# https://www.ncbi.nlm.nih.gov/refseq/refseq_select/
# ------------------------------------------------------------------------------


download_refseq <- function(outfile, 
                            species) {  # Human|Mouse
  
  if (file.exists(outfile)) return(message(outfile, " already exists!"))
  
  if(species == "Mouse") {
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
  
  if(species == "Mouse") {
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
    filter = "chromosome_name",
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


download_refseq(outfile = ref_hg, species = "Human")
download_refseq(outfile = ref_mm, species = "Mouse")

download_ensembl_pcoding(outfile = ens_hg, species = "Human")
download_ensembl_pcoding(outfile = ens_mm, species = "Mouse")


# 2) List of TFs http://bioinfo.life.hust.edu.cn/AnimalTFDB/
# ------------------------------------------------------------------------------


tf_url_hg <- "http://bioinfo.life.hust.edu.cn/static/AnimalTFDB3/download/Homo_sapiens_TF"
tf_url_mm <- "http://bioinfo.life.hust.edu.cn/static/AnimalTFDB3/download/Mus_musculus_TF"


if (!file.exists(tf_hg)) {
  download.file(url = tf_url_hg, destfile = tf_hg)
}


if (!file.exists(tf_mm)) {
  download.file(url = tf_url_mm, destfile = tf_mm)
}


# 3) Filter/organize high-confidence 1:1 orthologs between mouse and human
# https://www.flyrnai.org/diopt
# Recommended heuristic (passed on by Sanja): keep scores >= 5, require mutual 
# bestscore, and remove symbols with more than one match
# ------------------------------------------------------------------------------


# Raw DIOPT table for filtering mouse/human orthologs. Note that this was a 
# provided data dump on Pavlab servers
diopt <- read.delim("/space/grp/DIOPT/DIOPTvs8_export_Sanja Rogic.txt", stringsAsFactors = FALSE)

# Protein coding genes
pc_hg <- read.delim(ref_hg, stringsAsFactors = FALSE)
pc_mm <- read.delim(ref_mm, stringsAsFactors = FALSE)

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


if (!file.exists(ortho)) {
  write.table(symbols,
              sep = "\t",
              quote = FALSE,
              file = ortho)
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


download_blacklist(bl_hg, bl_url_hg, chr_hg)
download_blacklist(bl_mm, bl_url_mm, chr_mm)



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


download_ccre(ccre_hg, ccre_url_hg, chr_hg)
download_ccre(ccre_mm, ccre_url_mm, chr_mm)



# 6) DE prior ranking (likelihood of being measured as DE across studies)
# Original PNAS version human only) and Nathaniel updated/expanded version
# https://doi.org/10.1073/pnas.1802973116
# ------------------------------------------------------------------------------


# Unpublished at time of writing, but using same algo. as in PNAS paper but more
# and diverse platforms
depr_in_hg <- "/home/nlim/MDE/RScripts/Chapter_4/SPACE_RDATA/Analysis/human/Final/Rare_Prior.RDS"
depr_in_mm <- "/home/nlim/MDE/RScripts/Chapter_4/SPACE_RDATA/Analysis/mouse/Final/Rare_Prior.RDS"


format_depr <- function(infile, outfile) {
  
  if (!file.exists(outfile)) {
    
    depr <- readRDS(infile)
    
    colnames(depr) <- str_replace(colnames(depr), "[\\.]", "_")
    colnames(depr)[colnames(depr) == "gene_Name"] <- "Symbol"
    colnames(depr)[colnames(depr) == "deStrict_mfxRank"] <- "DE_Prior_Rank"
    
    write.table(depr,
                sep = "\t",
                quote = FALSE,
                row.names = FALSE,
                outfile)
  }
}


format_depr(depr_in_hg, depr_out_hg)
format_depr(depr_in_mm, depr_out_mm)



# Human from PNAS paper
# NOTE: Now gives 503 error on download request. Website works fine
depr_url <- "https://www.pnas.org/doi/suppl/10.1073/pnas.1802973116/suppl_file/pnas.1802973116.sd02.txt"


# if (!file.exists(depr_out_pnas)) {
#   
#   download.file(depr_url, destfile = depr_out_pnas)
#   table <- read.delim(depr_out_pnas, stringsAsFactors = FALSE)
#   colnames(table)[colnames(table) == "Gene_Name"] <- "Symbol"
#   
#   write.table(table,
#               sep = "\t",
#               quote = FALSE,
#               row.names = FALSE,
#               outfile_hg_pnas)
# }
