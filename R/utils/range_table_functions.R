## Code to load/work with ChIP-seq narrowPeak/range/BED tables as data frames
## and as GenomicRanges objects


library(GenomicRanges)
library(tidyverse)
library(parallel)


# Read peak tables
# ------------------------------------------------------------------------------


# Read peak table generated from ENCODE pipeline.
# Assumes that dir has the default structure output by the pipeline.
# peakset corresponds to IDR or overlap reproducible peak sets
# type can either be optimal or conservative - suggest optimal.

read_encpeak <- function(dir,
                         peakset = "idr",
                         type = "optimal",
                         path = "/cosmos/data/pipeline-output/chipseq-encode-pipeline/chip/") {
  
  
  stopifnot(peakset %in% c("idr", "overlap") | type %in% c("optimal", "conservative"))
  
  path <- 
    paste0(dir, "/peak/", peakset, "_reproducibility/", peakset, ".", type, "_peak.narrowPeak.gz")
  
  
  range_table <- read.delim(
    path,
    header = FALSE,
    stringsAsFactors = FALSE,
    row.names = NULL,
    col.names = c(
      "Chromosome",
      "Start",
      "End",
      "ID",
      "Score",
      "Strand",
      "Foldchange",
      "Neglog10pval",
      "Neglog10qval",
      "Summit_from_start"
    )
  )
  
  return(order_chroms(range_table))
  
}


# Manipulation of range tables
# ------------------------------------------------------------------------------


# returns the table ordered by chromosome and start

order_chroms <- function(range_table) {
  
  stopifnot(c("Chromosome", "Start") %in% names(range_table))
  
  chroms <- as.character(c(1:22, "MT", "X", "Y"))
  
  order_ix <-
    order(match(str_replace(range_table$Chromosome, "chr", ""), chroms),
          range_table$Start)
  
  range_table <- range_table[order_ix,]
  rownames(range_table) <- NULL
  
  return(range_table)
}


# Adds a column for the peak summit position, and strips the original 
# MACS column which just shows distance relative to peak start

get_summit_position <- function(range_table) {
  
  stopifnot(c("Summit_from_start", "Start") %in% names(range_table))
  
  range_table$Peak_summit <- range_table$Start + range_table$Summit_from_start
  range_table <- range_table[, colnames(range_table) != "Summit_from_start"]
  
  return(range_table)
}


# remove 'chr' prefix of chromosome identifiers, coerce mitochondrial to 
# 'MT' and only keep standard autosomal and sex chromosomes

get_standard_chr <- function(range_table) {
  
  stopifnot("Chromosome" %in% names(range_table))
  
  range_table$Chromosome <- str_replace(range_table$Chromosome, "^chr", "")
  range_table$Chromosome <- str_replace(range_table$Chromosome, "^M$", "MT")
  range_table <- filter(range_table, Chromosome %in% c(1:22, "MT", "X", "Y"))
  return(range_table)
}


# change strand information from 1/-1 to +/-

strand_to_plusminus <- function(range_table) {
  
  stopifnot("Strand" %in% names(range_table))
  
  range_table$Strand <- str_replace(as.character(range_table$Strand), "^1$", "+")
  range_table$Strand <- str_replace(as.character(range_table$Strand), "^-1$", "-")
  return(range_table)
}


# Coerce start and end coordinates to be the peak summit - used to fix a peak
# as a single point when performing overlaps

startend_to_summit <- function(range_table) {
  
  stopifnot(c("Start", "End", "Peak_summit") %in% names(range_table))
  
  range_table$Start <- range_table$End <- range_table$Peak_summit
  
  return(range_table)
}


# Coerce start and end of gene annotation table to just the TSS - used to fix
# TSS in protein coding tables as a single point when performing overlaps

startend_to_tss <- function(range_table) {
  
  stopifnot(c("Start", "End", "Transcription_start_site") %in% names(range_table))
  
  range_table$Start <- range_table$End <- range_table$Transcription_start_site
  
  return(range_table)
}


# GRanges convenience functions
# ------------------------------------------------------------------------------


# Used to format a peak table to be compatible for conversion to a GRanges
# object for overlap with a gene annotation table. 

format_peak_attributes <- function(range_table) {
  
  range_table <- range_table %>% 
    get_summit_position() %>% 
    get_standard_chr() %>% 
    startend_to_summit() %>% 
    dplyr::select(-Strand)  # not useful for our ChIP-seq - all are "."
  return(range_table) 
}


# Convert a peak table to a GR object. If format is TRUE, the table will
# first be processed to be consistent with standard used for overlapping with
# a gene annotation table. 

peak_to_gr <- function(range_table, format = TRUE) {
  
  if (format) {
    range_table <- format_peak_attributes(range_table)
  }
  
  gr <- makeGRangesFromDataFrame(range_table,
                                 keep.extra.columns = TRUE, 
                                 ignore.strand = TRUE)
  return(gr)
}



# Convert a gene annotation table to a GR object. If TSS is TRUE, will fix
# gene start and end coordinates to be the 1bp TSS


pc_to_gr <- function(range_table, TSS = TRUE) {
  
  range_table <- strand_to_plusminus(range_table)
  
  if (TSS) {
    range_table <- startend_to_tss(range_table)
  }

  gr <- makeGRangesFromDataFrame(range_table,
                                 keep.extra.columns = TRUE,
                                 ignore.strand = FALSE)
  
  return(gr)
}


# Convert a blacklisted regions table to a GR object

bl_to_gr <- function(range_table) {
  
  stopifnot(c("Chromosome", "Start", "End") %in% names(range_table))
  
  gr <- makeGRangesFromDataFrame(range_table, 
                                 keep.extra.columns = TRUE, 
                                 ignore.strand = TRUE)
  return(gr)
  
}


# Convert a GRanges object back to a dataframe

gr_to_df <- function(gr) {
  
  df <- GenomicRanges::as.data.frame(gr) %>% 
    dplyr::rename(Chromosome = seqnames, 
                  Start = start, 
                  End = end, 
                  Strand = strand)
  return(df)
}



# Fix range of GR to the summit position and pad with a fixed window

summit_window <- function(gr, window_size) {
  
  stopifnot("Summit_from_start" %in% names(gr@elementMetadata))
  
  start(gr) <- end(gr) <- start(gr) + gr$Summit_from_start
  gr <- gr + window_size
  return(gr)
}


# Given GR objects corresponding to a peak table and a blacklist table,
# return peak_gr with any ranges overlapping bl_gr removed

filter_blacklist <- function(peak_gr, bl_gr) {
  
  hits <- suppressWarnings(
    findOverlaps(
      query = peak_gr,
      subject = bl_gr,
      ignore.strand = TRUE,
      type = "any",
      select = "all"
    )
  )
  if (length(hits) > 0) {
    peak_gr <- peak_gr[-hits@from]
  }
  
  return(peak_gr)
}


# Binding score gene assignment using GRanges
# ------------------------------------------------------------------------------


# Get a binary vector indicating whether genes had a proximal binding event
#
# pc_gr: A GR object of protein coding gene with range fixed to the 1bp TSS
# peak_gr: A GR object of the peaks with range fixed to the 1bp summit
# max_dist: An integer of the distance cutoff in bps
# returns a numeric vector the length of pc_gr of the binary binding status


binary_scores <- function(pc_gr, peak_gr, max_dist = 25e3) {
  
  stopifnot(class(pc_gr) == "GRanges", class(peak_gr) == "GRanges")
  
  hits <- suppressWarnings(
    findOverlaps(
      query = pc_gr,
      subject = peak_gr,
      ignore.strand = TRUE,
      type = "any",
      select = "all",
      maxgap = max_dist
    )
  )
  
  gene_vec <- ifelse(1:length(pc_gr) %in% hits@from, 1, 0)
  names(gene_vec) <- pc_gr$Symbol
  
  return(gene_vec)
}


# Generate a gene binding score of peak-gene distance using an exponential decay
# function proposed in Ouyang et al., 2009 https://www.pnas.org/content/106/51/21521
# The original formulation scaled the score by the MACS2 score, which is 
# excluded here following contemporary practices
#
# distance: A vector of integers of the basepairs between a gene TSS and peak summits.
# decay_constant: An integer controlling how steeply the score decreases
# returns: An integer


ouyang <- function(distance, decay_constant = 5e3) {
  
  stopifnot(is.numeric(distance), length(distance) > 0)
  
  scores <- lapply(distance, function(x) {
    exp(-(abs(x) / decay_constant))
  })
  
  return(sum(unlist(scores)))
}



# Generate a gene binding score of peak-gene distance using formulation proposed
# in Wang et al., 2013 https://www.nature.com/articles/nprot.2013.150 
#
# distance: A vector of integers of the basepairs between a gene TSS and peak summits.
# decay_constant: An integer controlling how steeply the score decreases
# returns: An integer


beta <- function(distance, base = 1e5) {
  
  stopifnot(is.numeric(distance), length(distance) > 0)
  
  scores <- lapply(distance, function(x) {
    exp(-(0.5 + 4 * abs(x) / base))
  })
  
  return(sum(unlist(scores)))
}



# Generate binding scores for every gene from a provided peak GR object.
#
# pc_gr: A GR object of protein coding gene with range fixed to the 1bp TSS
# peak_gr: A GR object of the peaks with range fixed to the 1bp summit
# method: A character indicating which scoring method to use
# max_dist: An integer specifying the cutoff (in bps) of peaks to consider
# ncore: An integer of how many cores parallel will use
# ...: Additional arguments supplied to the scoring methods
# returns a numeric vector the length of pc_gr of the gene binding scores


binding_scores <- function(pc_gr, 
                           peak_gr, 
                           method, 
                           max_dist = 1e6, 
                           ncore = 1,
                           ...) {
  
  stopifnot(class(pc_gr) == "GRanges", class(peak_gr) == "GRanges")
  stopifnot(method %in% c("Ouyang", "Beta"))
  
  score_l <- mclapply(1:length(pc_gr), function(x) {
    
    dist <- GenomicRanges::distance(pc_gr[x], peak_gr, select = "all")
    dist <- dist[!is.na(dist) & dist < max_dist]
    
    if (length(dist) == 0) {
      return(0)
    }
    
    if (method == "Ouyang") {
      score <- ouyang(dist, ...)
    } else {
      score <- beta(dist, ...)
    }
  }, mc.cores = ncore)
  
  names(score_l) <- pc_gr$Symbol
  
  return(unlist(score_l))
}
