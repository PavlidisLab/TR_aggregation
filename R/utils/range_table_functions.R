## Code to load/work with ChIP-seq narrowPeak/range/BED tables and GRanges

library(GenomicRanges)
library(tidyverse)
library(parallel)
library(rtracklayer)


# Read peak tables
# ------------------------------------------------------------------------------


read_encpeak <- function(dir,
                         peakset = "idr",
                         type = "optimal",
                         path = "/cosmos/data/pipeline-output/chipseq-encode-pipeline/chip/") {
  # Read peak table generated from ENCODE pipeline.
  # Assumes that input dir has the standard dir structure output by the pipeline.
  # peakset corresponds to IDR or overlap reproducible peak sets
  # type can either be optimal or conservative - suggest optimal.
  
  stopifnot(peakset %in% c("idr", "overlap") | type %in% c("optimal", "conservative"))
  
  path <- 
    paste0(dir, "/peak/", peakset, "_reproducibility/", peakset, ".", type, "_peak.narrowPeak.gz")
  
  
  # if (chip == "TF") {
  #   path <- 
  #     paste0(dir, "/peak/idr_reproducibility/idr.", type, "_peak.narrowPeak.gz")
  # } else {
  #   path <- 
  #     paste0(dir, "/peak/overlap_reproducibility/overlap.", type, "_peak.narrowPeak.gz")
  # }
  
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


order_chroms <- function(range_table) {
  # returns the table ordered by chromosome and start
  
  stopifnot(c("Chromosome", "Start") %in% names(range_table))
  
  chroms <- as.character(c(1:22, "MT", "X", "Y"))
  
  order_ix <-
    order(match(str_replace(range_table$Chromosome, "chr", ""), chroms),
          range_table$Start)
  
  range_table <- range_table[order_ix,]
  rownames(range_table) <- NULL
  return(range_table)
}


get_summit_position <- function(range_table) {
  # Adds a column for the peak summit position, and strips the original 
  # MACS column which just shows distance relative to peak start
  
  stopifnot(c("Summit_from_start", "Start") %in% names(range_table))
  
  range_table$Peak_summit <- range_table$Start + range_table$Summit_from_start
  range_table <- range_table[, colnames(range_table) != "Summit_from_start"]
  return(range_table)
}



get_standard_chr <- function(range_table) {
  # remove 'chr' prefix of chromosome identifiers, coerce mitochondrial to 
  # 'MT' and only keep standard chromosomes
  
  stopifnot("Chromosome" %in% names(range_table))
  
  range_table$Chromosome <- str_replace(range_table$Chromosome, "^chr", "")
  range_table$Chromosome <- str_replace(range_table$Chromosome, "^M$", "MT")
  range_table <- filter(range_table, Chromosome %in% c(1:22, "MT", "X", "Y"))
  return(range_table)
}


strand_to_plusminus <- function(range_table) {
  # change strand information from 1/-1 to +/-
  
  stopifnot("Strand" %in% names(range_table))
  
  range_table$Strand <- str_replace(as.character(range_table$Strand), "^1$", "+")
  range_table$Strand <- str_replace(as.character(range_table$Strand), "^-1$", "-")
  return(range_table)
}


startend_to_summit <- function(range_table) {
  # Coerce start and end coordinates to be the peak summit - used to fix a peak
  # as a single point when performing overlaps
  
  stopifnot(c("Start", "End", "Peak_summit") %in% names(range_table))
  
  range_table$Start <- range_table$End <- range_table$Peak_summit
  return(range_table)
}


startend_to_tss <- function(range_table) {
  # Coerce start and end of gene annotation table to just the TSS - used to fix
  # TSS as a single point when performing overlaps
  
  stopifnot(c("Start", "End", "Transcription_start_site") %in% names(range_table))
  
  range_table$Start <- range_table$End <- range_table$Transcription_start_site
  return(range_table)
}


# GRanges convenience functions
# ------------------------------------------------------------------------------


format_peak_attributes <- function(range_table) {
  # Used to format a peak table to be compatible for conversion to a GRanges
  # object for overlap with a gene annotation table. 
  
  range_table <- range_table %>% 
    get_summit_position() %>% 
    get_standard_chr() %>% 
    startend_to_summit() %>% 
    dplyr::select(-Strand)  # not useful for our ChIP-seq - all are "."
  return(range_table) 
}


peak_to_gr <- function(range_table, format = TRUE) {
  # Convert a peak table to a GR object. If format is TRUE, the table will
  # first be processed to be consistent with standard used for overlapping with
  # a gene annotation table. 
  
  if (format) {
    range_table <- format_peak_attributes(range_table)
  }
  
  gr <- makeGRangesFromDataFrame(range_table,
                                 keep.extra.columns = TRUE, 
                                 ignore.strand = TRUE)
  return(gr)
}


pc_to_gr <- function(range_table) {
  # Prepare a gene annotation table for overlap and convert to a GR object
  
  range_table <- range_table %>% 
    strand_to_plusminus() %>% 
    startend_to_tss()
  
  gr <- makeGRangesFromDataFrame(range_table,
                                 keep.extra.columns = TRUE, 
                                 ignore.strand = FALSE)
  
  return(gr)
}


bl_to_gr <- function(range_table) {
  # Wrapper to convert blacklist range table to a GR object
  
  stopifnot(c("Chromosome", "Start", "End") %in% names(range_table))
  
  gr <- makeGRangesFromDataFrame(range_table, 
                                 keep.extra.columns = TRUE, 
                                 ignore.strand = TRUE)
  return(gr)
  
}


gr_to_df <- function(gr) {
  # Wrapper to convert a GRanges object back to a dataframe
  
  df <- GenomicRanges::as.data.frame(gr) %>% 
    dplyr::select(-width) %>%  # remove added width column
    dplyr::rename(Chromosome = seqnames, 
                  Start = start, 
                  End = end, 
                  Strand = strand)
  return(df)
}


liftover_peak <- function(peak_table, chain) {
  # Uses rtracklayer::liftOver, which expects a Granges object. Converts
  # peak_table to a GRanges object and after liftover will remove redundant
  # strand column and convert back to a data frame
  
  # format = FALSE keeps range information consistent with liftover formatting
  # format = TRUE for keeping range information consistent with ensembl table
  peak_gr <- peak_to_gr(peak_table, format = FALSE)
  peak_gr_lo <- unlist(rtracklayer::liftOver(peak_gr, chain))
  peak_gr_lo$Strand <- NULL
  peak_table_lo <- gr_to_df(peak_gr_lo)
  return (peak_table_lo)
}



# Binary annotation/nearest with GRanges
# ------------------------------------------------------------------------------


filter_blacklist <- function(peak_gr, bl_gr) {
  # Given GR objects corresponding to a peak table and a blacklisted table,
  # return peak_gr with any ranges overlapping bl_gr removed
  
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


get_gene_names <- function(range_table) {
  # Given a gene annotation range table, return a vector of length equal to the
  # number of rows in range_table, where each element is a concat of the gene 
  # symbol, chromosome, start, end, and transcript ID. 
  # Used for making unique names for repetitively named symbols/IDs
  
  stopifnot(
    c(
      "Gene_ID",
      "Symbol",
      "Chromosome",
      "Start",
      "End",
      "Transcript_ID"
    ) %in% names(range_table)
  )
  
  # Set blanks to NULL
  range_table$Symbol <- 
    ifelse(range_table$Symbol == "", "NULL", range_table$Symbol)
  
  range_table$Transcript_ID <- 
    ifelse(range_table$Transcript_ID == "", "NULL", range_table$Transcript_ID)
  
  
  gene_names <- paste(range_table$Gene_ID,
                      range_table$Symbol,
                      range_table$Chromosome,
                      range_table$Start,
                      range_table$End,
                      range_table$Transcript_ID,
                      sep = "_")
  
  stopifnot(n_distinct(gene_names) == nrow(range_table))
  return(gene_names)
}


subset_to_symbol <- function(bind_vec, range_table) {
  # Assumes bind_vec is a named vector whose elements correspond to bind scores
  # and range_table is a gene annotation coding. Names of bind_vec are assumed
  # to be of form [Gene_ID]_[Symbol]_[Chromosome]_[Start]_[End]_[Transcript_ID].
  # Returns a named binary vector corresponding to the unique symbols in the 
  # names of bind_vec
  
  stopifnot(all(names(bind_vec) == "") | "Symbol" %in% names(range_table))
  
  symbols <- unique(range_table$Symbol)
  symbols <- symbols[symbols != ""]
  
  symbol_vec <- rep(0, length(symbols))
  names(symbol_vec) <- symbols
  
  sub_vec <- bind_vec[bind_vec == 1]
  
  if (length(sub_vec) > 0) {
    bind_symbols <- unique(str_split(names(sub_vec), "_", simplify = TRUE)[, 2])
    symbol_vec[symbols %in% bind_symbols] <- 1
  }
  return(symbol_vec)
}



distance_anno <- function(peak_table, 
                          gene_table, 
                          bl_table = NULL, 
                          distance,
                          unique_symbols = TRUE) {
  # 
  
  # Create GR objects
  peak_gr <- peak_to_gr(peak_table)
  pc_gr <- pc_to_gr(gene_table)
  
  # Remove blacklisted regions
  if (hasArg(bl_table)) {
    bl_gr <- bl_to_gr(bl_table)
    peak_gr <- filter_blacklist(peak_gr, bl_gr)
  }
  
  # Generate GR object of overlaps within distance. Silence verbose warnings
  hits <- suppressWarnings(  
    findOverlaps(
      query = pc_gr,
      subject = peak_gr,
      ignore.strand = TRUE,
      type = "any",
      select = "all",
      maxgap = distance
    )
  )
  
  # Create a gene vector from gene table and fill if hits were found
  gene_vec <- ifelse(1:nrow(gene_table) %in% hits@from, 1, 0)
  names(gene_vec) <- get_gene_names(gene_table) 
  
  if (unique_symbols) {
    gene_vec <- subset_to_symbol(gene_vec, gene_table)
  }
  return(gene_vec)
}


# get_nearest_symbol <- function(peak_gr, pc_gr) {
#   # Given a GR of a chip/narrowpeak table and a protein coding GR object, return
#   # a vector size of length(peak_gr), where each element is the symbol of the 
#   # nearest gene
#   stopifnot(class(peak_gr) == "Granges" | class(pc_gr) == "GRanges")
#   
#   nearest <- pc_gr@elementMetadata$Symbol[nearest(peak_gr, pc_gr)]
#   return(nearest)
# }


# distance_to_nearest <- function(gr1, gr2) {
#   # gr1 and gr2 as Grange objects. for each range in gr1, find the nearest
#   # range in gr2. returns a vector of integers equal to length of gr1/
#   
#   stopifnot(class(gr1) == "Granges" | class(gr2) == "GRanges")
#   
#   distance <- distanceToNearest(gr1, gr2)@elementMetadata$distance
#   return(distance)
# }


# Code to generate the gene assignment scores proposed in Ouyang 2009. 
# Adapted from TFTargetCaller's implementation
# https://www.pnas.org/content/106/51/21521
# https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1003342
# ------------------------------------------------------------------------------


get_dist_mat <- function(peak_table, gene_table, cores = 8) {
  
  # Take in a peak table (typically narrowPeak from MACS2) and a protein coding
  # annotation table. Return a list of equal length to their chromosomes. Each
  # element of the list is a matrix, with rows as gene symbols and columns
  # as unique peaks. Each element of the matrix is the distance between the TSS
  # and the summit of the peak
  
  chromosomes <- intersect(peak_table$Chromosome, gene_table$Chromosome)
  stopifnot(all(unique(gene_table$Strand) %in% c(1, -1)))
  
  dist_mat_per_chr <- mclapply(chromosomes, function(chr) {
    
    gene_chr <- gene_table[which(gene_table$Chromosome == chr),]
    peak_chr <- peak_table[which(peak_table$Chromosome == chr),]
    
    distance <- matrix(
      sapply(peak_chr$Peak_summit, function(x) {
        x - gene_chr$Transcription_start_site
      }), nrow = nrow(gene_chr), ncol = nrow(peak_chr))
    
    distance <- distance * gene_chr$Strand
    
    rownames(distance) <- get_gene_names(gene_chr)
    colnames(distance) <- peak_chr$ID
    
    return(distance)
  }, mc.cores = cores)
  names(dist_mat_per_chr) <- paste0("chr", chromosomes)
  return(dist_mat_per_chr)
}


match_score_symbols <- function(bind_scores, gene_table) {
  # Note that this assumes that scores is named with gene IDs. Fills
  # out the array of unique gene IDs with 0s for any genes that were not
  # represented on the peak table
  
  scores_vec <- rep(0, nrow(gene_table))
  names(scores_vec) <- get_gene_names(gene_table)
  where <- match(names(bind_scores), names(scores_vec))
  scores_vec[where] <- bind_scores
  return(scores_vec)
}


get_ouyang_score <- function(intensity, distance, decay_constant = 5000) {
  # Exponential decay relationship proposed in Ouyang et al., 2009
  # Computes the score for a single peak/gene, where intensity is the
  # MACS2 score for a peak and distance is the peak summit to the TSS of interest
  return (intensity * exp(-(abs(distance)/decay_constant)))
}


get_beta_score <- function(distance, base = 1e5) {
  # Exponential decay relationship proposed in XX
  # Computes the score for a single peak/gene
  return (exp(-0.5 - (4 * abs(distance)/base)))
}


all_ouyang_scores <- function(range_table, 
                              gene_table, 
                              decay_constant = 5e3, 
                              cores = 8,
                              use_intensity = TRUE,
                              dist_threshold = 1e6) {
  
  # Generates individual Ouyang scores for each peak/TSS combo, and then 
  # sums up all peak scores for a gene, returning a vector of scores of length
  # equal to the rows of gene_table
  
  chromosomes <- intersect(range_table$Chromosome, gene_table$Chromosome)
  
  if (use_intensity) {
    peak_intensity_list <- lapply(chromosomes, function(chr) {
      return(range_table[range_table$Chromosome == chr, "Score"])
    })
  } else {
    peak_intensity_list <- lapply(chromosomes, function(chr) {
      return(rep(1, nrow(range_table[range_table$Chromosome == chr, ])))
    })
  }
  
  dist_matrix_list <- get_dist_mat(range_table, gene_table)
  
  ouyang_scores <- unlist(mcmapply(function(a, b) { 
    
    a[abs(a) > dist_threshold] <- NA
    
    result <- apply(a, 1, function(x) {
      y <- get_ouyang_score(b, x, decay_constant)
      return(sum(y, na.rm = TRUE))
    })
    
    return(result) 
  },
  a = dist_matrix_list,
  b = peak_intensity_list,
  mc.cores = cores), 
  use.names = FALSE)
  
  names(ouyang_scores) <- unlist(lapply(dist_matrix_list, row.names))
  ouyang_scores <- match_score_symbols(ouyang_scores, gene_table) 
  return(ouyang_scores)
}


all_beta_scores <- function(range_table,
                            gene_table,
                            cores = 8,
                            dist_base = 1e5) {
  
  # Generates individual Beta scores for each peak/TSS combo, and then 
  # sums up all peak scores for a gene, returning a vector of scores of length
  # equal to the rows of gene_table
  
  chromosomes <- intersect(range_table$Chromosome, gene_table$Chromosome)
  
  dist_matrix_list <- get_dist_mat(range_table, gene_table)

  beta_scores <- unlist(mclapply(dist_matrix_list, function(a) { 
    
    a[abs(a) > dist_base] <- NA

    result <- apply(a, 1, function(x) {
      y <- get_beta_score(x, dist_base)
      return(sum(y, na.rm = TRUE))
    })
    
    return(result) 
  }, mc.cores = cores))
  
  names(beta_scores) <- unlist(lapply(dist_matrix_list, row.names))
  beta_scores <- match_score_symbols(beta_scores, gene_table)
  return(beta_scores)
}



max_score_per_symbol <- function(score, gene_table) {
  # Only take the highest score per symbol.
  
  gene_df <- data.frame(
    Gene = str_split(names(score), pattern = "_", simplify = TRUE)[, 2],
    Score = score
  )
  
  max_gene_df <- gene_df %>% 
    group_by(Gene) %>% 
    arrange(desc(Score), .by_group = TRUE) %>% 
    filter(row_number() == 1)
  
  scores <- max_gene_df$Score
  names(scores) <- max_gene_df$Gene
  scores[order(match(names(scores), gene_table$Symbol))]
}
