## This script tracks packages used throughout the manuscript
## NOTE: Most packages were individually installed over the course of this
## project, not as a single batch as organized below. User beware if installing
## all from scratch.
## -----------------------------------------------------------------------------


packages <- c(
  "tidyverse",
  "BiocManager",
  "parallel",
  "pheatmap",
  "DescTools",
  "cowplot",
  "RColorBrewer",
  "googlesheets4",
  "RCurl",
  "rjson",
  "scales",
  "Matrix",
  "WGCNA",
  "VennDiagram",
  "egg",
  "viridis",
  "factoextra",
  "ggsci",
  "viridisLite",
  "ggExtra",
  "ggrepel",
  "ROCR",
  "data.table"
)


installed_packages <- packages %in% rownames(installed.packages())


if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}


if (!require("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}


BiocManager::install("biomaRt")
BiocManager::install("preprocessCore")
BiocManager::install("GenomicRanges")
BiocManager::install("limma")
BiocManager::install("edgeR")
BiocManager::install("GeneOverlap")


# Trackplot + bwtools
remotes::install_github(repo = "poisonalien/trackplot")
# installing bwtools: https://gist.github.com/PoisonAlien/e19b482ac6146bfb03142a0de1c4fbc8

