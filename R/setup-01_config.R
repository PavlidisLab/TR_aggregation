## This script establishes pathing of data and plots used throughout the project



date <- "Apr2022"  # Data freeze to use for analysis
cores <- 8  # For use in parallel


# Hard-coded data paths on Pavlab
# ------------------------------------------------------------------------------


# Expression platform info from Nathaniel Lim
platform_path <- "/space/grp/nlim/CronGemmaDump/AD_Dump.TSV"




# Googlesheets IDs
# ------------------------------------------------------------------------------


# Perturbation metadata
gsheets_perturb <- "1oXo1jfoPYcX94bq2F6Bqm1Es1glc8g9mnJvYAO37Vw4"


# ChIP-seq metadata
# gsheets_chip


# Metadata and other genomic tables
# ------------------------------------------------------------------------------


meta_dir <- "/home/amorin/Data/Metadata/"





# ChIP-seq 
# ------------------------------------------------------------------------------





# Perturbation
# ------------------------------------------------------------------------------


# Where the result set files live
rs_dir <- "/home/amorin/Data/Expression_files/Gemma/Resultsets/"

# Where to save the list of processed perturbation resultsets
expr_dir <- "/home/amorin/Data/Expression_files/Gemma/"

# Where to save the perturb effect size matrices
pmat_dir <- "/home/amorin/Data/Expression_files/Perturb_matrix/"

# Top level of perturb files
pplot_dir <- "/home/amorin/Plots/TF_perturb/"
