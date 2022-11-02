## This script establishes pathing of data and plots used throughout the project



date <- "Apr2022"  # Data freeze to use for analysis
cores <- 8  # For use in parallel


# Hard-coded data paths: big data on Pavlab
# ------------------------------------------------------------------------------



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

# Where to save the processed objects
expr_dir <- "/home/amorin/Data/Expression_files/Gemma/"