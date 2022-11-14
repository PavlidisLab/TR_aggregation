## This script establishes pathing of data and plots used throughout the project
## TODO: consolidate expression dirs (and parts of scratch)
## TODO: collapsing meta paths/subdirs? (eg meta/chipseq/atlas append)
## TODO: schematic of ChIP-seq workflow
## TODO: input for pipeline install and fastq download better organized
## TODO: ask about local vs utils functions
## TODO: pathing of pipeout chip/


# Variables of interest
# ------------------------------------------------------------------------------



date <- "Apr2022"  # Data freeze to use for analysis
cores <- 8  # For use in parallel
min_peaks <- 100  # Min peak filter for keeping ChIP-seq experiments



# Hard-coded data paths on Pavlab
# ------------------------------------------------------------------------------

# ENCODE pipeline installation
encode_dir <- "/home/amorin/Projects/encode-pipeline/"

# Where downloaded fastq files live
fastq_dir <- "/cosmos/data/downloaded-data/chipseq"

# Where the ENCODE pipeline output lives
pipeout_dir <-  "/cosmos/data/pipeline-output/chipseq-encode-pipeline/"

# Expression platform info from Nathaniel Lim
platform_path <- "/space/grp/nlim/CronGemmaDump/AD_Dump.TSV"

# Scratch location where outputs were variably saved
scratch_dir <- "/home/amorin/scratch/R_objects/"



# Googlesheets IDs
# ------------------------------------------------------------------------------


# Perturbation metadata
gsheets_perturb <- "1oXo1jfoPYcX94bq2F6Bqm1Es1glc8g9mnJvYAO37Vw4"


# ChIP-seq metadata
gsheets_chip <- "1rGVnLL0eXHqr97GM1tloiWwwrJUUmj_ZjW5UOHFN1cc"



# Metadata and other genomic tables
# ------------------------------------------------------------------------------


meta_dir <- "/home/amorin/Data/Metadata/"

# ENCODE blacklists
bl_path_hg <- "/home/amorin/Data/Chromosome_info/blacklist_hg38.tsv"
bl_path_mm <- "/home/amorin/Data/Chromosome_info/blacklist_mm10.tsv"

# ENCODE candidate cis regulatory elements (cCREs)
ccre_path_hg <- "/home/amorin/Data/Chromosome_info/cCREs_V3_hg38.bed"
ccre_path_mm <- "/home/amorin/Data/Chromosome_info/cCREs_V3_mm10.bed"



# ChIP-seq 
# ------------------------------------------------------------------------------


# Gene x experiment bind score matirces
cmat_dir <- "/home/amorin/Data/Annotated_objects/Bind_matrices/Encpipe/"

# ChIP-seq genomic range objects
gr_dir <- "/home/amorin/Data/Annotated_objects/GRanges/"

# Top level of ChIP-seq plot files
cplot_dir <- "/home/amorin/Plots/Chipseq/"


# For trackplot
# installing bwtools: https://gist.github.com/PoisonAlien/e19b482ac6146bfb03142a0de1c4fbc8
# saved this in my own /home/user/bin, then added this dir to R path:

bwtool_path <- "/home/amorin/bin/bwtool/"
if (!grepl(Sys.getenv("PATH"), bwtool_path)) {
  Sys.setenv(PATH = paste(Sys.getenv("PATH"), bwtool_path, sep = ":"))
}


# c("GRanges/", "Trackplots/")


# Perturbation
# ------------------------------------------------------------------------------


# Where the result set files live
rs_dir <- "/home/amorin/Data/Expression_files/Gemma/Resultsets/"

# Where to save the list of processed perturbation resultsets
expr_dir <- "/home/amorin/Data/Expression_files/Gemma/"

# Where to save the perturb effect size matrices
pmat_dir <- "/home/amorin/Data/Expression_files/Perturb_matrix/"

# Top level of perturb plot files
pplot_dir <- "/home/amorin/Plots/TF_perturb/"


# c("Meta_sample_matrix/", "Effect_size/", "Describe_FDR_counts/", "Experiment_similarity/")
