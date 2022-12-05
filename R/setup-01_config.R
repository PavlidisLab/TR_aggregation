## This script establishes pathing of data and plots used throughout the project
## TODO: consolidate expression dirs (and parts of scratch)
## TODO: collapsing meta paths/subdirs? (eg meta/chipseq/atlas append)
## TODO: schematic of ChIP-seq workflow
## TODO: input for pipeline install and fastq download better organized
## TODO: ask about local vs utils functions (eg, int01 ranking function)
## TODO: pathing of pipeout chip/


# Variables of interest
# ------------------------------------------------------------------------------

date <- "Apr2022"  # Data freeze to use for analysis
cores <- 8  # For use in parallel

min_peaks <- 100  # Min peak filter for keeping ChIP-seq experiments
binary_dist <- 25e3  # distance threshold used for binary gene scores

fdr <- 0.1  # false discovery rate for perturbation experiments to be diff expr. 


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


# c("GRanges/", "Trackplots/", "Preprocess/", "Binding_similarity/", "Binding_summary/")


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


# Intersect
# ------------------------------------------------------------------------------


# Top level of plot files
iplot_dir <- "/home/amorin/Plots/Intersect/"

# c("Experiment_similarity/", "Gene_rankings")



# Figure plot paths
# ------------------------------------------------------------------------------


# F1 B
# "/home/amorin/Plots/TF_perturb/Meta_sample_matrix/batch1_tf_perturb_counts_Apr2022.png"
# "/home/amorin/Plots/Chipseq/Describe_meta/batch1_chip_experiment_counts_all_Apr2022.png"

# FS1 A, B, C
# "/home/amorin/Plots/Chipseq/Describe_meta/batch1_chip_count_tf-idr_mecp2_overlap_peaks_by_symbol_Apr2022.png"
# "/home/amorin/Plots/Chipseq/Describe_meta/batch1_chip_count_tf-idr_mecp2_overlap_peaks_by_input_Apr2022.png"
# "/home/amorin/Plots/Chipseq/Describe_meta/batch1_chip_count_tf-idr_mecp2_overlap_peaks_vs_avg_exp_mreadsApr2022.png"

# F2 A, B, C, D, E
# "/home/amorin/Plots/Chipseq/Binding_similarity/Densplot_all_human_Apr2022.png"
# "/home/amorin/Plots/Chipseq/Binding_similarity/Densplot_all_mouse_Apr2022.png"
# "/home/amorin/Plots/Chipseq/Binding_similarity/Densplot_all_ortho_Apr2022.png"
# "/home/amorin/Plots/Chipseq/Binding_summary/Bind_specificity_heatmap_human_Apr2022.png"
# "/home/amorin/Plots/Chipseq/Binding_summary/VBplot_bindscore_human_Apr2022.png"

# FS2 A, B, C
# "/home/amorin/Plots/Chipseq/Binding_similarity/Vbplot_tf_human_intersect_Apr2022.png"
# "/home/amorin/Plots/Chipseq/Binding_similarity/Vbplot_tf_mouse_intersect_Apr2022.png"
# "/home/amorin/Plots/Chipseq/Binding_similarity/Vbplot_tf_ortho_intersect_Apr2022.png"

# FS3 A, B
# "/home/amorin/Plots/Chipseq/Binding_similarity/Densplot_cor_runx1_ortho_Apr2022.png"
# "/home/amorin/Plots/Chipseq/Binding_similarity/Densplot_inter_runx1_ortho_Apr2022.png"

# FS4 A, B
# "/home/amorin/Plots/Chipseq/Binding_summary/Bind_specificity_heatmap_mouse_Apr2022.png"
# "/home/amorin/Plots/Chipseq/Binding_summary/VBplot_bindscore_mouse_Apr2022.png"

# FS5 A, B
# "/home/amorin/Plots/Chipseq/Binding_summary/Human_top_bound.png"
# "/home/amorin/Plots/Chipseq/Binding_summary/Mouse_top_bound.png"

# F3 A, B (+/- legend, as ended up labeling axis in Adobe due to clipping)
# "/home/amorin/Plots/Chipseq/GRanges/cCRE_overlap_by_exp_human_legend_Apr2022.png"
# "/home/amorin/Plots/Chipseq/GRanges/cCRE_overlap_by_exp_human_nolegend_Apr2022.png"
# "/home/amorin/Plots/Chipseq/GRanges/cCRE_allprop_human_Apr2022.png"
# "/home/amorin/Plots/Chipseq/GRanges/cCRE_overlap_by_exp_mouse_legend_Apr2022.png"
# "/home/amorin/Plots/Chipseq/GRanges/cCRE_overlap_by_exp_mouse_nolegend_Apr2022.png"
# "/home/amorin/Plots/Chipseq/GRanges/cCRE_allprop_mouse_Apr2022.png"

# F6 A, B (+/- overlaid text, as ended labeling in Adobe as trackplot text too small)
# "/home/amorin/Plots/Chipseq/Trackplots/Human_ASCL1_titles_chr9:38022988-38023319.png"
# "/home/amorin/Plots/Chipseq/Trackplots/Human_ASCL1_notitles_chr9:38022988-38023319.png"
# "/home/amorin/Plots/Chipseq/Trackplots/Human_sample-TR_titles_chr9:38022988-38023319.png"
# "/home/amorin/Plots/Chipseq/Trackplots/Human_sample-TR_notitles_chr9:38022988-38023319.png"
# "/home/amorin/Plots/Chipseq/Trackplots/ASCL1_SHB_trackplot_experimentIDs.tsv"

# FS7 A, B, C, D, E, F
# "/home/amorin/Plots/TF_perturb/Meta_sample_matrix/Gene_measure_count_human_Apr2022.png"
# "/home/amorin/Plots/TF_perturb/Meta_sample_matrix/Gene_measure_count_mouse_Apr2022.png"
# "/home/amorin/Plots/TF_perturb/Meta_sample_matrix/human_sample_matrix_clustered_Apr2022.png"
# "/home/amorin/Plots/TF_perturb/Meta_sample_matrix/mouse_sample_matrix_clustered_Apr2022.png"
# "/home/amorin/Plots/TF_perturb/Meta_sample_matrix/batch1_tf_perturb_tech_counts_Apr2022.png"
# "/home/amorin/Plots/TF_perturb/Meta_sample_matrix/batch1_tf_perturb_counts_Apr2022.png"

# F4 A, B (+/- legend to better fit main plot in figure), C 
# "/home/amorin/Plots/TF_perturb/Effect_size/DEG_counts_boxplot_Apr2022.png"
# "/home/amorin/Plots/TF_perturb/Effect_size/DEG_counts_byPertandSpecies_Apr2022.png"
# "/home/amorin/Plots/TF_perturb/Effect_size/DEG_counts_byPertandSpecies_nolegend_Apr2022.png"
# "/home/amorin/Plots/TF_perturb/Experiment_similarity/cor_heatmap_ortho_Apr2022.png"

# FS8 A, B, C, D, E
# "/home/amorin/Plots/TF_perturb/Experiment_similarity/Densplot_pval-intersect_all_human_Apr2022.png"
# "/home/amorin/Plots/TF_perturb/Experiment_similarity/Densplot_pval-intersect_all_mouse_Apr2022.png"
# "/home/amorin/Plots/TF_perturb/Experiment_similarity/Densplot_pval-intersect_all_ortho_Apr2022.png"
# "/home/amorin/Plots/TF_perturb/Experiment_similarity/Vbplot_all_human_Apr2022.png"
# "/home/amorin/Plots/TF_perturb/Experiment_similarity/Vbplot_all_mouse_Apr2022.png"

# FS9 A, B, C
# "/home/amorin/Plots/TF_perturb/Effect_size/Perturbed_TF_FC_Apr2022.png"
# "/home/amorin/Plots/TF_perturb/Effect_size/DEG_counts_vs_FC_Apr2022.png"
# "/home/amorin/Plots/TF_perturb/Effect_size/DEG_counts_vs_PercRankFC_Apr2022.png"

# FS10 A, B, C
# "/home/amorin/Plots/TF_perturb/Effect_size/DEG_counts_boxplot_byperturb_allexp_Apr2022.png"
# "/home/amorin/Plots/TF_perturb/Effect_size/DEG_counts_boxplot_byperturb_exclnodeg_Apr2022.png"
# "/home/amorin/Plots/TF_perturb/Effect_size/hist_countde_propde_FDR=0.1_Apr2022.png"

# F5 A, B (both +/- legend to fit better in Adobe), C, D, E
# "/home/amorin/Plots/TF_perturb/Describe_FDR_counts/prior_decount_bin_human_FDR=0.1_Apr2022.png"
# "/home/amorin/Plots/TF_perturb/Describe_FDR_counts/prior_decount_bin_human_nolegend_FDR=0.1_Apr2022.png"
# "/home/amorin/Plots/TF_perturb/Describe_FDR_counts/prior_decount_bin_mouse_FDR=0.1_Apr2022.png"
# "/home/amorin/Plots/TF_perturb/Describe_FDR_counts/prior_decount_bin_mouse_nolegend_FDR=0.1_Apr2022.png"
# "/home/amorin/Plots/TF_perturb/Describe_FDR_counts/TF_hist_decounts_FDR=0.1_Apr2022.png"
# "/home/amorin/Plots/TF_perturb/Describe_FDR_counts/Human_ASCL1_hist_decounts_FDR=0.1_Apr2022.png"
# "/home/amorin/Plots/TF_perturb/Describe_FDR_counts/Mouse_Ascl1_hist_decounts_FDR=0.1_Apr2022.png"
# "/home/amorin/Plots/TF_perturb/Describe_FDR_counts/Human_ASCL1_FC_heatmap_Apr2022.png"
# "/home/amorin/Plots/TF_perturb/Describe_FDR_counts/Human_ASCL1_DEprior_binary_heatmap_Apr2022.png"
# "/home/amorin/Plots/TF_perturb/Describe_FDR_counts/Mouse_Ascl1_FC_heatmap_Apr2022.png"
# "/home/amorin/Plots/TF_perturb/Describe_FDR_counts/Mouse_Ascl1_DEprior_binary_heatmap_Apr2022.png"

# FS11
# "/home/amorin/Plots/TF_perturb/Describe_FDR_counts/human_max_countde_vs_deprior_Apr2022.png"

# FS12 A, B
# "/home/amorin/Plots/TF_perturb/Describe_FDR_counts/Purity_vs_decount_bplot_mincount_FDR=0.1_Apr2022.png"
# "/home/amorin/Plots/TF_perturb/Describe_FDR_counts/Human_ASCL1_example_purity_heatmap_Apr2022.png"

# FS13 A, B, C, D, E
#
#
#
#
#

# F6 B, C, D, E, F, G
#
#
#
#
#

# FS14

# FS15

# FS16 A, B, C
#
#
#


