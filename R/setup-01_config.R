## This script establishes pathing of data and plots used throughout the project
## TODO: consolidate expression dirs (and parts of scratch)
## TODO: collapsing meta paths/subdirs? (eg meta/chipseq/atlas append)
## TODO: schematic of ChIP-seq workflow
## TODO: input for pipeline install and fastq download better organized
## TODO: ask about local vs utils functions
## TODO: pathing of pipeout chip/
## TODO: bulky plot functions local or utils?
## TODO: ask about bulky config!
## TODO: filepath instead of paste0

# Variables of interest
# ------------------------------------------------------------------------------

date <- "Apr2022"  # Data freeze to use for analysis
cores <- 8  # For use in parallel

min_peaks <- 100  # Min peak filter for keeping ChIP-seq experiments
binary_dist <- 25e3  # Absolute distance threshold used for binary gene scores

fdr <- 0.1  # False discovery rate for perturbation experiments to be diff expr. 


# TODO: better data pathing
# Scratch location where outputs were variably saved
scratch_dir <- "/space/scratch/amorin/R_objects/"



# Hard-coded data paths on Pavlab
# ------------------------------------------------------------------------------


# ENCODE pipeline installation
encode_dir <- "/space/grp/amorin/Projects/encode-pipeline/"

# Where downloaded fastq files live
fastq_dir <- "/cosmos/data/downloaded-data/chipseq"

# Where the ENCODE pipeline output lives
pipeout_dir <-  "/cosmos/data/pipeline-output/chipseq-encode-pipeline/"

# Expression platform info from Nathaniel Lim
platform_path <- "/space/grp/nlim/CronGemmaDump/AD_Dump.TSV"

# DE prior from Nathaniel Lim
depr_in_hg <- "/home/nlim/MDE/RScripts/Chapter_4/SPACE_RDATA/Analysis/human/Final/Rare_Prior.RDS"
depr_in_mm <- "/home/nlim/MDE/RScripts/Chapter_4/SPACE_RDATA/Analysis/mouse/Final/Rare_Prior.RDS"

# DIOPT data dump for orthologous genes provided by Sanja Rogic
diopt_path <- "/space/grp/DIOPT/DIOPTvs8_export_Sanja Rogic.txt"


# Googlesheets IDs for curated sheets
# ------------------------------------------------------------------------------


# Perturbation metadata
gsheets_perturb <- "1oXo1jfoPYcX94bq2F6Bqm1Es1glc8g9mnJvYAO37Vw4"

# ChIP-seq metadata
gsheets_chip <- "1rGVnLL0eXHqr97GM1tloiWwwrJUUmj_ZjW5UOHFN1cc"

# Curated targets - note that this is a copy of the master sheet taken on July 
# 4th 2022 as a data freeze
gsheets_curated <- "1ngjKoRGaOgF-8BlxUPK7o7XRg7wimTxYYQkSokSVYUM"


# Metadata and other genomic tables
# ------------------------------------------------------------------------------


meta_dir <- "/space/grp/amorin/Metadata/"

# Final curated perturbation experiment metadata
perturb_meta_path <- paste0(meta_dir, "batch1_tfperturb_meta_final_", date, ".tsv")

# ENCODE blacklists
bl_path_hg <- "/space/grp/amorin/Chromosome_info/blacklist_hg38.tsv"
bl_path_mm <- "/space/grp/amorin/Chromosome_info/blacklist_mm10.tsv"

# ENCODE candidate cis regulatory elements (cCREs)
ccre_path_hg <- "/space/grp/amorin/Chromosome_info/cCREs_V3_hg38.bed"
ccre_path_mm <- "/space/grp/amorin/Chromosome_info/cCREs_V3_mm10.bed"

# Ensembl and refseq select protein coding genes
ref_path_hg <- paste0(meta_dir, "refseq_select_hg38.tsv")
ref_path_mm <- paste0(meta_dir, "refseq_select_mm10.tsv")
ens_path_hg <- paste0(meta_dir, "ensembl_human_protein_coding_V98.tsv")
ens_path_mm <- paste0(meta_dir, "ensembl_mouse_protein_coding_V98.tsv")

# Tables of transcription factors
tf_path_hg <- paste0(meta_dir, "human_tfs.tsv")
tf_path_mm <- paste0(meta_dir, "mouse_tfs.tsv")

# Mapping of 1:1 high-confidence orthologous genes
ortho_path <- paste0(meta_dir, "hg_mm_1to1_ortho_genes_DIOPT-v8.tsv")

# Output of DE prior rankings
depr_path_hg <- paste0(meta_dir, "DE_prior_hg.tsv")
depr_path_mm <- paste0(meta_dir, "DE_prior_mm.tsv")

# Raw input of Chu2021 curated targes
chu2021_path_records <- paste0(meta_dir, "Chu2021_records_DTRI.tsv")
chu2021_path_all <- paste0(meta_dir, "Chu2021_all_DTRI.tsv")

# Output of formatted curated targets
curated_path_all <- paste0(meta_dir, "Curated_targets_all_July2022.tsv")
curated_path_pavlab <- paste0(meta_dir, "Curated_targets_pavlab_July2022.tsv")



# ChIP-seq 
# ------------------------------------------------------------------------------


# Gene x experiment bind score matirces
cmat_dir <- "/space/grp/amorin/Annotated_objects/Bind_matrices/Encpipe/"

# ChIP-seq genomic range objects
gr_dir <- "/space/grp/amorin/Annotated_objects/GRanges/"


# For trackplot
# installing bwtools: https://gist.github.com/PoisonAlien/e19b482ac6146bfb03142a0de1c4fbc8
# saved this in my own /home/user/bin, then added this dir to R path:

bwtool_path <- "/space/grp/amorin/bin/bwtool/"
if (!grepl(Sys.getenv("PATH"), bwtool_path)) {
  Sys.setenv(PATH = paste(Sys.getenv("PATH"), bwtool_path, sep = ":"))
}


# Perturbation
# TODO: reconsider this pathing structure
# ------------------------------------------------------------------------------


# Where the result set files live
rs_dir <- "/space/grp/amorin/Expression_files/Gemma/Resultsets/"

# Where to save the list of processed perturbation resultsets
expr_dir <- "/space/grp/amorin/Expression_files/Gemma/"

# Where to save the perturb effect size matrices
pmat_dir <- "/space/grp/amorin/Expression_files/Perturb_matrix/"

# RDS list objects of the filtered and unfiltered perturbation result sets
rs_filt_path <- paste0(expr_dir, "TF_perturb_batch1_rslist_", date, ".RDS")
rs_unfilt_path <- paste0(expr_dir, "TF_perturb_batch1_unfiltered_rslist_", date, ".RDS")

# List of matrices of perturbation effect sizes
pmat_path_hg <- paste0(pmat_dir, "human_list_perturb_matrix_", date, ".RDS")
pmat_path_mm <- paste0(pmat_dir, "mouse_list_perturb_matrix_", date, ".RDS")
pmat_path_ortho <- paste0(pmat_dir, "ortho_list_perturb_matrix_", date, ".RDS")

# Output list of DE counts, grouped by TR or for all experiment combined
prank_path_group <- paste0(expr_dir, "TF_perturb_DE_counts_list_by_TF_FDR01_", date, ".RDS")
prank_path_all <- paste0(expr_dir, "TF_perturb_DE_counts_list_all_FDR01_", date, ".RDS")

# R object of similarity between perturbation experiments
perturb_sim_path <- paste0(scratch_dir, "perturb_similarity_", date, ".RDS")


# Intersect
# ------------------------------------------------------------------------------





# Plot paths
# ------------------------------------------------------------------------------


# Top level
plot_dir <- "/home/amorin/Plots"

# ChIP-seq plot files
cplot_dir <- file.path(plot_dir, "Chipseq")
cplot_subdir <- c("GRanges/", "Trackplots/", "Preprocess/", "Binding_similarity/", "Binding_summary/")

# Perturb plot files
pplot_dir <- file.path(plot_dir, "TF_perturb")
pplot_subdir <-  c("Meta_sample_matrix/", "Effect_size/", "Describe_FDR_counts/", "Experiment_similarity/")

# Integrated data plot paths
iplot_dir <- file.path(plot_dir, "Intersect/")
iplot_subdir <- c("Experiment_similarity/", "Gene_rankings")


# for (subdir in cplot_subdir) {
#   dir.create(file.path(cplot_dir, subdir), showWarnings = FALSE)
# }
# 
# 
# for (subdir in pplot_subdir) {
#   dir.create(file.path(pplot_dir, subdir), showWarnings = FALSE)
# }
# 
# 
# for (subdir in iplot_subdir) {
#   dir.create(file.path(iplot_dir, subdir), showWarnings = FALSE)
# }



# Figure exact pathing:

# F1 B
# paste0(plot_dir, "TF_perturb/Meta_sample_matrix/batch1_tf_perturb_counts_Apr2022.png")
# paste0(plot_dir, "Chipseq/Describe_meta/batch1_chip_experiment_counts_all_Apr2022.png")

# FS1 A, B, C
# paste0(plot_dir, "Chipseq/Describe_meta/batch1_chip_count_tf-idr_mecp2_overlap_peaks_by_symbol_Apr2022.png")
# paste0(plot_dir, "Chipseq/Describe_meta/batch1_chip_count_tf-idr_mecp2_overlap_peaks_by_input_Apr2022.png")
# paste0(plot_dir, "Chipseq/Describe_meta/batch1_chip_count_tf-idr_mecp2_overlap_peaks_vs_avg_exp_mreadsApr2022.png")

# F2 A, B, C, D, E
# paste0(plot_dir, "Chipseq/Binding_similarity/Densplot_all_human_Apr2022.png")
# paste0(plot_dir, "Chipseq/Binding_similarity/Densplot_all_mouse_Apr2022.png")
# paste0(plot_dir, "Chipseq/Binding_similarity/Densplot_all_ortho_Apr2022.png")
# paste0(plot_dir, "Chipseq/Binding_summary/Bind_specificity_heatmap_human_Apr2022.png")
# paste0(plot_dir, "Chipseq/Binding_summary/VBplot_bindscore_human_Apr2022.png")

# FS2 A, B, C
# paste0(plot_dir, "Chipseq/Binding_similarity/Vbplot_tf_human_intersect_Apr2022.png")
# paste0(plot_dir, "Chipseq/Binding_similarity/Vbplot_tf_mouse_intersect_Apr2022.png")
# paste0(plot_dir, "Chipseq/Binding_similarity/Vbplot_tf_ortho_intersect_Apr2022.png")

# FS3 A, B
# paste0(plot_dir, "Chipseq/Binding_similarity/Densplot_cor_runx1_ortho_Apr2022.png")
# paste0(plot_dir, "Chipseq/Binding_similarity/Densplot_inter_runx1_ortho_Apr2022.png")

# FS4 A, B
# paste0(plot_dir, "Chipseq/Binding_summary/Bind_specificity_heatmap_mouse_Apr2022.png")
# paste0(plot_dir, "Chipseq/Binding_summary/VBplot_bindscore_mouse_Apr2022.png")

# FS5 A, B
# paste0(plot_dir, "Chipseq/Binding_summary/Human_top_bound.png")
# paste0(plot_dir, "Chipseq/Binding_summary/Mouse_top_bound.png")

# F3 A, B (+/- legend, as ended up labeling axis in Adobe due to clipping)
# paste0(plot_dir, "Chipseq/GRanges/cCRE_overlap_by_exp_human_legend_Apr2022.png")
# paste0(plot_dir, "Chipseq/GRanges/cCRE_overlap_by_exp_human_nolegend_Apr2022.png")
# paste0(plot_dir, "Chipseq/GRanges/cCRE_allprop_human_Apr2022.png")
# paste0(plot_dir, "Chipseq/GRanges/cCRE_overlap_by_exp_mouse_legend_Apr2022.png")
# paste0(plot_dir, "Chipseq/GRanges/cCRE_overlap_by_exp_mouse_nolegend_Apr2022.png")
# paste0(plot_dir, "Chipseq/GRanges/cCRE_allprop_mouse_Apr2022.png")

# F6 A, B (+/- overlaid text, as ended labeling in Adobe as trackplot text too small)
# paste0(plot_dir, "Chipseq/Trackplots/Human_ASCL1_titles_chr9:38022988-38023319.png")
# paste0(plot_dir, "Chipseq/Trackplots/Human_ASCL1_notitles_chr9:38022988-38023319.png")
# paste0(plot_dir, "Chipseq/Trackplots/Human_sample-TR_titles_chr9:38022988-38023319.png")
# paste0(plot_dir, "Chipseq/Trackplots/Human_sample-TR_notitles_chr9:38022988-38023319.png")
# paste0(plot_dir, "Chipseq/Trackplots/ASCL1_SHB_trackplot_experimentIDs.tsv")

# FS7 A, B, C, D, E, F
# paste0(plot_dir, "TF_perturb/Meta_sample_matrix/Gene_measure_count_human_Apr2022.png")
# paste0(plot_dir, "TF_perturb/Meta_sample_matrix/Gene_measure_count_mouse_Apr2022.png")
# paste0(plot_dir, "TF_perturb/Meta_sample_matrix/human_sample_matrix_clustered_Apr2022.png")
# paste0(plot_dir, "TF_perturb/Meta_sample_matrix/mouse_sample_matrix_clustered_Apr2022.png")
# paste0(plot_dir, "TF_perturb/Meta_sample_matrix/batch1_tf_perturb_tech_counts_Apr2022.png")
# paste0(plot_dir, "TF_perturb/Meta_sample_matrix/batch1_tf_perturb_counts_Apr2022.png")

# F4 A, B (+/- legend to better fit main plot in figure), C 
# paste0(plot_dir, "TF_perturb/Effect_size/DEG_counts_boxplot_Apr2022.png")
# paste0(plot_dir, "TF_perturb/Effect_size/DEG_counts_byPertandSpecies_Apr2022.png")
# paste0(plot_dir, "TF_perturb/Effect_size/DEG_counts_byPertandSpecies_nolegend_Apr2022.png")
# paste0(plot_dir, "TF_perturb/Experiment_similarity/cor_heatmap_ortho_Apr2022.png")

# FS8 A, B, C, D, E
# paste0(plot_dir, "TF_perturb/Experiment_similarity/Densplot_pval-intersect_all_human_Apr2022.png")
# paste0(plot_dir, "TF_perturb/Experiment_similarity/Densplot_pval-intersect_all_mouse_Apr2022.png")
# paste0(plot_dir, "TF_perturb/Experiment_similarity/Densplot_pval-intersect_all_ortho_Apr2022.png")
# paste0(plot_dir, "TF_perturb/Experiment_similarity/Vbplot_all_human_Apr2022.png")
# paste0(plot_dir, "TF_perturb/Experiment_similarity/Vbplot_all_mouse_Apr2022.png")

# FS9 A, B, C
# paste0(plot_dir, "TF_perturb/Effect_size/Perturbed_TF_FC_Apr2022.png")
# paste0(plot_dir, "TF_perturb/Effect_size/DEG_counts_vs_FC_Apr2022.png")
# paste0(plot_dir, "TF_perturb/Effect_size/DEG_counts_vs_PercRankFC_Apr2022.png")

# FS10 A, B, C
# paste0(plot_dir, "TF_perturb/Effect_size/DEG_counts_boxplot_byperturb_allexp_Apr2022.png")
# paste0(plot_dir, "TF_perturb/Effect_size/DEG_counts_boxplot_byperturb_exclnodeg_Apr2022.png")
# paste0(plot_dir, "TF_perturb/Effect_size/hist_countde_propde_FDR=0.1_Apr2022.png")

# F5 A, B (both +/- legend to fit better in Adobe), C, D, E
# paste0(plot_dir, "TF_perturb/Describe_FDR_counts/prior_decount_bin_human_FDR=0.1_Apr2022.png")
# paste0(plot_dir, "TF_perturb/Describe_FDR_counts/prior_decount_bin_human_nolegend_FDR=0.1_Apr2022.png")
# paste0(plot_dir, "TF_perturb/Describe_FDR_counts/prior_decount_bin_mouse_FDR=0.1_Apr2022.png")
# paste0(plot_dir, "TF_perturb/Describe_FDR_counts/prior_decount_bin_mouse_nolegend_FDR=0.1_Apr2022.png")
# paste0(plot_dir, "TF_perturb/Describe_FDR_counts/TF_hist_decounts_FDR=0.1_Apr2022.png")
# paste0(plot_dir, "TF_perturb/Describe_FDR_counts/Human_ASCL1_hist_decounts_FDR=0.1_Apr2022.png")
# paste0(plot_dir, "TF_perturb/Describe_FDR_counts/Mouse_Ascl1_hist_decounts_FDR=0.1_Apr2022.png")
# paste0(plot_dir, "TF_perturb/Describe_FDR_counts/Human_ASCL1_FC_heatmap_Apr2022.png")
# paste0(plot_dir, "TF_perturb/Describe_FDR_counts/Human_ASCL1_DEprior_binary_heatmap_Apr2022.png")
# paste0(plot_dir, "TF_perturb/Describe_FDR_counts/Mouse_Ascl1_FC_heatmap_Apr2022.png")
# paste0(plot_dir, "TF_perturb/Describe_FDR_counts/Mouse_Ascl1_DEprior_binary_heatmap_Apr2022.png")

# FS11
# paste0(plot_dir, "TF_perturb/Describe_FDR_counts/human_max_countde_vs_deprior_Apr2022.png")

# FS12 A, B
# paste0(plot_dir, "TF_perturb/Describe_FDR_counts/Purity_vs_decount_bplot_mincount_FDR=0.1_Apr2022.png")
# paste0(plot_dir, "TF_perturb/Describe_FDR_counts/Human_ASCL1_example_purity_heatmap_Apr2022.png")

# FS13 A, B, C, D, E
# paste0(plot_dir, "Intersect/Experiment_similarity/Densplot_pval-intersect_all_human_Apr2022.png")
# paste0(plot_dir, "Intersect/Experiment_similarity/Densplot_pval-intersect_all_mouse_Apr2022.png")
# paste0(plot_dir, "Intersect/Experiment_similarity/Densplot_pval-intersect_all_ortho_Apr2022.png")
# paste0(plot_dir, "Intersect/Experiment_similarity/Vbplot_tf_human_intersect_Apr2022.png")
# paste0(plot_dir, "Intersect/Experiment_similarity/Vbplot_tf_mouse_intersect_Apr2022.png")

# F6 B, C, D, E, F, G (+/- legend to fit in Adobe)
# paste0(plot_dir, "Intersect/Gene_rankings/Human_ASCL1_precisionrecall.png")
# paste0(plot_dir, "Intersect/Gene_rankings/Human_ASCL1_sample_AUPRC.png")
# paste0(plot_dir, "Intersect/Gene_rankings/Human_ASCL1_single_experiment_AUPRC.png")
# paste0(plot_dir, "Intersect/Gene_rankings/Human_sampled_AUPRC_gt_observed_table.png")
# paste0(plot_dir, "Intersect/Gene_rankings/Human_aggregate_AUPRC_percentile_table.png")

# FS14
# paste0(plot_dir, "Intersect/Gene_rankings/aggregated_scores_by_curated_human_Apr2022.png")

# FS15
# paste0(plot_dir, "Intersect/Gene_rankings/aggregated_scores_by_curated_mouse_Apr2022.png")

# FS16 A, B, C
# paste0(plot_dir, "Intersect/Gene_rankings/count_curated_all.png")
# paste0(plot_dir, "Intersect/Gene_rankings/Human_AUPRC_table.png)
# paste0(plot_dir, "Intersect/Gene_rankings/Human_sampled_AUPRC_gt_observed_table.png")
# paste0(plot_dir, "Intersect/Gene_rankings/Human_aggregate_AUPRC_percentile_table.png")
# paste0(plot_dir, "Intersect/Gene_rankings/"Human_AUPRC_table.png")
# paste0(plot_dir, "Intersect/Gene_rankings/Mouse_sampled_AUPRC_gt_observed_table.png")
# paste0(plot_dir, "Intersect/Gene_rankings/Mouse_aggregate_AUPRC_percentile_table.png")
# paste0(plot_dir, "Intersect/Gene_rankings/Mouse_AUPRC_table.png")

# FS17
# paste0(plot_dir, "Chipseq/GRanges/RUNX1_Kasumi1_cCRE_overlap.png")
