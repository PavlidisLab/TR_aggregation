## "Characterizing the targets of transcription regulators by aggregating

## ChIP-seq and perturbation expression data sets"

<https://www.biorxiv.org/content/10.1101/2022.08.30.505909v2>

Note that multiple aspects of the analysis came from rounds of curation on google sheets, and drawing on data available on the Pavlab servers that are not readily shareable (like the aligned signal tracks for every ChIP-seq experiment). Some data paths are hard coded to point to these sources. Summarized data and curated sheets used for analysis are provided wherever possible. Please contact if additional information is required.

### setup-

Preliminary organization of necessary information.

    01_config.R : Sets global variables and defines output paths.

    02_install_packages.R : R packages used for analysis.

    03_download_genomic_tables.R : Download and format information like protein coding genes and ENCODE cCREs.

    04_organize_curated_targets.R : Load, format, and export table of literature curated targets from gsheets.

### chipmeta-

These scripts collect metadata from the Chip Atlas database for further curation so that the underlying samples fastq files can be downloaded and submitted to the ENCODE pipeline. Information from the resulting pipeline outputs are then organized/consolidated with the master metadata. Ideally the processes of meta organization, fastq download, pipeline submission, qc organization etc would be more modular with the ENCODE pipeline for ease of adapting to other experiments. However, figuring this process out was entwined with this paper and finalizing the ChIP-seq metadata, so I'm keeping the "logic" laid out as executed.

    01_download_and_process_chipatlas_metadata.R : Download latest metadata of all ChIP-seq experiments on Chip Atlas and isolate TF experiments.

    02_save_batch1_chip_meta_for_curation.R : Only keep TFs of interest for 
    export for human curation in gsheets.

    03_save_curated_input_for_encode_pipeline.R : Hacky script that uses the curated ChIP-seq metadata to exports text tables of fastq files to download, and another text table for sample->experiment structure for input to the ENCODE pipeline. Candidate to improve.

    04_organize_encode_pipeline_output.R : This associates experiment IDs with the respective output dir from the ENCODE pipeline, checks that the expected samples in the metadata.json produced by the pipeline match master metadata samples for that experiment, subsets meta to successfully ran experiments, and organizes the location of qc.json files produced for each experiment for input to the qc2tsv tool.

    05_organize_encode_pipeline_qc.R : Takes the output of the qc2tsv tool and further processes. Also averages sample QC within experiments to get an 
    experiment level summary. This information is than appended to create the final metadata.

    06_describe_metadata_and_qc.R : Gives overview statistics of the assembled
    ChIP-seq experiments and saves out plots.

### chipmatrix-

These scripts score each ChIP-seq experiment into a gene by experiment matrix, which is then used for analysis.

    01_build_encpipe_bindingscore_mat.R : Create a gene by experiment matrix of continuous binding scores for export.

    02_build_encpipe_binary_mat.R : Create a gene by experiment matrix of binary binding scores for export.

    03_explore_and_savepreprocess.R : Filter and process the binding matrices for export, and plot effect of processing steps.

    04_experiment_similarity.R : Describe and plot several different measures of similarity between the ChIP-seq experiments.

    05_binding_summary.R : Summarize binding across experiments and fit a model that looks for TF-specific binding. Saves out associated plots and the data summaries.

### chipgr-

These scripts organize and analyze ChIP-seq experiments at the loci level as GenomicRange objects.

    01_save_grlist.R : Exports all experiments as GRList objects grouped by species.

    02_explore_overlap.R : Summarize peak sizes and distributions before/after re-sizing and reducing, and generate matrices that tally TR overlap by regions.

    03_overlap_ccre.R : Explores the overlap of ChIP-seq peaks with ENCODE candidatecis regulatory elements.

### trackplot-

For visualizing the aligned signal tracks of ChIP-seq experiments.

    01_ASCL1_SHB.R : Exports a plot of an intronic region in SHB bound by every human ASCL1 experiment but rarely bound by the other TRs

### gemma-

These scripts are for organizing the expression perturbation data. These scripts are reliant on a metadata sheet tracked in Googlesheets that recorded experiments. Note that this process was started before Gemma.R functionality was finalized, and thus experiments of interest were downloaded on Pavlab servers and more or less manually tracked in a sheet. Future efforts would be wise to access the data via Gemma.R

    01_match_resultsets_to_experiments.R : Reads curated table from Gsheets, flags experiments that need to be curated/paired to a resultset ID, and saves hard copy of metadata.

    02_save_tfperturb_rds.R : Reads curated table that matches resultset IDs to experiments. Loads, processes, and saves these experiments as a list.

### perturbmatrix-

These scripts construct effect size matrices for analysis from the list of perturbation experiments.

    01_save_perturb_effect_size_matrices.R : Create tstat/fc/pval/fdr matrices for each species as well as combined matrix of ortho genes.

    02_describe_perturb_meta_and_mat_coverage.R : Gives overview statistics of the assembled perturbation data and saves plots of features like gene coverage.

    03_describe_perturb_effect_sizes.R : Gives overview statistics and creates plots of features like the relationship between the FC magnitude of the perturbed TR and the count of DEGs.

    04_describe_and_save_FDR_counts.R : Summarizes gene effect sizes across experiments for export and plotting.

    05_perturb_experiment_similarity: Describe and plot several different measures of similarity between the perturbation experiments.

### intersect-

These scripts takes the summarized outputs from the perturbation and ChIP-seq experiments to compare, combine, and ultimately rank gene targets.

    01_saved_merged_rankings.R : Merge the summary data frames of the ChIP-seq and perturbation experiments generated for each TR into a single table of gene rankings.

    02_chip_perturb_similarity : Describe and plot several different measures of similarity between the perturbation and ChIP-seq experiments, and export an object containing the overlapping top scoring genes.

    03_gene_rankings_curated_evaluation.R : Test and plot the ability of the aggregated rankings to recover literature curated targets.

    04_explore_rankings.R : For interactive exploration and plotting of the ranking tables.
