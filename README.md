## This repo contains the code for "Characterizing the targets of transcription regulators by aggregating ChIP-seq and perturbation expression data sets"
https://www.biorxiv.org/content/10.1101/2022.08.30.505909v2


A note about the structure/organization of the scripts, which follows the logic of the paper:

Multiple aspects of the analysis came from rounds of curation on google sheets, and drawing on data
available on the Pavlab servers that are not readily shareable (like the aligned signal tracks
for every ChIP-seq experiment). Some data paths are hard coded to point to these sources. Summarized 
data and curated sheets used for analysis are provided wherever possible. Please contact if additional 
information is required.


### setup-

````
01
````


### chipmeta-
These scripts collect metadata from the Chip Atlas database for further curation
so that the underlying samples (fastq files) can be downloaded and submitted
to the ENCODE pipeline. Information from the resulting pipeline outputs are then
organized/consolidated with the master metadata. 
NOTE: Ideally the processes of meta org, fastq download, pipeline submission, qc
organization etc would be more modular with the ENCODE pipeline for ease of
adapting to other experiments. However, figuring this process out was entwined 
with this project/paper and the associated ChIP-seq metadata, so I'm keeping the 
"logic" laid out as executed.

```
01_download_and_process_chipatlas_metadata.R : Download latest metadata of all
ChIP-seq experiments on Chip Atlas and isolate TF ChIP-seq

02_save_batch1_chip_meta_for_curation.R : Only keep TFs of interest for export
for human curation in gsheets.

03_save_curated_input_for_encode_pipeline.R : This is a rather hacky script that
takes the curated ChIP-seq metadata, and exports text tables of fastq files to
download, and another text table for sample->experiment structure for input to
the ENCODE pipeline. 

04_organize_encode_pipeline_output.R : This associates experiment IDs with the
respective output dir from the ENCODE pipeline, checks that the expected samples
in the metadata.json produced by the pipeline match master metadata samples for
that experiment, subsets meta to successfully ran experiments, and organizes
the location of qc.json files produced for each experiment for input to the 
qc2tsv tool.

05_organize_encode_pipeline_qc.R : Takes the output of the qc2tsv tool and 
further processes. Also averages sample QC within experiments to get an 
experiment level summary. This information is than appended to create the final
metadata.

06_describe_metadata_and_qc.R : Gives overview statistics of the assembled
ChIP-seq experiments and saves out plots.
```

### chipmatrix-

```
01
```

### chipgr-

```
01
```

### trackplot-

```
01
````

### gemma-
These scripts are for organizing the expression perturbation data. These scripts
are reliant on a metadata sheet tracked in Googlesheets that recorded 
experiments. Note that this process was start before Gemma.R functionality was 
finalized, and thus experiments of interest were downloaded on Pavlab servers 
and more or less manually tracked in a sheet. Future efforts would be wise to 
access the data via Gemma.R

```
01_match_resultsets_to_experiments.R :  Reads curated table from Gsheets, flags
experiments that need to be curated/paired to a resultset ID, and saves hard
copy of metadata

02_save_tfperturb_rds.R : Reads curated table that matches resultset IDs to 
experiments. Loads, processes, and saves these experiments to a list RDS object
```

### perturbmatrix-
These scripts construct effect size matrices for analysis from the list of
perturbation experiments 

```
01_save_perturb_effect_size_matrices.R : Create tstat/fc/pval/fdr matrices for
human and mouse experiments, as well as combined matrix of ortho genes

02_describe_perturb_meta_and_mat_coverage.R : Gives overview statistics of the
assembled perturbation data and saves plots of features like gene coverage

03_describe_perturb_effect_sizes.R : Gives overview statistics and creates plots
of features like the relationship between the FCmagnitude of the perturbed TR 
and the count of DEGs

04_describe_and_save_FDR_counts.R : This arguably overstuffed script generates
information about features like the count of times a gene was DE and the 
consistency of the change of direction. This information is organized as lists
and exported out for later downstream analysis. In addition, this script provides
an overview and plots for the generated statistics.

05_perturb_experiment_similarity: Generates the pearson correlation of fold
changes between experiments as well as various metrics of the top n overlap
of genes. Summarizes in-versus-out by TR groups.
```

### intersect-
```
01
```
