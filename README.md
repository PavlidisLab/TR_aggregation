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

```
01
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
experiments. Loads, processes, and saves these experiments to an RDS object
```

### perturbmatrix-

```
01
```

### intersect-
```
01
```
