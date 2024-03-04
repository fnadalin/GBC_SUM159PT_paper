# Documentation

## File `Run_VariableFeatures.R`

### Description
This script takes a Seurat object as input, normalizes the data and calculates the variable features for each modality. Normalization for RNA is standard, while for ATAC is TF-IDF

### Usage
`Usage: Run_VariableFeatures.R [options]`

Options:

	--Seurat=SEURAT
		[REQUIRED] Seurat Object path

	--out_RDS=OUT_RDS
		[OPTIONAL] Output file name, otherwise input file name will be used

	--RNA=RNA
		[OPTIONAL] Number of RNA features to include in the set [default=3000]

	--ATAC=ATAC
		[OPTIONAL] Number of ATAC features to include in the set [default=5000]


### Notes
Part of the pipeline ran in the shell script `shell/1{C,E}_ARC_varfeat_preproc_clust.sh`

## File `Run_clustering.R`

### Description
Takes as input a Seurat object and perform a clustering of the data (or a reduction of the data itself)

### Usage
`Usage: Run_clustering.R [options]`

Options:

	--Seurat=SEURAT
		[REQUIRED] Seurat Object path

	--assay=ASSAY
		[OPTIONAL] Assay from which select variable features. Only active when --features is active [default=RNA]

	--out_RDS=OUT_RDS
		[OPTIONAL] Output file name, otherwise input file name will be used

	--reduction=REDUCTION
		[OPTIONAL] Reduction to be used to construct the SNN graph with FindNeighbours function (PCA, MOFA, cisTopic ecc). Set to NULL if --features is used [default=pca]

	--algorithm=ALGORITHM
		[OPTIONAL] Value for the 'algorithm' parameter in the FindClusters Seurat function [default=3]

	--resolution=RESOLUTION
		[OPTIONAL] Value for the 'resolution' parameter in the FindClusters Seurat function [default=0.5]

	--features=FEATURES
		Features to select for clustering. Only works if reduction is NULL [choose between geneBasis or VariableFeatures



### Notes

Part of the pipeline ran in the shell script `shell/1{C,E}_ARC_varfeat_preproc_clust.sh`


## File `Run_genebasis.R`

### Description
Takes as input a Seurat object and outputs a Seurat object with the results of geneBasis saved as a reduction slot

### Usage
Usage: Run_genebasis.R [options]

Options:

	--Seurat=SEURAT
		[REQUIRED] Seurat Object path

	--out_RDS=OUT_RDS
		[OPTIONAL] Output file name, otherwise input file name will be used

	--n_genes=N_GENES
		[OPTIONAL] Maximum number of genes to search [default=100]
### Notes
This script has not been utilized in the final analysis, but works well and was tested. It is quite time and resource intensive

## File `Seurat_load.R`

### Description
Takes a CellRanger output matrix folder as input and generate a Seurat object as output

### Usage
`Usage: Seurat_load.R [options]`


Options:

	--matrix_dir=MATRIX_DIR
		[REQUIRED] path to the input matrix directory

	--GBC=GBC
		[REQUIRED] path to the input genomic barcodes (GBC) file

	--out_RDS=OUT_RDS
		[REQUIRED] output seurat object containing the cell barcodes from the matrices

	--frag_path=FRAG_PATH
		[REQUIRED] fragment file path

	--min_genes=MIN_GENES
		[OPTIONAL] minimum number of detected genes to keep a cell [default=0]

	--min_cells=MIN_CELLS
		[OPTIONAL] minimum number of cells where a gene is detected to keep it [default=0]


### Notes


## File `Seurat_preprocess.R`

### Description
Takes a Seurat object as input and performs subsetting, normalization, scaling

### Usage
`Usage: Seurat_preprocess.R [options]`

Options:

	--Seurat=SEURAT
		[REQUIRED] Seurat Object path

	--subset=SUBSET
		[OPTIONAL] logical expression with the subset of cells 
              to filter


### Notes
Part of the pipeline ran in the shell script `shell/1{C,E}_ARC_varfeat_preproc_clust.sh`


## File `bubble_plot_AUC_top_matching_topics.R`

### Description
Plot file. Bubble plot with all clusters for a single replicate and for both replicates for S1 S2 and S3 with consistency. Produces two bubble plots with NO consistency information of the first 30 topics sorted by reproducibility score and a bubble plot with just the three conserved clusters.
### Usage
`Usage: bubble_plot_AUC_top_matching_topics.R [options]`


Options:

	--matchingtopics=MATCHINGTOPICS
		[REQUIRED] File containing the list of matching topics with repr. score and AUC information.

	--out=OUT
		[REQUIRED]  Output file names separated by a comma (3 png files, replicate1, replicate2, both)
### Notes
Suited exclusively for file inputs from cistopic_idr.R --> calculate_IDR_repr_score.R --> merge_sizelist_cistopic_AUC.R. Requires a library that I did not install on the cluster (ggnewscale)

## File `calculate_IDR_repr_score.R`

### Description
It takes as input the result file of cistopic_idr.R (size_list) and return the same file with a new column (named repr.score) attached.
### Usage
`Usage: calculate_IDR_repr_score.R [options]`

Options:

	--sizelist=SIZELIST
		[REQUIRED] Output file of cistopic_idr.R
### Notes


## File `calculate_circos_plots.R`

### Description
Plot script to produce circos plots.
### Usage
Only runs locally on my machine, needs a small fixing (and setting up an actually scriptable implementation)
### Notes


## File `calculate_heatmap_IDR_score_1C_1E.R`

### Description
Plot Topic-Topic heatmap showing Reproducibility score

### Usage
`Usage: calculate_heatmap_IDR_score_1C_1E.R [options]`


Options:

	--cisTopic1=CISTOPIC1
		[REQUIRED] cisTopic1 Object path

	--cisTopic2=CISTOPIC2
		[REQUIRED] cisTopic1 Object path

	--idr_file=IDR_FILE
		[REQUIRED] Path to the result of calculate_IDR_repr_score.R

	--meta1=META1
		[REQUIRED] Path to the metadata file of experiment 1

	--meta2=META2
		[REQUIRED] Path to the metadata file of experiment 2
    
### Notes
This script actually calculates correlation with coverage on-the-spot and works only for 1C and 1E (maybe I will produce a twin script for 1C-1C and 1E-1E)

## File `calculate_heatmap_region_probabilities_topics.R`

### Description
This script produces the region probability heatmap with signatures of the 3 subpopulations.
### Usage
`Usage: calculate_heatmap_region_probabilities_topics.R`
### Notes
paths to bed files are specified in the script. Signatures file too (this may change). Signature file is produced by running `calculate_signatures.R`


## File `calculate_marker_proximity_IDR.R`

### Description
Adds marker distances to bed file with reproducible regions between topic pairs
### Usage
`Usage: calculate_marker_proximity_IDR.R [options]`


Options:

	--idr=IDR
		[REQUIRED] BED file containing regions that passed IDR filtering

	--signatures=SIGNATURES
		[REQUIRED] Signature file obtained by calculate_signatures.R

	--out=OUT
		[REQUIRED] Output file name

### Notes



## File `calculate_region_based_AUC_marker_proximity.R`

### Description
Locally working Bubble plot stub. Not usable for analysis.
### Usage
Runs locally on my machine only.
### Notes

## File `calculate_roc_curves_topics_clusters.R`

### Description
Plot script to produce ROC curves associated to a set of epigenetic modules
### Usage
Runs locally on my machine
### Notes
This needs to become a script.


## File `calculate_signatures.R`

### Description
Script to fuse all signatures and get their coordinates into a single csv file
### Usage
`Usage: calculate_signatures.R [options]`


Options:

	--signatures=SIGNATURES
		[REQUIRED] signature file names as comma-separated list

	--names=NAMES
		[REQUIRED] Signatures names as comma-separated list

	--out=OUT
		[OPTIONAL] Output filename

### Notes


## File `calculate_umap_top_regions_tfidf_topics.R`

### Description
Calculate UMAPs for top regions in each topic by visualizing TF-IDF
### Usage
Runs locally on my machine only
### Notes
This needs to become a script

## File `cisTopic_plotting.R`

### Description
Returns a folder with UMAPs colored by topic z-score for each topic.
### Usage
`Usage: cisTopic_plotting.R [options]`


Options:

	--Seurat=SEURAT
		[REQUIRED] Seurat Object path
### Notes
Run cisTopic on the Seurat object you're gonna pass to this script.

## File `cisTopic_training.R`

### Description
Takes a Seurat object and perform cisTopic analysis on it, extract the cell-topics matrix and attaches it to the original seurat object as a reduction. The full cistopic object is stored in the misc slot of the reduction of the seurat object.

### Usage
`Usage: cisTopic_training.R [options]`


Options:

	--Seurat=SEURAT
		[REQUIRED] Seurat Object path

	--Seurat2=SEURAT2
		[REQUIRED] Seurat Object path with which select overlapping regions only

	--out_RDS=OUT_RDS
		[OPTIONAL] Output file name, otherwise input file name will be used

	--topics=TOPICS
		[OPTIONAL] Range and/or comma-separated of topics number to test [start:end=20:50]

	--n_cpus=N_CPUS
		[OPTIONAL] Number of CPUs to be utilized during model training

	--min_cutoff=MIN_CUTOFF
		[OPTIONAL] Minimum cutoff for regions to be included in the cisTopic analysis (available as variable features, input 'none' if Variable Features)

	--n_iter=N_ITER
		[OPTIONAL] Total number of iterations for the algorithm

	--n_burnin=N_BURNIN
		[OPTIONAL] Number of burn-in iterations

	--var_feat=VAR_FEAT
		[OPTIONAL] Number of variable features to be included in the cisTopic analysis (available as variable features)

	--clean_regions=CLEAN_REGIONS
		[OPTIONAL] Removes the lowly accessible regions included when running feature selection


### Notes
This is time and resource intensive. For some weird reasons, allocate on the cluster **double** the CPUs you'll be passing to the script, otherwise it will crash.


## File `cistopic_calculate_auc.R`

### Description
Script to calculate AUC between topic-cell probabilities and clusters
### Usage
`Usage: cistopic_calculate_auc.R [options]`

Options:

	--cisTopic1=CISTOPIC1
		[REQUIRED] cisTopic Object 1 path (must have same topic search space size as cisTopic2)

	--cisTopic2=CISTOPIC2
		[REQUIRED] cisTopic Object 2 path (must have same topic search space size as cisTopic1)

	--metadata1=METADATA1
		[REQUIRED] File containing Seurat metadata information for object cisTopic1

	--metadata2=METADATA2
		[REQUIRED] File containing Seurat metadata information for object cisTopic2

	--clusters1=CLUSTERS1
		[REQUIRED] Column name of Seurat clustering for object cisTopic1

	--clusters2=CLUSTERS2
		[REQUIRED] Column name of Seurat clustering for object cisTopic2

### Notes
It iterates over all possible topic numbers combinatorially and looks up files from `cistopic_match_experiments.R`. Needs urgent fixing.

## File `cistopic_calculate_auc_IDR.R`
### Description
File to calculate AUC for cisTopic topic-cell probabilities against clusters

### Usage
`Usage: cistopic_calculate_auc_IDR.R [options]`

Options:

	--cisTopic1=CISTOPIC1
		[REQUIRED] cisTopic Object 1 path (must have same topic search space size as cisTopic2)

	--cisTopic2=CISTOPIC2
		[REQUIRED] cisTopic Object 2 path (must have same topic search space size as cisTopic1)

	--metadata1=METADATA1
		[REQUIRED] File containing Seurat metadata information for object cisTopic1

	--metadata2=METADATA2
		[REQUIRED] File containing Seurat metadata information for object cisTopic2

	--clusters1=CLUSTERS1
		[REQUIRED] Column name of Seurat clustering for object cisTopic1

	--clusters2=CLUSTERS2
		[REQUIRED] Column name of Seurat clustering for object cisTopic2

	--out=OUT
		[REQUIRED] Output file name
### Notes
This is the current version of the script to obtain AUC from clusters.
This file can be merged with the output of the script that matches topics based on reproducibility score to obtain the Bubble plot.


## File `cistopic_enrichment_heatmaps.R`

### Description
Old script to measure correlation between topics of two experiments.
### Usage
`Usage: cistopic_enrichment_heatmaps.R [options]`


Options:

	--Seurat1=SEURAT1
		[REQUIRED] Seurat Object 1 path

	--Seurat2=SEURAT2
		[REQUIRED] Seurat Object 2 path

	--clustering_name1=CLUSTERING_NAME1
		[REQUIRED] name of the metadata column of the clusters

	--clustering_name2=CLUSTERING_NAME2
		[REQUIRED] name of the metadata column of the clusters
### Notes

## File `cistopic_idr.R`

### Description
Calculate IDR scores between two cisTopic experiments. Takes as input two cisTopic objects and return for each pair of topics a bed file containing IDR values
### Usage
`Usage: cistopic_idr.R [options]`


Options:

	--cisTopic1=CISTOPIC1
		[REQUIRED] cisTopic1 Object path

	--cisTopic2=CISTOPIC2
		[REQUIRED] cisTopic1 Object path

	--idr_threshold=IDR_THRESHOLD
		[OPTIONAL] IDR threshold for lowly reproducible regions

	--ncpus=NCPUS
		[OPTIONAL] Number of cores to use in the IDR calculation [default=4]

	--skip_bed_writing=SKIP_BED_WRITING
		[OPTIONAL] Logical to indicate whether to write topic-specific bed files with associated region probabilities

	--skip_IDR_calculations=SKIP_IDR_CALCULATIONS
		[OPTIONAL] Logical to indicate whether to skip all calculations and prepare only the sizelist
### Notes
Time intensive, does not run on MacOS machines

## File `cistopic_match_experiments.R`

### Description
Old script to match experiments based on overlap.
### Usage
`Usage: cistopic_match_experiments.R [options]`


Options:

	--cisTopic1=CISTOPIC1
		[REQUIRED] cisTopic Object 1 path (must have same topic search space size as cisTopic2)

	--cisTopic2=CISTOPIC2
		[REQUIRED] cisTopic Object 2 path (must have same topic search space size as cisTopic1)

	--markers_file=MARKERS_FILE
		[REQUIRED] csv file containing marker genes genomic coordinates

	--thrP=THRP
		[OPTIONAL] Probability threshold for the Gamma Fit of regions for each topic in the binarization step [default=0.95]

	--metadata1=METADATA1
		[REQUIRED] File containing Seurat metadata information for object cisTopic1

	--metadata2=METADATA2
		[REQUIRED] File containing Seurat metadata information for object cisTopic2
### Notes


## File `clean_corr_coverage_topics.R`

### Description
clean sizelist from topics that correlate with coverage over a specific threshold
### Usage
`Usage: clean_corr_coverage_topics.R [options]`


Options:

	--idr_file=IDR_FILE
		[REQUIRED] Path to the result of cistopic_idr.R

	--cisTopic1=CISTOPIC1
		[REQUIRED] Topics numbers from replicate 1 that do correlate with coverage

	--cisTopic2=CISTOPIC2
		[REQUIRED] Topics numbers from replicate 1 that do correlate with coverage

	--out=OUT
		[REQUIRED] Output file name

	-h, --help
		Show this help message and exit

### Notes

## File `cluster_specific_umap_and_markers_heatmap.R`

### Description
Poster-specific UMAPs
### Usage
Runs locally
### Notes


## File `extract_cisTopic_regions.R`

### Description
Binarize cisTopic analysis and attaches results to a Seurat Object
### Usage
`Usage: extract_cisTopic_regions.R [options]`


Options:

	--Seurat=SEURAT
		[REQUIRED] Seurat Object path

### Notes

## File `extract_cluster_specific_signatures.R`

### Description
Obtain a set of files with signatures from a single cluster.
### Usage
`Usage: extract_cluster_specific_signatures.R [options]`


Options:

	--signatures_file=SIGNATURES_FILE
		[REQUIRED] Path to the result of FindMarkerFeatures function, saved as a .csv file

### Notes


## File `find_marker_genes.R`

### Description
Calculates markers for each cluster.
### Usage
`Usage: find_marker_genes.R [options]`


Options:

	--Seurat=SEURAT
		[REQUIRED] Seurat Object path

	--clustering_name=CLUSTERING_NAME
		[REQUIRED] name of the metadata column of the clusters

	--assay=ASSAY
		[OPTIONAL] Specifiy assay to select features (genes or regions)

	--test_use=TEST_USE
		[OPTIONAL] name of the statistical test to be used
### Notes


## File `merge_IDR_sizelist_cistopic_AUC.R`

### Description
Merge IDR sizelist with cisTopic AUC
### Usage
`Usage: merge_IDR_sizelist_cistopic_AUC.R [options]`


Options:

	--sizelist=SIZELIST
		[REQUIRED] Output file of cistopic_idr.R

	--auc=AUC
		[REQUIRED] Output file of cistopic_calculate_auc_IDR.R

	--out=OUT
		[REQUIRED] Output file name
### Notes
