# GBC SUM159PT paper

This repository contains the code to reproduce the analysis reported in the article:

Nadalin *et al.* Multi-omic lineage tracing predicts the transcriptional, epigenetic and genetic determinants of cancer evolution. *biorxiv*. https://doi.org/10.1101/2023.06.28.546923

[Description](##description)  
[System requirements](##system-requirements)  
[Installation](##installation)  
[Instructions for reproducing the analysis](##instructions-for-reproducing-the-analysis)  
[Instructions for generating the figures](##instructions-for-generating-the-figures)  
[Contact](##contact) 

## Description

Data, scripts and pipelines are organised in three folders:

* analysis/ is divided into sub-directories, each referring to a specific analysis step. They contain the commands used to run the pipelines and information on the path containing the raw and processed data.
* scripts/ contains the pipeline scripts. They are also grouped by analysis type. 
* data/ contains external data and information required for the analysis.

General-purpose custom code is stored in separate repositories (see below).

## System requirements

To reproduce the analysis, the following software packages are required:

* CellRanger v6.0
* seqkit v2.1.0
* idr v2.0.3
* R v4.0.3
* R packages: biovizBase (v1.41.0), BSgenome.Hsapiens.UCSC.hg38 (v1.4.4), EnsDb.Hsapiens.v86 (v2.99.0), Signac (v1.4.0), cisTopic (v0.3.0)
RcisTarget (v1.13.0), AUCell (v1.15.0), ComplexHeatmap (v2.9.4), MAST (v1.19.0), clusterProfiler (v4.1.4), ReactomePA (v1.37.0), org.Hs.eg.db (v3.14.0), limma (v3.49.5), Seurat (v4.0.5), Scissor (v2.0.0)

Additionally, the following packages should be downloaded and installed (see [Installation](#installation)):

* https://gitlab.ebi.ac.uk/francesca/barcode_groups
* https://gitlab.ebi.ac.uk/francesca/fastqsplit
* https://gitlab.ebi.ac.uk/francesca/r_scripts

The code was executed on a GNU/Linux machine, using 4 cores and 64GB RAM for cell calling and 1 core and ~25GB RAM for downstream analysis.

## Installation

To run the analysis, download and save the present repository and the required packages as follows:

```
$ mkdir <MYDIR>
$ cd <MYDIR>
$ git clone https://gitlab.ebi.ac.uk/francesca/gbc_sum159pt_paper.git
$ mkdir git
$ cd git
$ git clone https://gitlab.ebi.ac.uk/francesca/barcode_groups.git
$ git clone https://gitlab.ebi.ac.uk/francesca/fastqsplit.git
$ git clone https://gitlab.ebi.ac.uk/francesca/r_scripts.git
```

where ```<MYDIR>``` is the directory where analysis and scripts will be stored.
Make sure that both R and Rscript executables are in your PATH and install the required libraries as specified in [System requirements](#system-requirements).

All but steps 0 and 2 were executed on HPC infrastructures.
Steps 1 and 3 were run on torque architecture, whereas and 4 to 9 on lsf, hence the syntax of the submission commands is slightly different.

## Instructions for reproducing the analysis

The analysis is organised in 10 steps. There are dependencies between them, so the analysis should be run in the same order as reported below.
We report all the steps, including the code used for generating the feature-barcode matrices. If those are provided, go directly to step 4.

### 0. Data download

Data accession codes are reported in our manuscript and will be publicly available upon publication.
If you plan to use the feature-barcode matrices we generated, these should be stored on specific locations, for example:

```
$ cd gbc_sum159pt_paper/analysis/cell_calling
$ mkdir -p 1B_GEX/T0_1/outs/filtered_feature_bc_matrix
$ cd 1B_GEX/T0_1/outs/filtered_feature_bc_matrix
$ wget <URL>
$ ls
barcodes.tsv.gz  features.tsv.gz  matrix.mtx.gz 
```

where ```<URL>``` is the url pointing to the files.
The above command downloads the feature-barcode expression matrix for baseline (T0) cells, replicate 1. To proceed with the analysis, go directly to step 4.
If instead you want to recreate the gene expression profiles from the raw data, you should first download the FASTQ files and then provide the directory name to step 1-3 in a specific file (see below).


### 1. GBC calling

GBC sequences are not known in advance and hence have to be determined. On infected cells, we performed both scRNA-Seq, which captures boh endogenous and GBC-carrying transcripts, and targeted DNA-seq of the reporter locus.

Samples are grouped by experiment.
FASTQ files should be downloaded and their location provided in a tab-separated file for each experiment, for example:

```
$ cd gbc_sum159pt_paper/analysis/GBC_calling
$ cat files_raw_exp_1B_DNA.tsv
D0	/path/to/Sample_S22108_D0
D3	/path/to/Sample_S22109_D3
D9	/path/to/Sample_S22110_D9
D13	/path/to/Sample_S22111_D13
D17	/path/to/Sample_S22112_D17
D24	/path/to/Sample_S22113_D24
```

where ```/path/to/``` should be replaced with the actual directory name.
```files_raw_exp_1B_DNA.tsv``` contains the sample ID on the first column and the location of the FASTQ files on the second column.
File ```files_counts_exp_1B_DNA_amplicon_v2.tsv``` is provided and contains the name of the output files for that specific input; it should not be modified unless you wish to rename the samples.

The GBC calling pipeline works as follows.
First, starting from the FASTQ files (DNA or RNA), the set of sequences lying in the GBC locus is determined.
Second, sequences are error-corrected, *i.e.*, sequences that differ from each other by < d sequencing errors are grouped together. The top abundant sequence in each group represents the original GBC and is kept for the subsequent analysis steps:

```
$ cd gbc_sum159pt_paper/analysis/GBC_calling
$ for i in 1 1_invivo_treated 4 5 
> do
> bash RUN_barcode_detection_exp_${i}_amplicon_v2.sh
> bash RUN_barcode_groups_exp_${i}_amplicon_v2.sh
> done
```

The experiment number represent the infection, so that GBC IDs from samples in the same experiment refer to the same GBC sequence.
Note that the order of the experiments matters ("1" must be executed before "1\_invivo\_treated").

The output will contain the list of GBC sequences and the associated statistics (see [barcode_groups](https://gitlab.ebi.ac.uk/francesca/barcode_groups) for details).

### 2.  Clone estimate (bulk samples)

Estimate the relative clone content from bulk DNA-Seq data from GBC fraction in samples. Here, we assume that a GBC is equivalent to a clone. Then, we define the tumour-initiating clones (TICs) and the drug-tolerant clones (DTC) according on their occurrence in tumours (for TICs) and > 13 days post-treatment (for DTCs).
The pipeline is run as follows:

```
$ cd gbc_sum159pt_paper/analysis/clone_estimate_bulk
$ bash RUN_clone_estimate.sh
$ bash RUN_DRC_TIC_prediction.sh
```

### 3. Cell calling and sample demultiplexing

Do cell calling using cellranger (for scRNA-Seq) and cellranger-arc (for Multiome). This is supposed to be the first step to the definition of which barcodes are actual cells and which ones are not. See 4. and 5. for post-processing details.

The GBC list determined at step 1 is used as a feature barcoding library. 
Note that some samples have been multiplexed (MULTI-Seq), so a demultiplexing step is added to properly assign cells to samples:

```
$ cd gbc_sum159pt_paper/analysis/cell_calling
$ bash RUN_cellranger.sh
$ bash RUN_demultiplex.sh
```

Given the high noise level observed in some samples and in order to recover more cells, demultiplexing is performed with a slightly modified version of [deMULTIplex](https://github.com/chris-mcginnis-ucsf/MULTI-seq). The difference between the two approaches essentially lies in how the MULTI-Seq barcode coverage across cells is evaluated.

Only cell calling, for MDA-MB-231:

```
$ cd gbc_sum159pt_paper/analysis/cell_calling_MDA
$ bash RUN_cellranger.sh
```

### 4. Clone calling (single-cell samples)

From step 4 on, we will use the Seurat package. Feature-cell expression matrices are stored in Seurat objects.

Call clones from GBC expressed in single cells. A clone here refers to a set of GBCs and cells infected with the same set of GBCs are assigned to the same clonal population. The clone content of each sample is calculated as the number of cells assigned to that clone. Note that a cellular barcode can encode for 0, 1 or > 1 clones.

First, the set of expressed GBCs in a cell are determined by removing low-abundance GBCs, which likely belong to the soup and not to the cell's genome. This step discriminates cells infected by 0 or > 0 GBCs. 
Cellular barcodes associated to > 1 GBCs are then classified into co-infection (one cell carrying multiple GBCs) and doublet events (identical CB assigned to distinct cells).
Finally, clone content is calculated for each sample and compared across samples.

The pipeline is run as follows:

```
$ cd gbc_sum159pt_paper/analysis/clone_calling_sc
$ bash RUN_clone_calling.sh
```

### 5. Barcoded scRNA-seq data analysis

Analysis of scRNA-Seq data in GBC-infected samples.
First, doublets and cells carrying no GBCs (as determined at step 4. above) are discarded. On baseline samples, we also regress out the cell-cycle score.

Feature selection, dimensional reduction and clustering are done on the gene-cell expression matrices obtained. Multiple clustering solutions are generated, by varying the input parameters. These are evaluated using the silhouette score. Several plots are generated: sample composition, UMAP, TSNE, PCA plots.

One clustering solution is determined at each iteration, and for each we remove the cluster with much lower UMI count compared to its complement, until the UMI count statistics is comparable across clusters (which is fine since data come from a cell line).
This concludes the cell calling procedure.

We do differential expression on each cluster compared with its complement and possibly do a functional annotation on its top up-regulated genes.
Finally, clones are grouped into lineages using an unsupervised approach that computes clone-clone gene expression similarity. Lineages are defined as sets of clones that are grouped together, by transcriptional similarity, in indipendent replicates. Finally, we determine a gene expression signature for each lineage.

The pipeline is run as follows:

```
$ cd gbc_sum159pt_paper/analysis/gene_expression
$ bash RUN_clustering_ccRegress.sh
$ bash RUN_clustering_ccRegress_filt.sh
$ bash RUN_clustering_ccRegress_filt_2ndRound.sh
$ bash RUN_DEA_ccRegress.sh
$ bash RUN_lineage.sh
```

For MDA-MB-231:

```
$ cd gbc_sum159pt_paper/analysis/gene_expression_MDA
$ bash RUN_clustering_ccRegress.sh
$ bash RUN_clustering_ccRegress_filt.sh
$ bash RUN_DEA_ccRegress.sh
```

### 6. Barcoded multiome data analysis

Multi-omic (RNA+ATAC) analysis of clonal populations at single-nucleus resolution. 
First, clone information is added to each nucleus:

```
$ cd gbc_sum159pt_paper/analysis/multiome
$ bash RUN_clone_assignment.sh
```

Then, RNA and ATAC are analysed separately. See 

```
gbc_sum159pt_paper/analysis/multiome/workflow.md 
```

for details.

### 7. Identification of transcriptionally stable clones

Computation of the clone sharedness score between clusters across samples and extraction of the transcriptionally stable clones. Pairs of clusters that share a high fraction of clones between the two replicates at baseline are called S1, S2, and S3. We obtain a gene expression signature for each, which is then annotated (gene ontology, pathway enrichment analysis).

The distance between and within pre-defined groups (phenotypes) is computed in specific samples. Finally, we compare clone abundance as measured in bulk (GBC fraction, see 2.) and in single-cell samples (cell count, see 4.).

The pipeline is run as follows:

```
$ cd gbc_sum159pt_paper/analysis/sister_cells
$ bash RUN_sister_similarity_ccRegress.sh
$ bash RUN_sister_similarity_treated.sh
$ bash RUN_transcriptionally_identifiable.sh
$ bash RUN_phenotype_distance.sh
$ bash RUN_DRC_TIC_analysis.sh
```

### 8. Whole-exome sequencing data analysis

Using CNVkit, we obtained the list of predicted copy-number variations (CNV) between treated and untreated conditions. A consensus CNV is obtained across multiple replicates, asking that the CNV locus should be shared by > 80%. We then compute the read coverage on CNV loci in scATAC-Seq samples:

```
$ cd gbc_sum159pt_paper/analysis/wes
$ bash RUN_CNV_analysis.sh
$ bash RUN_CNV_double_treatment_analysis.sh
```

### 9. Analysis of public datasets

#### 9.0. scRNA-Seq of primary breast tumours [Pal *et al.*, 2021]

This steps aims to detect S1, S2 and S3 signatures in scRNA-Seq data from primary TNBC patients from [Pal *et al.*, 2021]. We analyse samples 0106, 0114, 0125 and 0126.

The dataset should be downloaded from GEO (accession: GSE161529) and stored in ```gbc_sum159pt_paper/data/Visvader``` according to the paths in ```gbc_sum159pt_paper/analysis/Visvader/matrix_dir_list.tsv```, as follows:

```
$ ls gbc_sum159pt_paper/data/Visvader/0106
barcodes.tsv.gz  features.tsv.gz  matrix.mtx.gz
```

Using Scissor, we detect "phenotype-positive" for each signature S1, S2 and S3 respectively:

```
$ cd gbc_sum159pt_paper/analysis/Visvader
$ bash RUN_Visvader_analysis.sh
```

#### 9.1. Curated Cancer Cell Atlas (3CA) meta-programs [Gavish *et al.*, 2023]

We downloaded the top 50 significant genes for the tumour meta-programs from [Gavish *et al.*]. The table is stored in data/3CA.

To run the analysis:

```
$ cd gbc_sum159pt_paper/analysis/3CA
$ bash RUN_3CA_analysis.sh
```

Similarly, for MDA-MB-231:

```
$ cd gbc_sum159pt_paper/analysis/3CA_MDA
$ bash RUN_3CA_analysis.sh
```

#### 9.2. TCGA

```
TODO
```

## Instructions for generating the figures

First, we need to download the Seurat objects and the tables:

```
TODO
```


## Contact

For feedback or questions about this repository, please contact [Francesca Nadalin](mailto:francesca@ebi.ac.uk). 
