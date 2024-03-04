# Workflows

---

## From cellranger output to UMAP coloured with clusters

### 1. Seurat_load.R

#### Input
```
MATRIX_DIR=../cell_calling/1{C,E}_ARC/P0_Multiome/outs/filtered_feature_bc_matrix/
GBC=1{C,E}_ARC/CB_metadata.tsv
FRAGPATH=../cell_calling/1{C,E}_ARC/P0_Multiome/outs/atac_fragments.tsv.gz
```
#### Script
```
Rscript ../../scripts/multiome/Seurat_load.R --matrix_dir $MATRIX_DIR --out_RDS $OUTFILE --GBC $GBC --frag_path $FRAGPATH
```

#### Output
`OUTFILE=1{C,E}_ARC_seurat.Rds`

---

### 2. Run_VariableFeatures.R

#### Input
`SEURAT=1{C,E}_seurat.Rds`

#### Script
`Rscript ../../scripts/multiome/Run_VariableFeatures.R --Seurat $SEURAT`

#### Output
`OUTFILE=1{C,E}_ARC_seurat.Rds` 

---

### 3. Seurat_preprocess.R

#### Input
`SEURAT=1{C,E}_seurat.Rds`

#### Script
`Rscript ../../scripts/multiome/Seurat_preprocess.R --Seurat $SEURAT`

#### Output
`OUTFILE=1{C,E}_ARC_seurat.Rds` 

---

### 4. Run_clustering.R

#### Input
`SEURAT=1{C,E}_ARC_seurat.Rds`

#### Script
```
Rscript ../../scripts/multiome/Run_clustering.R --Seurat $SEURAT --reduction "pca_rna" \
--resolution "0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.2"
Rscript ../../scripts/multiome/Run_clustering.R --Seurat $SEURAT --reduction "lsi_atac" \
--resolution "0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.2"
```

#### Output
```
OUTFILE=1{C,E}_ARC_seurat.Rds 
plots/1{C,E}_ARC_seurat_clustering_{pca_rna, lsi_atac}_{RESOLUTION}_clusters.png
plots/1{C,E}_ARC_seurat_clustering_{pca_rna, lsi_atac}_{RESOLUTION}_clusters_silhouette_plot.png
```

---

## From Seurat object to epigenetic modules

---

### 1. cisTopic_training.R

#### Input
```
SEURAT=cisTopic/cisTopic_all/1C_ARC_seurat.Rds
SEURAT2=cisTopic/cisTopic_all/1E_ARC_seurat.Rds
NCORES=5
```

#### Script
```
Rscript ../../scripts/multiome/cisTopic_training.R --Seurat $SEURAT \
--Seurat2 $SEURAT2 --n_cpus $N_CORES --topics 10,20,30,40,50 --n_iter 1000 \
--n_burnin 750 --min_cutoff 0
```

#### Output
```
cisTopic/cisTopic_all/1{C,E}_ARC_seurat.Rds
cisTopic/cisTopic_all/1{C,E}_ARC_seurat_cisTopicObject.Rds
cisTopic/cisTopic_all/1{C,E}_ARC_seurat_cisTopic_loglikelihood_plot.png
cisTopic/cisTopic_all/1{C,E}_ARC_seurat_cisTopic_model_selection.png
```

---

### 2. cisTopic_plotting.R

#### Input

`cisTopic_all/1{C,E}_ARC_seurat.Rds`

#### Script

`Rscript ../../scripts/multiome/cisTopic_plotting.R --Seurat $SEURAT`

#### Output

`cisTopic/cisTopic_all/plots/1{C,E}_ARC_seurat_40_Topics_featureplot_{TOPIC_NUM}_umap_rna.png`

#### Notes

This script was only run locally (not on the cluster).

---

### 3. cistopic_idr.R

#### Input
```
SEURAT=cisTopic/cisTopic_all/1C_ARC_seurat_cisTopicObject.Rds
SEURAT2=cisTopic/cisTopic_all/1E_ARC_seurat_cisTopicObject.Rds
NCORES=10
```

#### Script

```
Rscript ../../scripts/multiome/cistopic_idr.R --cisTopic1 $SEURAT --cisTopic2 $SEURAT2 \
--ncpus $N_CORES --idr_threshold 0.05
```

#### Output

```
CISTOPICDIR=cisTopic/cisTopic_all
$CISTOPICDIR/1{C,E}_ARC_seurat_cisTopicObject_topic_{TOPICNUM}_prob.bed
$CISTOPICDIR/IDR_files/1C_ARC_seurat_cisTopicObject_1E_ARC_seurat_cisTopicObject_{TOPIC_NUM_1C}_{TOPIC_NUM_1E}_IDR_0.05.bed
$CISTOPICDIR/1C_ARC_seurat_cisTopicObject_1E_ARC_seurat_cisTopicObject_IDR_0.05_size_list.csv
```

---

### 4. calculate_IDR_repr_score.R

#### Input

`CISTOPICDIR=cisTopic/cisTopic_all`
`SIZELIST=$CISTOPICDIR/IDR_files/1C_ARC_seurat_cisTopicObject_1E_ARC_seurat_cisTopicObject_IDR_0.05_size_list.csv`

#### Script

`Rscript ../../scripts/multiome/calculate_IDR_repr_score.R --sizelist $SIZELIST`

#### Output

`$CISTOPICDIR/IDR_files/1C_ARC_seurat_cisTopicObject_1E_ARC_seurat_cisTopicObject_IDR_0.05_size_list.csv`

---

### 5. calculate_heatmap_IDR_score_1C_1E.R

#### Input

```
CISTOPICDIR=cisTopic/cisTopic_all
CISTOPIC1=$CISTOPICDIR/1C_ARC_seurat_cisTopicObject.Rds
CISTOPIC2=$CISTOPICDIR/1E_ARC_seurat_cisTopicObject.Rds
META1=1C_ARC_metadata.csv
META2=1E_ARC_metadata.csv
SIZELIST=$CISTOPICDIR/IDR_files/1C_ARC_seurat_cisTopicObject_1E_ARC_seurat_cisTopicObject_IDR_0.05_size_list.csv
```

#### Script

```
Rscript ../../scripts/multiome/calculate_heatmap_IDR_score_1C_1E.R --cisTopic1 $CISTOPIC1 \
--cisTopic2 $CISTOPIC2 --meta1 $META1 --meta2 $META2 --idr_file $SIZELIST
```

#### Output
```
cisTopic/cisTopic_all/plots/heatmap_repr_score_1C_1E_full.png
cisTopic/cisTopic_all/plots/heatmap_repr_score_1C_1E.png
```

#### Notes
This script was only run locally (not on the cluster).

---

### 6. clean_corr_coverage_topics.R

#### Input

```
SIZELIST=$CISTOPICDIR/IDR_files/1C_ARC_seurat_cisTopicObject_1E_ARC_seurat_cisTopicObject_IDR_0.05_size_list.csv
TOPICS1=1,7,19,30,33,38,39,40
TOPICS2=17,21,25,30,36
```

#### Script

`../../scripts/multiome/clean_corr_coverage_topics.R --idr_file $SIZELIST --cisTopic1 $TOPICS1 --cisTopic2 $TOPICS2`

#### Output

`$CISTOPICDIR/IDR_files/1C_ARC_seurat_cisTopicObject_1E_ARC_seurat_cisTopicObject_IDR_0.05_size_list_CLEAN.csv`

---

### 7. cistopic_calculate_auc_IDR.R

#### Input

```
CISTOPIC1=cisTopic/cisTopic_all/1C_ARC_seurat_cisTopicObject.Rds
CISTOPIC2=cisTopic/cisTopic_all/1E_ARC_seurat_cisTopicObject.Rds
META1=1C_ARC_seurat_metadata.csv
META21E_ARC_seurat_metadata.csv
OUTFILE=cisTopic/cisTopic_all/1C_ARC_pca_rna_0.4_clusters_1E_ARC_pca_rna_0.6_clusters_cisTopic_zscores_AUC.csv
```

#### Script

```
Rscript ../../scripts/multiome/cistopic_calculate_auc_IDR.R --cisTopic1 $CISTOPIC1  \
--cisTopic2 $CISTOPIC2 --metadata1 $META1  --metadata2 $META2 --clusters1 pca_rna_0.4_clusters \
--clusters2 pca_rna_0.6_clusters --out $OUTFILE
```

#### Output

`cisTopic/cisTopic_all/1C_ARC_pca_rna_0.4_clusters_1E_ARC_pca_rna_0.6_clusters_cisTopic_zscores_AUC.csv`

---


### 8. merge_IDR_sizelist_cistopic_AUC.R

#### Input

```
SIZELIST=cisTopic/cisTopic_all/1C_ARC_seurat_cisTopicObject_1E_ARC_seurat_cisTopicObject_IDR_0.05_size_list_CLEAN.csv
AUC=cisTopic/cisTopic_all/1C_ARC_pca_rna_0.4_clusters_1E_ARC_pca_rna_0.6_clusters_cisTopic_zscores_AUC.csv
OUTFILE=cisTopic/cisTopic_all/1C_ARC_seurat_cisTopicObject_1E_ARC_seurat_cisTopicObject_IDR_0.05_size_list_CLEAN_clusters_AUC.csv
```

#### Script
```
Rscript ../../scripts/multiome/merge_IDR_sizelist_cistopic_AUC.R --sizelist $SIZELIST \
--auc $AUC --out $OUTFILE
```

#### Output

```
cisTopic/cisTopic_all/1C_ARC_seurat_cisTopicObject_1E_ARC_seurat_cisTopicObject_IDR_0.05_size_list_CLEAN_clusters_AUC.csv
```

---

### 9.

#### Input

```
SIGDIR=signatures

$SIGDIR/1C_ARC_seurat.Rds_pca_rna_0.4_clusters_wilcox_marker_genes_cluster6.tsv
$SIGDIR/1C_ARC_seurat.Rds_pca_rna_0.4_clusters_wilcox_marker_genes_cluster3.tsv
$SIGDIR/1C_ARC_seurat.Rds_pca_rna_0.4_clusters_wilcox_marker_genes_cluster4.tsv

$SIGDIR/1E_ARC_seurat.Rds_pca_rna_0.6_clusters_wilcox_marker_genes_cluster5.tsv
$SIGDIR/1E_ARC_seurat.Rds_pca_rna_0.6_clusters_wilcox_marker_genes_cluster3.tsv
$SIGDIR/1E_ARC_seurat.Rds_pca_rna_0.6_clusters_wilcox_marker_genes_cluster4.tsv
```

#### Script

```

Rscript ../../scripts/seurat_analysis/deg_signature_up.R 0 0.05 $SIGDIR/1C_1E_signatures_S1.txt $SIGDIR/1C_ARC_seurat.Rds_pca_rna_0.4_clusters_wilcox_marker_genes_cluster6.tsv $SIGDIR/1E_ARC_seurat.Rds_pca_rna_0.6_clusters_wilcox_marker_genes_cluster5.tsv

Rscript ../../scripts/seurat_analysis/deg_signature_up.R 0 0.05 $SIGDIR/1C_1E_signatures_S2.txt $SIGDIR/1C_ARC_seurat.Rds_pca_rna_0.4_clusters_wilcox_marker_genes_cluster3.tsv $SIGDIR/1E_ARC_seurat.Rds_pca_rna_0.6_clusters_wilcox_marker_genes_cluster3.tsv

Rscript ../../scripts/seurat_analysis/deg_signature_up.R 0 0.05 1C_1E_signatures_S3.txt 1C_ARC_seurat.Rds_pca_rna_0.4_clusters_wilcox_marker_genes_cluster4.tsv 1E_ARC_seurat.Rds_pca_rna_0.6_clusters_wilcox_marker_genes_cluster4.tsv

```

#### Output

```
$SIGDIR/1C_1E_signatures_S1.txt
$SIGDIR/1C_1E_signatures_S2.txt
$SIGDIR/1C_1E_signatures_S3.txt
```
---


### 10. calculate_signatures.R

#### Input

```
SIGDIR=signatures
CLUSTERS=$SIGDIR/1C_1E_signatures_S1.txt,$SIGDIR/1C_1E_signatures_S2.txt,$SIGDIR/1C_1E_signatures_S3.txt
OUTFILE=$SIGDIR/1C_1E_signatures_coords.txt
```

#### Script

`Rscript ../../scripts/multiome/calculate_signatures.R --signatures $CLUSTERS --names S1,S2,S3 --out $OUTFILE`

#### Output

`$SIGDIR/1C_1E_signatures_coords.txt`

---

### 11. calculate_marker_proximity_IDR.R

#### Input

```
CISTOPICDIR=cisTopic/cisTopic_all
IDR=$CISTOPICDIR/IDR_files/1C_ARC_seurat_cisTopicObject_1E_ARC_seurat_cisTopicObject_{TOPIC_NUM_1C}_{TOPIC_NUM_1E}_IDR_0.05.bed
SIGNATURES=signatures/1C_1E_signatures_coords.txt
OUTFILE=$CISTOPICDIR/IDR_files/1C_ARC_seurat_cisTopicObject_1E_ARC_seurat_cisTopicObject_{TOPIC_NUM_1C}_{TOPIC_NUM_1E}_IDR_0.05_with_markers.bed
```

#### Script 

`../../scripts/multiome/calculate_marker_proximity_IDR.R --idr $IDR --signatures $SIGNATURES --out $OUTFILE`

#### Output

`$CISTOPICDIR/IDR_files/1C_ARC_seurat_cisTopicObject_1E_ARC_seurat_cisTopicObject_{TOPIC_NUM_1C}_{TOPIC_NUM_1E}_IDR_0.05_with_markers.bed`

---

