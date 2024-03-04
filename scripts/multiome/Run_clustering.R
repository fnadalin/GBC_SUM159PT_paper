
# Takes as input a Seurat object and perform a clustering of the data (or a reduction of the data itself)

# what options do we want to choose?

# which view/modality/reduction to cluster on
# if and which features should we use to cluster (variable features, set of genebasis genes, etc)
# what algorithms do we want to test
# what resolution do we want to use

set.seed(1234)
ALGORITHM <- 3 # SLM algorithm is set as default
RESOLUTION <- 0.5 # clustering resolution is set to 0.5
############################## OPTIONS MENU ###################################

library("optparse")

# specify our desired options in a list
# by default OptionParser will add an help option equivalent to 
# make_option(c("-h", "--help"), action="store_true", default=FALSE, 
#               help="Show this help message and exit")

cat("\n")
option_list <- list(
  make_option("--Seurat", type="character", 
              help="[REQUIRED] Seurat Object path"),
  make_option("--assay", type="character", default="RNA",
              help="[OPTIONAL] Assay from which select variable features. Only active when --features is active [default=%default]"),
  make_option("--out_RDS", type="character",
              help="[OPTIONAL] Output file name, otherwise input file name will be used"),
  make_option("--reduction", type="character", default="pca", 
              help="[OPTIONAL] Reduction to be used to construct the SNN graph with FindNeighbours function (PCA, MOFA, cisTopic ecc). Set to NULL if --features is used [default=%default]"),
  make_option("--algorithm", type="integer", default=3,
              help="[OPTIONAL] Value for the 'algorithm' parameter in the FindClusters Seurat function [default=%default]"),
  make_option("--resolution", type="character", default="0.5",
              help="[OPTIONAL] Value for the 'resolution' parameter in the FindClusters Seurat function [default=%default]"),
  make_option("--features", type="character", 
              help="Features to select for clustering. Only works if reduction is NULL [choose between geneBasis or VariableFeatures]"))


############################# PARSE OPTIONS ###################################

opt <- parse_args(OptionParser(option_list=option_list))

if (is.null(opt$Seurat)){
  write("Option --Seurat is required\nTry --help for help", stderr())
  q()
} else {
  SEURAT <- opt$Seurat
}

if (is.null(opt$out_RDS)) {
  OUT_FILE <- opt$Seurat
}

if (!is.null(opt$reduction)) {
  REDUCTION <- opt$reduction
  FEATURES <- NULL
}

if (!is.null(opt$algorithm)) {
  ALGORITHM <- as.integer(opt$algorithm)
}

if (!is.null(opt$resolution)) {
  RESOLUTION <- as.numeric(unlist(strsplit(opt$resolution, ",")))
}

if (!is.null(opt$features)) {
  FEATURES <- opt$features
}

if (!is.null(opt$assay)) {
  ASSAY <- opt$assay
}

############################### EXECUTION #####################################
library(Seurat)
library(Signac)
library(ggplot2)

EXP_NAME <- sub(".rds", "", basename(SEURAT), ignore.case = T)
object <- readRDS(SEURAT)

# check if the object has no silhouette compartment in the misc section
if (is.null(object@misc$silhouette)) {
  Misc(object, slot = "silhouette") <- list()
}

if (!is.null(FEATURES)) {
  # if --feature is active, that means that we're using the VariableFeatures or the geneBasis gene list
  FEATURES <- switch(EXPR = FEATURES,
                     VariableFeatures=VariableFeatures(object[[ASSAY]]),
                     geneBasis=object@misc$geneBasis$gene)
  
  object <- FindNeighbors(object, features = FEATURES)
} else {
  # regardless of the reduction, we include all dimensions
  DIMS <- seq(dim(object@reductions[[REDUCTION]]@cell.embeddings)[2])
  if (REDUCTION == "lsi_atac") DIMS <- seq(2,dim(object@reductions[[REDUCTION]]@cell.embeddings)[2])
  GRAPH_NAME_nn <- paste(REDUCTION, "nn", sep='_')
  GRAPH_NAME_snn <- paste(REDUCTION, "snn", sep='_')
  object <- FindNeighbors(object, reduction = REDUCTION, dims = DIMS, graph.name = c(GRAPH_NAME_nn,GRAPH_NAME_snn))
  
  # distance matrix to calculate the Silhouette plot
  reduction.dist <- dist(object@reductions[[REDUCTION]]@cell.embeddings)
}
dir.create(paste(dirname(SEURAT),"/plots", sep=''), recursive=T)

if (length(RESOLUTION) > 1) {
  for (res in RESOLUTION) {
    
    CLUST_NAME <- paste(REDUCTION,res,"clusters",sep="_")
    clus <- FindClusters(object, algorithm = ALGORITHM, resolution = res, graph.name = GRAPH_NAME_snn)
    object@meta.data[[CLUST_NAME]] <- clus$seurat_clusters
    
    # RNA plot
    p1 <- DimPlot(object, reduction="umap_rna", group.by = CLUST_NAME, label = TRUE) + NoLegend() + ggtitle(paste("RNA UMAP - Clusters from ", REDUCTION, "- Res.", res))
    # ATAC plot
    p2 <- DimPlot(object, reduction="umap_atac", group.by = CLUST_NAME, label = TRUE) + ggtitle("ATAC UMAP")
    
    p1+p2
    ggsave(paste(paste(dirname(SEURAT),"/plots/", sep=''), EXP_NAME, "_clustering_",CLUST_NAME,".png", sep=''), width=15,height=5)
    
    # calculating silhouette plot
    silhouette <- cluster::silhouette(as.numeric(clus$seurat_clusters), reduction.dist)
    
    Misc(object, slot="silhouette")[[CLUST_NAME]] <- silhouette
    
    png(paste(paste(dirname(SEURAT),"/plots/", sep=''), EXP_NAME, "_clustering_",CLUST_NAME,"_silhouette_plot", ".png", sep=''), 
           width=15,height=10, units = "cm", res=300)
    plot(cluster::sortSilhouette(silhouette), 
         border=NA, 
         col="lightsalmon",
         main=paste("Silhouette plot for", REDUCTION, "- Res. ", res))
    dev.off()
  }
} else {
  CLUST_NAME <- paste(REDUCTION,RESOLUTION,"clusters",sep="_")
  clus <- FindClusters(object, algorithm = ALGORITHM, resolution = RESOLUTION, graph.name = GRAPH_NAME_snn)
  object@meta.data[[CLUST_NAME]] <- clus$seurat_clusters
  
  # RNA plot
  p1 <- DimPlot(object, reduction="umap_rna", group.by = CLUST_NAME, label = TRUE) + NoLegend() + ggtitle(paste("RNA UMAP - Clusters from ", REDUCTION, "- Res.", RESOLUTION))
  # ATAC plot
  p2 <- DimPlot(object, reduction="umap_atac", group.by = CLUST_NAME, label = TRUE) + ggtitle("ATAC UMAP")
  
  p1+p2
  ggsave(paste(paste(dirname(SEURAT),"/plots/", sep=''), EXP_NAME, "_clustering_",CLUST_NAME,".png", sep=''), width=15,height=5)
  
  # calculating silhouette plot
  silhouette <- cluster::silhouette(as.numeric(clus$seurat_clusters), reduction.dist)
  
  Misc(object, slot="silhouette")[[CLUST_NAME]] <- silhouette
  
  png(paste(paste(dirname(SEURAT),"/plots/", sep=''), EXP_NAME, "_clustering_",CLUST_NAME,"_silhouette_plot", ".png", sep=''), 
         width=10,height=15, units = "cm", res=300)
  plot(cluster::sortSilhouette(silhouette), 
       border=NA, 
       col="lightsalmon",
       main=paste("Silhouette plot for", REDUCTION, "- Res. ", res))
  dev.off()
}

saveRDS(object, file = OUT_FILE)
q()