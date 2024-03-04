# Take a Seurat object and perform cisTopic analysis on it, 
# extract the cell-topics matrix and attaches it to the original seurat object as a reduction
# the full cistopic object is stored in the misc slot of the reduction of the seurat object

set.seed(1234)

############################## OPTIONS MENU ###################################

library("optparse")

# specify our desired options in a list
# by default OptionParser will add an help option equivalent to 
# make_option(c("-h", "--help"), action="store_true", default=FALSE, 
#               help="Show this help message and exit")

cat("\n")
option_list <- list(
  make_option("--Seurat1", type="character", 
              help="[REQUIRED] Seurat Object 1 path"),
  make_option("--Seurat2", type="character", 
              help="[REQUIRED] Seurat Object 2 path"),
  make_option("--clustering_name1", type="character",
              help="[REQUIRED] name of the metadata column of the clusters"),
  make_option("--clustering_name2", type="character",
              help="[REQUIRED] name of the metadata column of the clusters"))

############################## PARSE OPTIONS ##################################


# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults
opt <- parse_args(OptionParser(option_list=option_list))

if (is.null(opt$Seurat1)) {
  write("Option --Seurat1 is required\nTry --help for help", stderr())
  q()
} else {
  SEURAT1 <- opt$Seurat1
}


if (is.null(opt$Seurat2)) {
  write("Option --Seurat2 is required\nTry --help for help", stderr())
  q()
} else {
  SEURAT2 <- opt$Seurat2
}

if (is.null(opt$clustering_name1)) {
  write("Option --clustering_name1 is required\nTry --help for help", stderr())
  q()
} else {
  CLUST_NAME1 <- opt$clustering_name1
}

if (is.null(opt$clustering_name2)) {
  write("Option --clustering_name2 is required\nTry --help for help", stderr())
  q()
} else {
  CLUST_NAME2 <- opt$clustering_name2
}

################################# EXECUTION ###################################

library(Seurat)
library(Signac)
library(ggplot2)
library(dplyr)
library(ComplexHeatmap)
library(cisTopic)
library(EnsDb.Hsapiens.v86)

# we load the Seurat object of our experiment
object1 <- readRDS(SEURAT1)

seurat1.name <- sub("_seurat.Rds","",basename(SEURAT1))
cisTopicObject1 <- object1@reductions$cisTopic@misc$cisTopicObject

cor.cells.1 <- cor(object1@reductions$cisTopic@cell.embeddings)

png(file.path(dirname(SEURAT1), paste(seurat1.name, "topic_correlation_matrix_cellbased.png", sep="_")), 
    width=16, 
    height=15,
    units="cm",
    res=300)
Heatmap(cor.cells.1,
        column_title=paste(seurat1.name, "Topic Correlation Matrix (cell-based)"),
        row_names_gp = gpar(fontsize=10),
        column_names_gp = gpar(fontsize=10),heatmap_legend_param = list("title"="Pearson's R")
        )
dev.off()

topic.regions.1 <- modelMatSelection(cisTopicObject1, 
                                     method="Z-score", 
                                     target = "region", 
                                     all.regions = TRUE)

cor.regions.1 <- cor(t(topic.regions.1))

png(file.path(dirname(SEURAT1),paste(seurat1.name, "topic_correlation_matrix_regionsbased.png", sep="_")), 
    width=16, 
    height=15, 
    units="cm",
    res=300)
Heatmap(cor.regions.1,
        column_title=paste(seurat1.name, "Topic Correlation Matrix (region-based)"),
        row_names_gp = gpar(fontsize=10),
        column_names_gp = gpar(fontsize=10),heatmap_legend_param = list("title"="Pearson's R")
)
dev.off()

object2 <- readRDS(SEURAT2)

seurat2.name <- sub("_seurat.Rds","",basename(SEURAT2))
cisTopicObject2 <- object2@reductions$cisTopic@misc$cisTopicObject

cor.cells.2 <- cor(object2@reductions$cisTopic@cell.embeddings)

png(file.path(dirname(SEURAT2),paste(seurat2.name, "topic_correlation_matrix_cellbased.png", sep="_")), 
    width=16, 
    height=15,
    units="cm",
    res=300)
Heatmap(cor.cells.2,
        column_title=paste(seurat2.name, "Topic Correlation Matrix (cell-based)"),
        row_names_gp = gpar(fontsize=10),
        column_names_gp = gpar(fontsize=10),heatmap_legend_param = list("title"="Pearson's R")
)
dev.off()

topic.regions.2 <- modelMatSelection(cisTopicObject2, 
                                     method="Z-score", 
                                     target = "region", 
                                     all.regions = TRUE)

cor.regions.2 <- cor(t(topic.regions.2))

png(file.path(dirname(SEURAT2),paste(seurat2.name, "topic_correlation_matrix_regionsbased.png", sep="_")), 
    width=16, 
    height=15, 
    units="cm",
    res=300)
Heatmap(cor.regions.2,
        column_title=paste(seurat2.name, "Topic Correlation Matrix (region-based)"),
        row_names_gp = gpar(fontsize=10),
        column_names_gp = gpar(fontsize=10),heatmap_legend_param = list("title"="Pearson's R")
)
dev.off()

############################# EXPERIMENT COMPARISONS ##########################

overlaps <- findOverlaps(query=cisTopicObject1@region.ranges,subject=cisTopicObject2@region.ranges)

topic.regions.1 <- t(topic.regions.1)[queryHits(overlaps),]
topic.regions.2 <- t(topic.regions.2)[subjectHits(overlaps),]

cross.cor.regions <- cor(topic.regions.1, topic.regions.2)

p1 <- Heatmap(cross.cor.regions,
        row_names_gp = gpar(fontsize=10),
        row_title_side = "left",
        row_title = seurat2.name,
        column_names_gp = gpar(fontsize=10),
        column_title_side = "bottom",
        column_title=seurat1.name,
        heatmap_legend_param = list("title"="Pearson's R"))

png(file.path(dirname(SEURAT1),paste(seurat1.name, seurat2.name, "topic_cross_correlation_matrix_regionsbased.png", sep="_")), 
        width=16, 
        height=15, 
        units="cm",
        res=300)
draw(p1, column_title=paste(seurat1.name, seurat2.name, "Topic Correlation Matrix (region-based)"))
dev.off()

