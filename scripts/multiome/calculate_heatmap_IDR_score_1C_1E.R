
############################## OPTIONS MENU ###################################

library("optparse")

# specify our desired options in a list
# by default OptionParser will add an help option equivalent to 
# make_option(c("-h", "--help"), action="store_true", default=FALSE, 
#               help="Show this help message and exit")

cat("\n")
option_list <- list(
  make_option("--cisTopic1", type="character", 
              help="[REQUIRED] cisTopic1 Object path"),
  make_option("--cisTopic2", type="character", 
              help="[REQUIRED] cisTopic1 Object path"),
  make_option("--idr_file", type="character",
              help="[REQUIRED] Path to the result of calculate_IDR_repr_score.R"),
  make_option("--meta1", type="character",
              help="[REQUIRED] Path to the metadata file of experiment 1"),
  make_option("--meta2", type="character",
              help="[REQUIRED] Path to the metadata file of experiment 2")
)

############################## PARSE OPTIONS ##################################


# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults
opt <- parse_args(OptionParser(option_list=option_list))

if (is.null(opt$cisTopic1)) {
  write("Option --cisTopic1 is required\nTry --help for help", stderr())
  q()
} else {
  CISTOPIC1 <- opt$cisTopic1
}

if (is.null(opt$cisTopic2)) {
  write("Option --cisTopic2 is required\nTry --help for help", stderr())
  q()
} else {
  CISTOPIC2 <- opt$cisTopic2
}

if (is.null(opt$meta2)) {
  write("Option --meta2 is required\nTry --help for help", stderr())
  q()
} else {
  META2 <- opt$meta2
}

if (is.null(opt$meta1)) {
  write("Option --meta1 is required\nTry --help for help", stderr())
  q()
} else {
  META1 <- opt$meta1
}

if (is.null(opt$idr_file)) {
  write("Option --idr_file is required\nTry --help for help", stderr())
  q()
} else {
  IDR.FILE <- opt$idr_file
}

################################# EXECUTION ###################################
suppressPackageStartupMessages({
  library(Signac)
  library(cisTopic)
  library(ComplexHeatmap)
  library(tidyverse)
  library(R.utils)
})

# Plot Topic-Topic heatmap showing Reproducibility score

# load list of Reproducibility scores
idr.list <- read.csv2(IDR.FILE, header = T, row.names = 1)

####################### Cleaned Heatmap #######################################

cisTopic1.blacklist <- c(1,7,16,19,30,33,38,39,40)
cisTopic2.blacklist <- c(17,20,21,25,30,36)

idr.list.clean <- subset(idr.list, !(cisTopic1 %in% cisTopic1.blacklist) & !(cisTopic2 %in% cisTopic2.blacklist))

order.1 <- paste("Topic", 
                 unique(idr.list.clean[order(idr.list.clean$repr.score, decreasing = T),"cisTopic1"]))
order.2 <- paste("Topic", 
                 unique(idr.list.clean[order(idr.list.clean$repr.score, decreasing = T),"cisTopic2"]))

idr_col_fun <- circlize::colorRamp2(breaks = c(0, 0.5, 1), colors = c("yellow", "orange", "red"))

# cast list into a matrix
idr.matrix <- reshape2::dcast(data = idr.list.clean, cisTopic1 ~ cisTopic2, value.var = "repr.score")

rownames(idr.matrix) <- paste("Topic", idr.matrix$cisTopic1)
idr.matrix$cisTopic1 <- NULL
colnames(idr.matrix) <- paste("Topic", colnames(idr.matrix))

idr.matrix[which(idr.matrix == 0, arr.ind = TRUE)] <- NA

png("heatmap_repr_score_1C_1E.png",
    height = 20,
    width = 25,
    units = "cm",
    res=100)
h <- Heatmap(as.matrix(idr.matrix[order.1, order.2]),
        cluster_rows = F, 
        cluster_columns = F,
        col = idr_col_fun,
        name="Reproducibility\nScore",
        row_names_gp = gpar(fontsize=18),
        column_names_gp = gpar(fontsize=18),
        row_title_gp = gpar(fontsize=18),
        column_title_gp = gpar(fontsize=18),
        column_title = "Replicate 2",
        row_title = "Replicate 1",
        heatmap_legend_param = list(
          title_gp = gpar(fontsize = 18,
                          fontface = "bold"), 
          labels_gp = gpar(fontsize = 18)
        )
        )
ComplexHeatmap::draw(h)
dev.off()

########################## Non-cleaned heatmap ################################
idr_col_fun <- circlize::colorRamp2(breaks = c(0, 0.5, 1), colors = c("yellow", "orange", "red"))

# cast list into a matrix
idr.matrix <- reshape2::dcast(data = idr.list, cisTopic1 ~ cisTopic2, value.var = "repr.score")

rownames(idr.matrix) <- paste("Topic", idr.matrix$cisTopic1)
idr.matrix$cisTopic1 <- NULL
idr.matrix[which(idr.matrix == 0, arr.ind = TRUE)] <- NA

colnames(idr.matrix) <- paste("Topic", colnames(idr.matrix))

order.1.idx <- unique(idr.list[order(idr.list$repr.score, decreasing = T),"cisTopic1"])
order.2.idx <- unique(idr.list[order(idr.list$repr.score, decreasing = T),"cisTopic2"])

order.1 <- paste("Topic", order.1.idx)
order.2 <- paste("Topic", order.2.idx)

meta.1C <- read.csv2(META1)
meta.1E <- read.csv2(META2)

# Find regions proximal to markers within 50kb
cisTopic1C <- readRDS(CISTOPIC1)
cisTopic1E <- readRDS(CISTOPIC2)

cisTopic1C <- selectModel(cisTopic1C,type = "maximum")
cisTopic1E <- selectModel(cisTopic1E,type = "maximum")

# calculate AUC between region-topic probability and marker proximality
cisTopic1C.emb <- modelMatSelection(cisTopic1C, target = "cell", method="Z-score", all.regions = TRUE)[order.1.idx,]
cisTopic1E.emb <- modelMatSelection(cisTopic1E, target = "cell", method="Z-score", all.regions = TRUE)[order.2.idx,]

corr_abs_ATAC_1C <- abs(cor(t(cisTopic1C.emb), meta.1C$nCount_ATAC)) # topics
corr_abs_ATAC_1E <- abs(cor(t(cisTopic1E.emb), meta.1E$nCount_ATAC)) # topics2
corr_abs_RNA_1C <- abs(cor(t(cisTopic1C.emb), meta.1C$nFeature_RNA)) # topics
corr_abs_RNA_1E <- abs(cor(t(cisTopic1E.emb), meta.1E$nFeature_RNA)) # topics2

cat(paste(length(corr_abs_ATAC_1C),length(corr_abs_ATAC_1E)))

pch_1C_ATAC <- rep("*",40)
pch_1C_ATAC[corr_abs_ATAC_1C < 0.5] <- NA
pch_1E_ATAC <- rep("*",40)
pch_1E_ATAC[corr_abs_ATAC_1E < 0.5] <- NA

pch_1C_RNA <- rep("*",40)
pch_1C_RNA[corr_abs_RNA_1C < 0.5] <- NA
pch_1E_RNA <- rep("*",40)
pch_1E_RNA[corr_abs_RNA_1E < 0.5] <- NA

col_fun_RNA <- circlize::colorRamp2(c(0, 0.5, 1), c("green", "white", "red"))
col_fun_ATAC <- circlize::colorRamp2(c(0, 0.5, 1), c("yellow", "white", "red"))

anno_1C <- rowAnnotation("ATAC Counts"=anno_simple(corr_abs_ATAC_1C, col=col_fun_ATAC, pch = pch_1C_ATAC),
                        "RNA Features"=anno_simple(corr_abs_RNA_1C, col=col_fun_RNA, pch= pch_1C_RNA))

anno_1E <- HeatmapAnnotation("ATAC Counts"=anno_simple(corr_abs_ATAC_1E, col=col_fun_ATAC, pch = pch_1E_ATAC),
                            "RNA Features"=anno_simple(corr_abs_RNA_1E, col=col_fun_RNA, pch= pch_1E_RNA))

lgd_corr_ATAC <- Legend(title = "Abs. correlation \nnCount_ATAC", col_fun = col_fun_ATAC, at=c(0,0.5,1))
lgd_corr_RNA <- Legend(title = "Abs. correlation \nnFeatures_RNA", col_fun = col_fun_RNA, at=c(0,0.5,1))
lgd_thresh <- Legend(pch = "*", type = "points", labels = "> 0.5")

heatmap.idr.full <- Heatmap(as.matrix(idr.matrix[order.1, order.2]),
        cluster_rows = F, 
        cluster_columns = F,
        col = idr_col_fun,
        name="Reproducibility\nScore",
        column_title = "Replicate 2",
        row_title = "Replicate 1",
        bottom_annotation = anno_1E,
        left_annotation =  anno_1C)

png("heatmap_repr_score_1C_1E_full.png",
    height = 20,
    width = 26,
    units = "cm",
    res=400)
ComplexHeatmap::draw(heatmap.idr.full, annotation_legend_list=list(lgd_corr_ATAC, lgd_corr_RNA, lgd_thresh))
dev.off()

########################## save datasets for reproducibility ##################
# 
# # If you want to re-run the algorithm, use this IDR list (complete and correct)
# write.csv2(idr.list, "hpc/cisTopic_all/1C_ARC_seurat_cisTopic_40_1E_ARC_seurat_cisTopic_40_cross_regions_overlap_IDR_AUC_list.csv",
#            quote = FALSE, row.names = FALSE)
# 
# write.csv2(idr.list.clean, "hpc/cisTopic_all/1C_ARC_seurat_cisTopic_40_1E_ARC_seurat_cisTopic_40_cross_regions_overlap_IDR_AUC_list_CLEAN.csv",
#            quote = FALSE, row.names = FALSE)
# 
# 
# idr.list <- read.csv2("hpc/cisTopic_all/1C_ARC_seurat_cisTopic_40_1E_ARC_seurat_cisTopic_40_cross_regions_overlap_IDR_AUC_list.csv")
# idr.list.clean <- read.csv2("hpc/cisTopic_all/1C_ARC_seurat_cisTopic_40_1E_ARC_seurat_cisTopic_40_cross_regions_overlap_IDR_AUC_list_CLEAN.csv")











