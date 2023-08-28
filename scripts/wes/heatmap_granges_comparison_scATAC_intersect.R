
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
    cat("\nUsage: <sc_name>\n")
    q()
}

SC_EXP_NAME <- args[1]

library("ComplexHeatmap")
library("circlize")
library("scales")
library("RColorBrewer")

OUTDIR <- paste0("out_", SC_EXP_NAME, "-WES_intersect")
    
NORM_COUNT_TABLE <- file.path(OUTDIR, "matrix.tsv")
COL_ANNOTATION <- file.path(OUTDIR, "cl_annotation.txt")
ROW_ANNOTATION <- file.path(OUTDIR, "logfc_annotation.tsv")
    
AUC_TABLE <- file.path(OUTDIR, "auc.txt")
    
NORM_COUNT_HEATMAP <- file.path(OUTDIR, "heatmap_norm_counts.pdf")
AUC_HEATMAP <- file.path(OUTDIR, "heatmap_auc.pdf")
    
# load data
    
auc <- as.matrix(read.table(AUC_TABLE, sep = "\t", header = TRUE))
M <- as.matrix(read.table(NORM_COUNT_TABLE, sep = "\t", header = TRUE))
colnames(auc) <- gsub("X","",colnames(auc))
row_annotation <- as.matrix(read.table(ROW_ANNOTATION))
colnames(row_annotation) <- paste0("WES",1:ncol(row_annotation))
col_annotation <- read.table(COL_ANNOTATION)[,1]
    
ordering <- order(col_annotation)
col_annotation <- col_annotation[ordering]
idx <- order(apply(row_annotation,1,mean),decreasing = TRUE)
M <- M[idx,ordering]
row_annotation <- row_annotation[idx,]
auc <- auc[idx,]

M <- t(scale(t(M)))

# generate AUC heatmap

col_fun_ann <- colorRamp2(c(min(row_annotation),0,max(row_annotation)), c("blue","gray90","red"))
col_fun <- colorRamp2(c(min(auc),0.5,max(auc)), c("white","gray90","gray10"))  

ha <- HeatmapAnnotation(log2fc = row_annotation, which = "row", col = list(log2fc = col_fun_ann), 
                       show_legend = TRUE, show_annotation_name = TRUE)
h <- Heatmap(auc, name = "AUC", col = col_fun, right_annotation = ha, cluster_rows = FALSE, cluster_columns = FALSE,
             column_names_side = "top", column_title = "Multiome clusters at day 0")
    
pdf(AUC_HEATMAP, width = 6, height = 5)
draw(h)
dev.off()
    
# generate norm count heatmap
  
cl <- unique(col_annotation)
ncl <- length(cl)
   
col_fun_cl <- hue_pal()(ncl)
names(col_fun_cl) <- cl
col_fun <- colorRamp2(c(-2,0,2), c("white","gray80","black"))
    
# NB: do not put DTC annotation because we don't see anything
    
row_split <- rep("depleted p.t.", nrow(row_annotation))
idx <- which(apply(row_annotation,1,mean) > 0)
row_split[idx] <- rep("amplified p.t.", length(idx))
   
labelfont <- rep("plain",nrow(M))
labelfont[1] <- "bold" 

hb <- HeatmapAnnotation(cluster = col_annotation, which = "column", col = list(cluster = col_fun_cl),
                        show_legend = TRUE, show_annotation_name = FALSE)
h <- Heatmap(M, name = "Z-score", col = col_fun, left_annotation = ha, top_annotation = hb,
             column_split = col_annotation, row_split = row_split, row_names_gp = gpar(fontface = labelfont),
             cluster_rows = FALSE, show_column_names = FALSE, column_title = "DNA accessibility at day 0",
             row_title = "reproducible CNV regions")
    
pdf(NORM_COUNT_HEATMAP, width = 7.5, height = 5)
draw(h)
dev.off()
    


q()

