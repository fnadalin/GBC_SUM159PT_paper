
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
    cat("\nUsage: <sc_exp_obj> <dir>\n")
    q()
}

SC_EXP_OBJ <- args[1]
DIR <- args[2]
   
library("Seurat")
library("ComplexHeatmap")
library("circlize")
library("scales")
library("ggrepel")
library("RColorBrewer")

cnv_ann <- read.table(file.path(DIR, "logfc_annotation.tsv"))
n_cnv <- nrow(cnv_ann)

object <- readRDS(SC_EXP_OBJ)

for (k in 1:n_cnv) {
    
    cnv_desc <- read.table(paste0(DIR,"/",k,"/cnv_desc.txt"), sep = "\t")[,1]
    features <- read.table(paste0(DIR,"/",k,"/feature_list.tsv"), sep = "\t", header = TRUE)[,1]
    idx <- match(features, rownames(object))
    norm_counts <- as.matrix(object@assays$RNA@data[idx,])

    NORM_COUNT_HEATMAP <- paste0(DIR,"/",k,"/heatmap_cnv_UMIcount.pdf")
    DENSITY_PLOT <- paste0(DIR,"/",k,"/density_cnv_UMIcount.pdf")

    M <- t(scale(t(norm_counts)))
    idx <- which(rowSums(is.na(M))==0)
    M <- M[idx,]

    v <- colSums(M)
    d <- density(v)
    
    pdf(DENSITY_PLOT)
    plot(d)
    dev.off()

    # select palette

    col_fun <- colorRamp2(c(-1,0,1), c("white","gray90","black"))

    h <- Heatmap(M, name = "Z-score", col = col_fun, column_title = cnv_desc, cluster_columns = TRUE,
                 cluster_rows = TRUE, show_column_names = FALSE, show_row_names = FALSE,
                 show_row_dend = FALSE, show_column_dend = FALSE, row_title = "gene loci spanned by CNV")
             
    pdf(NORM_COUNT_HEATMAP, width = 5, height = 4.5)
    draw(h)
    dev.off()

}


q()

