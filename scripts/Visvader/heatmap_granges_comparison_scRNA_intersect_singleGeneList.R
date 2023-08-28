
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
    cat("\nUsage: <sc_exp_obj> <gene_list> <out_dir>\n")
    q()
}

SC_EXP_OBJ <- args[1]
GENE_LIST <- args[2]
DIR <- args[3]
   
library("Seurat")
library("ComplexHeatmap")
library("circlize")
library("scales")
library("ggrepel")
library("RColorBrewer")

dir.create(DIR, showWarnings = FALSE, recursive = TRUE)

object <- readRDS(SC_EXP_OBJ)

features <- read.table(GENE_LIST, sep = "\t", header = TRUE)[,1]
features <- features[features %in% rownames(object)]
idx <- match(features, rownames(object))
norm_counts <- as.matrix(object@assays$RNA@data[idx,])

NORM_COUNT_HEATMAP <- file.path(DIR,"heatmap_cnv_UMIcount.pdf")
DENSITY_PLOT <- file.path(DIR,"density_cnv_UMIcount.pdf")

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

h <- Heatmap(M, name = "Z-score", col = col_fun, cluster_columns = TRUE, 
             cluster_rows = TRUE, show_column_names = FALSE, show_row_names = FALSE,
             show_row_dend = FALSE, show_column_dend = FALSE)
             
pdf(NORM_COUNT_HEATMAP, width = 5, height = 4.5)
draw(h)
dev.off()



q()

