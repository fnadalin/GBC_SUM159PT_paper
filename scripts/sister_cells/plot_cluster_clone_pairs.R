
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
        cat("\nUsage: Rscript plot_cluster_clone_pairs.R <out_dir>\n")
        q()
}

OUTDIR <- args[1]

library("ComplexHeatmap")
library("RColorBrewer")
library("circlize")
library("scales")

M_norm <- read.table(file.path(OUTDIR, "cluster_pairs_norm.tsv"), sep = "\t")
rownames(M_norm) <- 1:nrow(M_norm)-1
colnames(M_norm) <- 1:ncol(M_norm)-1

# plot matrix

m1 <- min(M_norm)
m2 <- max(M_norm)
col_fun <- colorRamp2(c(m1,m2), c("white", "red"))

col_row <- hue_pal()(nrow(M_norm))
col_col <- hue_pal()(ncol(M_norm))
names(col_row) <- rownames(M_norm)
names(col_col) <- colnames(M_norm)

ha <- HeatmapAnnotation(cluster = 1:nrow(M_norm)-1, which = "row", col = list(cluster = col_row))
hb <- HeatmapAnnotation(cluster = 1:ncol(M_norm)-1, col = list(cluster = col_col))
h <- Heatmap(M_norm, name = "normalized pairs count", col = col_fun, left_annotation = ha, top_annotation = hb,
		cluster_rows = FALSE, cluster_columns = FALSE, row_names_side = "left", column_names_side = "top")

pdf(file.path(OUTDIR,"heatmap_cluster_pairs_norm.pdf"), width = 5, height = 3)
draw(h)
dev.off()


q()


