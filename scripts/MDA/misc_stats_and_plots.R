
S1_EMT_III <- c("S100A6","HSPB1","TM4SF1","S100A4","LGALS3","CLIC1","ANXA1","KRT18")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
    cat("\nUsage: <DIR> <slot> <sign_file> <sign_name>\n")
    q()
}

DIR <- args[1]
SLOT <- args[2]
SIGN_FILE <- args[3]
SIGN_NAME <- args[4]

library("Seurat")
library("ComplexHeatmap")
library("ggrepel")

PCA_SLOT <- gsub("clusters_", "", gsub("_k.+", "", SLOT))
UMAP_SLOT <- paste0("umap_", PCA_SLOT)

OUTDIR <- file.path(DIR, "figures_signatures", SLOT)
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

sign <- read.table(SIGN_FILE, sep = "\t", header = TRUE)[,1]

object <- readRDS(file.path(DIR,"hvg_pca_clust/object.Rds"))
cl <- sort(unique(object@meta.data[[SLOT]]))

top_50_features <- c()
for (i in cl) {
    filename <- paste0(file.path(DIR, "DEA/"), SLOT, "/MAST/DEG_MAST_cl", i, "-all.tsv")
    df <- read.table(filename, sep = "\t", header = TRUE)
    df_filt <- df[df$avg_log2FC > 0,]
    comm <- sign[sign %in% df_filt$geneID]
    cat(paste0("cluster ", i, " ", SIGN_NAME, "+NQO1 signature genes: ", paste(comm, collapse = ","), "\n"))
    rank1 <- match(comm, sign)
    rank2 <- match(comm, df_filt$geneID)
    cat(paste0("rank of common genes in ", SIGN_NAME, " +NQO1 signature: ", paste(rank1, collapse = ","), "\n"))
    cat(paste0("rank of common genes in cluster ", i, ": ", paste(rank2, collapse = ","), "\n"))
    top_50_features <- c(top_50_features, df_filt$geneID[1:50])
}

M <- AverageExpression(object, features = c(sign[1:50], "NQO1"), group.by = SLOT)
M <- t(scale(t(M$RNA)))

pdf(file.path(OUTDIR, paste0("heatmap_", SIGN_NAME, "_features.pdf")), width = 4, height = 8)
Heatmap(M, cluster_rows = FALSE, cluster_columns = FALSE)
dev.off()

pdf(file.path(OUTDIR, paste0("umap_", SIGN_NAME, "_features.pdf")), width = 15, height = 12)
FeaturePlot(object, features = c(sign[1:20], "NQO1"), reduction = UMAP_SLOT, order = TRUE, ncol = 5, cols = c("yellow","blue"))
dev.off()

pdf(file.path(OUTDIR, paste0("vln_", SIGN_NAME, "_features.pdf")), width = 15, height = 12)
VlnPlot(object, features = c(sign[1:20], "NQO1"), group.by = SLOT, ncol = 5, pt.size = 0)
dev.off()

M <- AverageExpression(object, features = top_50_features, group.by = SLOT)
M <- t(scale(t(M$RNA)))

M <- M[rownames(M) %in% top_50_features,]
top_50_features <- top_50_features[top_50_features %in% rownames(M)]
M <- M[match(top_50_features, rownames(M)),] # order by significance 
comm <- rownames(M)[rownames(M) %in% sign]
idx_val <- match(comm, rownames(M))

gene_labels_row <- rowAnnotation(highlight = anno_mark(at = idx_val, labels = rownames(M)[idx_val], 
                                 side = 'right', labels_gp = gpar(fontsize = 10)))

pdf(file.path(OUTDIR, "heatmap_top50clusterFeatures.pdf"), width = 4, height = 8)
Heatmap(M, cluster_rows = FALSE, cluster_columns = FALSE, name = "Z-score", show_row_names = FALSE,
        column_title = "Top 50 cluster markers", right_annotation = gene_labels_row)
dev.off()


if (SIGN_NAME != "S1") { q() }

M <- AverageExpression(object, features = c(S1_EMT_III, "NQO1"), group.by = SLOT)
M <- t(scale(t(M$RNA)))

pdf(file.path(OUTDIR, "heatmap_S1_EMT_III_features.pdf"), width = 4, height = 4)
Heatmap(M, cluster_rows = FALSE, cluster_columns = FALSE)
dev.off()

pdf(file.path(OUTDIR, "umap_S1_EMT_III_features.pdf"), width = 15, height = 6)
FeaturePlot(object, features = c(S1_EMT_III, "NQO1"), reduction = UMAP_SLOT, order = TRUE, ncol = 5, cols = c("yellow","blue"))
dev.off()

pdf(file.path(OUTDIR, "vln_S1_EMT_III_features.pdf"), width = 15, height = 6)
VlnPlot(object, features = c(S1_EMT_III, "NQO1"), group.by = SLOT, ncol = 5, pt.size = 0)
dev.off()


q()

