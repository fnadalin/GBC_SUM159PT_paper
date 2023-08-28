
MAX_UMI <- 20000 

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
    cat("\nUsage: <in_obj> <out_dir>\n")
    q()
}

IN_OBJ <- args[1]
OUT_DIR <- args[2]

library("Seurat")

object <- readRDS(IN_OBJ)

norm_exp_epcam <- object@assays$RNA@data["EPCAM",]
d <- density(norm_exp_epcam)

which.2nd_max <- which.max(d$y[d$x > 0.5]) + head(which(d$x > 0.5),1)
which.2nd_min <- which.min(d$y[d$x < which.2nd_max])
which.2nd_min <- which.min(d$y[d$x < d$x[which.2nd_max]])

pdf(file.path(OUT_DIR, "density_EPCAM.pdf"))
plot(d)
abline(v = d$x[which.2nd_min], col = "magenta")
text(x = d$x[which.2nd_min], y = -0.2, label = round(d$x[which.2nd_min], digits = 2), col = "magenta")
dev.off()

sum(norm_exp_epcam > d$x[which.2nd_min])
sum(norm_exp_epcam < d$x[which.2nd_min])

is_epcam_pos <- rep(0, ncol(object))
idx <- which(norm_exp_epcam > d$x[which.2nd_min])
is_epcam_pos[idx] <- rep(1, length(idx))

df <- data.frame(is_epcam_pos = is_epcam_pos)
rownames(df) <- colnames(object)
object <- AddMetaData(object, metadata = df)

cells <- colnames(object)[object@meta.data$is_epcam_pos == 1]
write.table(cells, file = file.path(OUT_DIR, "cells_epcam_pos.tsv"), quote = FALSE, col.names = FALSE, row.names = FALSE)

object <- FindVariableFeatures(object)
object <- RunPCA(object)
object <- RunUMAP(object, reduction = "pca", dims = 1:10)

saveRDS(object, file = IN_OBJ)

pdf(file.path(OUT_DIR, "umap_all_is_EPCAM.pdf"))
DimPlot(object, dims = 1:2, reduction = "umap", cells.highlight = cells)
dev.off()

pdf(file.path(OUT_DIR, "umap_all_EPCAM.pdf"))
FeaturePlot(object, dims = 1:2, reduction = "umap",  features = "EPCAM")
dev.off()

object <- subset(object, cells = cells)
object <- subset(object, subset = nCount_RNA < MAX_UMI)

saveRDS(object, file = file.path(OUT_DIR, "object.Rds"))

q()

