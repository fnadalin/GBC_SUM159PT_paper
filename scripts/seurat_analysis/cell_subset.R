
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
    cat("\nUsage: <object> <mode> <cl> <cells>\n")
    cat("\n<object>        input object, containing <mode> metadata\n")
    cat("<mode>          cluster mode as stored in metadata\n")
    cat("<cl>            label of the cluster(s) to be filtered out\n")
    cat("<cells>         output list of cells retained\n\n")
    q()
}

OBJ <- args[1]
MODE <- args[2]
CL <- args[3]
CELLS <- args[4]

library("Seurat")

cls <- unlist(strsplit(CL, split=","))
dr <- gsub("_k.*", "", gsub("clusters_", "umap_", MODE))

object <- readRDS(OBJ)

dir <- dirname(CELLS)
pdf(file.path(dir, "umap_orig_UMIcount.pdf"))
FeaturePlot(object, reduction = dr, features = "nCount_RNA")
dev.off()

idx <- which(!(object@meta.data[[MODE]] %in% cls))
cells <- colnames(object)[idx]

write.table(cells, file = CELLS, col.names = FALSE, row.names = FALSE, quote = FALSE)


q()

