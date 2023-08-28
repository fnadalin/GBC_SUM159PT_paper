args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
    cat("\nUsage: <object> <lin>\n")
    cat("\n<object>  Seurat object\n")
    cat("<lin>     tsv file with clone ID and lineage ID\n")
    q()
}

OBJECT <- args[1]
LIN <- args[2]

library("Seurat")

object <- readRDS(OBJECT)
lineages <- read.table(LIN, sep = "\t")
classes <- unique(lineages[,2])

lineage_metadata <- rep(NA, ncol(object))
for (c in classes) {
    idx <- which(object@meta.data$expr.GBC.list %in% lineages[lineages[,2] == c,1])
    lineage_metadata[idx] <- rep(c, length(idx))
}

object <- AddMetaData(object, metadata = lineage_metadata, col.name = "lineage")

saveRDS(object, file = OBJECT)

q()

