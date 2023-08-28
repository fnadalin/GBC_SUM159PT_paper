
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
    cat("\nUsage: <sc_dataset> <pops> <out_DIR>\n")
    q()
}

OBJECT <- args[1]
POPS <- unlist(strsplit(args[2], split = ","))
OUT_DIR <- args[3]

library("Seurat")

object <- readRDS(OBJECT)

for (s in POPS) {
    meta <- paste0("scissor.", s)
    x <- as.numeric(object@meta.data[[meta]] == 1)
    colname <- paste0(meta, ".binary")
    object <- AddMetaData(object, metadata = x, col.name = colname)
    DEG <- FindMarkers(object, group.by = colname, ident.1 = 1)
    write.table(DEG, file = paste0(OUT_DIR, "/DEA_", s, ".tsv"), sep = "\t", quote = FALSE)
}


q()

