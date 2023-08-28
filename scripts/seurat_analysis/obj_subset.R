
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
    cat("\nUsage: <in_object> <cells> <out_object>\n")
    cat("\nsubset the object on <cells>\n\n")
    q()
}

IN_OBJ <- args[1]
CELLS <- args[2]
OUT_OBJ <- args[3]

library("Seurat")

object <- readRDS(IN_OBJ)
cells <- read.table(CELLS)[,1]

object <- subset(object, cells = cells)

saveRDS(object, file = OUT_OBJ)


q()

