
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
    cat("\nUsage: <object> <sc_dataset>\n")
    q()
}

OBJECT <- args[1]
SC_DATASET <- args[2]

library("Scissor")
library("Seurat")

object <- readRDS(OBJECT)
data <- object@assays$RNA@counts
sc_dataset <- Seurat_preprocessing(data, verbose = FALSE)

saveRDS(sc_dataset, SC_DATASET) 

q()

