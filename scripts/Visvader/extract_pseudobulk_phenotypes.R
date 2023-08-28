
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
    cat("\nUsage: <object> <meta_id> <meta_val> <out_prefix>\n")
    q()
}

OBJECT <- args[1]
META_ID <- args[2]
META_VALUES <- unlist(strsplit(args[3], split = ","))
OUT_PREFIX <- args[4]

library("Seurat")

object <- readRDS(OBJECT)

samples <- unique(object@meta.data$sample.name)
phenotypes <- c(1,0)

phenotype <- rep(phenotypes[2], ncol(object))
idx <- which(object@meta.data[[META_ID]] %in% META_VALUES)
phenotype[idx] <- rep(phenotypes[1], length(idx))

M <- matrix(NA, nrow = nrow(object), ncol = length(samples)*length(phenotypes))
rownames(M) <- rownames(object)
colnames(M) <- unlist(lapply(samples, function(x) paste(x, phenotypes, sep = "_")))

j <- 1
for (s in samples) { 
    for (p in phenotypes) {
        idx <- which(object@meta.data$sample.name == s & phenotype == p)
        M[,j] <- apply(object@assays$RNA@data[,idx], 1, mean)
        j <- j+1
    }
}

OUT <- paste0(OUT_PREFIX, "_geneExprTable.tsv")
write.table(M, file = OUT, sep = "\t", quote = FALSE)

OUT <- paste0(OUT_PREFIX, "_phenotype.tsv")
df <- data.frame(phenotype = rep(phenotypes, length(samples)))
write.table(df, file = OUT, sep = "\t", quote = FALSE)


q()

