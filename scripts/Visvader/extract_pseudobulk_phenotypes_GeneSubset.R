
# select only genes from the input subset, plus add the same number of random genes

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) {
    cat("\nUsage: <object> <meta_id> <meta_val> <subset> <out_prefix>\n")
    q()
}

OBJECT <- args[1]
META_ID <- args[2]
META_VALUES <- unlist(strsplit(args[3], split = ","))
SUBSET <- args[4]
OUT_PREFIX <- args[5]

library("Seurat")

object <- readRDS(OBJECT)

genes <- read.table(SUBSET, sep = "\t")[,1]
genes <- genes[genes %in% rownames(object)]
idx <- which(rownames(object) %in% genes)
v <- 1:nrow(object)
rand_genes <- sample(v[-idx], size = length(genes))

features <- c(genes, rand_genes)
object <- subset(object, features = features)

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

