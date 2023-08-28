
# create a clone-sample matrix containing the number of cells per clone

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
    cat("\nUsage: <metadata.tsv> <out>\n")
    cat("\n<metadata.tsv>  where cell metadata are stored; must contain expr.GBC.list\n")
    cat("<out>           output table\n\n")
    q()
}

METADATA <- args[1]
OUT <- args[2]

meta <- read.table(METADATA, sep = "\t", header = TRUE)

gbc_list <- unique(meta$expr.GBC.list)
samples <- unique(meta$sample.name)

M <- matrix(0, nrow = length(gbc_list), ncol = length(samples))
rownames(M) <- gbc_list
colnames(M) <- samples

i <- 1
for (s in samples) {
    freq <- table(meta$expr.GBC.list[meta$sample.name == s])
    M[match(names(freq), gbc_list),i] <- freq
    i <- i+1
}

write.table(M, file = OUT, sep = "\t", quote = FALSE)

q()
