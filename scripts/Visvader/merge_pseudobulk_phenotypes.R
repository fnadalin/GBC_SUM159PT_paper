
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
    cat("\nUsage: <in_prefixes> <out_prefix>\n")
    q()
}

IN_PREFIXES <- unlist(strsplit(args[1], split = ","))
OUT_PREFIX <- args[2]

nsamples <- 0
genes <- NULL
for (in_prefix in IN_PREFIXES) {

    filename <- paste0(in_prefix, "_geneExprTable.tsv")
    table <- read.table(filename, sep = "\t", header = TRUE)
    g <- rownames(table)
    if (length(genes) == 0) {
        genes <- g
    } else {
        genes <- genes[genes %in% g]
    }
    nsamples <- nsamples + ncol(table)
    
}

df_exp <- as.data.frame(matrix(NA, nrow = length(genes), ncol = 0))
df_phe <- as.data.frame(matrix(NA, nrow = nsamples, ncol = 0))
for (in_prefix in IN_PREFIXES) {

    filename <- paste0(in_prefix, "_geneExprTable.tsv")
    table <- read.table(filename, sep = "\t", header = TRUE)
    g <- rownames(table)
    idx <- match(genes, g)
    table <- table[idx,]
    df_exp <- cbind(df_exp, table)
    
    filename <- paste0(in_prefix, "_phenotype.tsv")
    phenotype <- read.table(filename, "\t", header = TRUE)
    df_phe <- rbind(df_phe, phenotype)
}
rownames(df_phe) <- 1:nrow(df_phe)

OUT <- paste0(OUT_PREFIX, "_geneExprTable.tsv")
write.table(df_exp, file = OUT, sep = "\t", quote = FALSE)

OUT <- paste0(OUT_PREFIX, "_phenotype.tsv")
write.table(df_phe, file = OUT, sep = "\t", quote = FALSE)


q()

