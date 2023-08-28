
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) {
    cat("\nUsage: <min_log2fc> <max_pval> <out> <dea1> <dea2> [<dea3> ... <deaN>] \n")
    cat("\n<min_log2fc> log fold-change cutoff\n")
    cat("<max_pval>   adj. p-value cutoff\n")
    cat("<out>        output file containing the signature\n")
    cat("<deaI>       input tsv file containing DEGs (avg_log2FC and geneID fields at least)\n\n")
    q()
}

MAX_LOG2FC <- as.numeric(args[1])
MAX_PVAL <- as.numeric(args[2])
SIGNATURE <- args[3]

genes <- c()

for (i in 4:length(args)) {

    DEA <- args[i]
    deg_table <- read.table(DEA, sep = "\t", header = TRUE)
    deg_table <- deg_table[deg_table$p_val_adj < MAX_PVAL,]
    idx <- which(deg_table$avg_log2FC < MAX_LOG2FC)
    
    if (length(genes) == 0) {
        genes <- deg_table$geneID[idx]
    } else {
        genes <- genes[genes %in% deg_table$geneID[idx]]
    }
}

write.table(genes, file = SIGNATURE, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)


q()
