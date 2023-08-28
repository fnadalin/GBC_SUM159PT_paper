
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
    cat("\nUsage: <sc_dataset> <pops> <out_prefix>\n")
    q()
}

OBJECT <- args[1]
POPS <- unlist(strsplit(args[2], split = ","))
OUT_PREFIX <- args[3]

library("Seurat")

object <- readRDS(OBJECT)

n <- length(POPS)
M_count <- matrix(0, nrow = n, ncol = n)
M_odds <- matrix(Inf, nrow = n, ncol = n)
M_pval <- matrix(0, nrow = n, ncol = n)
for (i in 1:(n-1)) {
    meta <- paste0("scissor.", POPS[i])
    x <- as.numeric(object@meta.data[[meta]] == 1)
    M_count[i,i] <- sum(x)
    for (j in (i+1):n) {
        meta <- paste0("scissor.", POPS[j])
        y <- as.numeric(object@meta.data[[meta]] == 1)
        M_count[j,j] <- sum(y)
        a <- sum(x & y)
        b <- sum(x & !y)
        c <- sum(!x & y)
        d <- sum(!x & !y)
        M <- matrix(c(a,b,c,d), nrow = 2, ncol = 2)
        f <- fisher.test(M, alternative = "greater")
        M_count[i,j] <- M_count[j,i] <- a
        M_odds[i,j] <- M_odds[j,i] <- f$estimate
        M_pval[i,j] <- M_pval[j,i] <- f$p.value
    }
}
rownames(M_count) <- colnames(M_count) <- POPS
rownames(M_odds) <- colnames(M_odds) <- POPS
rownames(M_pval) <- colnames(M_pval) <- POPS

write.table(M_count, file = paste0(OUT_PREFIX, "_counts.tsv"), sep = "\t", quote = FALSE)
write.table(M_odds, file = paste0(OUT_PREFIX, "_odds_ratio.tsv"), sep = "\t", quote = FALSE)
write.table(M_pval, file = paste0(OUT_PREFIX, "_pvalue.tsv"), sep = "\t", quote = FALSE)


q()

