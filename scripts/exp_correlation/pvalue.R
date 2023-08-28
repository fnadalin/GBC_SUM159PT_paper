
# shuffle clone labels to create null sets
# compute rank correlation coefficient between each gene's expression in sample S (averaged across cells) and clone expansion in sample S'
# only clones that are present in both S and S' are considered, otherwise a zero expression is confounded with no clone present

CORR_ABS <- 0.1
PVAL <- 0.05

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
    cat("\nUsage: <cl_dir> <corr_dir> <shuff_corr_dir>\n")
    cat("\n<cl_dir>          directory containing log2FC of clones\n")
    cat("<corr_dir>        where original correlations are saved\n")
    cat("<shuff_corr_dir>  where shuffled correlations are saved\n\n")
    q()
}

CL_DIR <- args[1]
CORR_DIR <- args[2]
SHUFF_CORR_DIR <- args[3]

LOG2FC <- file.path(CL_DIR, "log2fc.tsv")
log2fc <- read.table(LOG2FC, sep = "\t", header = TRUE)
samples_late <- colnames(log2fc)

for (s in samples_late) {

    corr <- read.table(file.path(CORR_DIR,paste0("rank_corr_clones_",s,".tsv")), sep = "\t")
    corr_perm <- read.table(file.path(SHUFF_CORR_DIR,paste0("rank_corr_clones_",s,"_permutest.tsv")), sep = "\t")
    
    nsamples <- ncol(corr_perm)

    idx <- which(!is.na(corr[,2]))
    corr_filt <- corr[idx,]
    corr_perm <- corr_perm[idx,]
    n <- nrow(corr_filt)
    m <- apply(corr_perm, 1, function(x) sum(!is.na(x))) # N.B. this mostly affects negative correlations, since genes are *low* in many cells

    # compute the fraction of permuted correlations above the top / below the bottom
    pval <- rep(1,n)
    idx <- which(corr_filt[,2] > 0) 
    corr_perm_bool <- (corr_perm[idx,] >= corr_filt[idx,2])
    count <- rowSums(corr_perm_bool, na.rm = TRUE)
    pval[idx] <- count / m[idx]
    idx <- which(corr_filt[,2] < 0)
    corr_perm_bool <- (corr_perm[idx,] <= corr_filt[idx,2])
    count <- rowSums(corr_perm_bool, na.rm = TRUE)
    pval[idx] <- count / m[idx]
    idx <- which(pval==0)
    if (length(idx) > 0) {
        pval[idx] <- rep(1/ncol(corr_perm),length(idx)) # for log-transform
    }

    filename <- file.path(SHUFF_CORR_DIR, paste0("dotplot_rank_corr_perm_",s,".pdf"))
    title <- "Gene correlation"
    subtitle <- paste0(nsamples, "-sample permutation test")
    
    m <- max(abs(corr_filt[,2]))
    xlim <- c(-m,m)
    
    df <- data.frame(corr_filt[,1], corr_filt[,2], pval)
    colnames(df) <- c("id", "cor", "p_value")
    df <- df[!is.na(df$p_value),]
    
    pdf(filename)
    plot(df$cor, -log10(df$p_value), main = title, xlab = "rho", ylab = "-log10(p-value)", col = "gray90", sub = subtitle, xlim = xlim)
    idx <- which(abs(df$cor) > 0.1 & -log10(df$p_value) > 2)
    points(df$cor[idx], -log10(df$p_value)[idx], col = "black")
    dev.off()

    filename <- file.path(SHUFF_CORR_DIR, paste0("rank_corr_pvalue_",s,".tsv"))
    write.table(df, file = filename, row.names = FALSE, quote = FALSE, sep = "\t")

}


q()
