
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
    cat("\nUsage: <indir> <out.txt>\n")
    cat("\n<indir>    where cpm_count_stat.tsv and count_sel_stat.tsv are stored\n")
    cat("<outdir>   where to save the gain values\n\n")
    q()
}

INDIR <- args[1]
OUTDIR <- args[2]

COUNT_SEL <- file.path(INDIR, "count_sel.tsv")
CPM_COUNT_STAT <- file.path(INDIR, "cpm_count_stat.tsv")

count_sel <- read.table(COUNT_SEL, sep = "\t")
cpm_cum <- read.table(CPM_COUNT_STAT, sep = "\t", header = TRUE)

# compute the fraction of cpm gained by including each GBC group (GBCs are grouped by count) for each sample
m <- nrow(cpm_cum)
n <- ncol(cpm_cum)
cpm <- cpm_cum
for (j in 1:n) {
    cpm[m,j] <- cpm_cum[m,j]
    for (i in 1:(m-1)) {
        cpm[i,j] <- cpm_cum[i,j] - cpm_cum[i+1,j]
    }
}

# compute the contribution to the fraction of GBC of each group for each sample
samples <- colnames(cpm_cum)
gbc <- matrix(NA, nrow = m, ncol = n)
for (j in 1:n) { 
    GBC_LIST <- file.path(INDIR, paste0(samples[j], "_sel.txt"))
    gbc_list <- read.table(GBC_LIST)[,1]
    for (i in 1:m) {
        gbc[i,j] <- sum(count_sel[count_sel[,2] == i,1] %in% gbc_list)
    }
}
gbc <- t(t(gbc)/colSums(gbc))

# gain = (average) gain in cpm coverage
# gain = 1 for each GBC <=> uniform coverage
gain <- cpm/gbc
colnames(gain) <- colnames(cpm)
rownames(gain) <- rownames(cpm)

GAIN <- file.path(OUTDIR, "gain.tsv")
write.table(gain, file = GAIN, sep = "\t", quote = FALSE)


q()

