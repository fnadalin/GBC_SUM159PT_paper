
# pseudo-count
PSC <- 0.5

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
    cat("\nUsage: <in.tsv> <out.tsv> <control> <case>\n")
    cat("\n<in.tsv>      table containing the library-normalized count (cpt)\n")
    cat("<out.tsv>     table containing the log2FC of each case vs the average of the controls\n")
    cat("<control>     comma-separated list of controls\n")
    cat("<case>        comma-separated list of cases\n\n")
    q()
}

IN <- args[1]
OUT <- args[2]
CTRL <- args[3]
CASE <- args[4]

cpt <- read.table(IN, header = TRUE, sep = "\t")

ctrl <- unlist(strsplit(CTRL, split=","))
case <- unlist(strsplit(CASE, split=","))

# compute cpt average across controls
idx <- match(ctrl, colnames(cpt))
if (length(idx) == 1) {
    cpt_P0_avg <- cpt[,idx]
} else {
    cpt_P0_avg <- apply(cpt[,idx], 1, mean)
}

# initialize the log2FC matrix
idx <- match(case, colnames(cpt))
log2FC <- matrix(NA, nrow = nrow(cpt), ncol = length(idx))
rownames(log2FC) <- rownames(cpt)
colnames(log2FC) <- colnames(cpt)[idx]

# populate the log2FC matrix
for (i in 1:length(idx)) {
    log2FC[,i] <- log2((cpt[,idx[i]]+PSC)/(cpt_P0_avg+PSC))
}

# write to file
write.table(log2FC, file = OUT, sep = "\t", quote = FALSE)

q()


