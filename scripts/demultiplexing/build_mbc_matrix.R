
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
    cat("\nUsage: <CR_matrix_dir> <bartable> <pseudo_CR_matrix_out>\n\n")
    q()
}

INDIR <- args[1]
BARTABLE <- args[2]
OUTDIR <- args[3]

library("Matrix")

M <- readMM(gzfile(file.path(INDIR,"matrix.mtx.gz")))
barcodes <- drop(as.matrix(read.table(gzfile(file.path(INDIR,"barcodes.tsv.gz")))))
features <- read.table(gzfile(file.path(INDIR,"features.tsv.gz")), sep = "\t")
barcodes <- gsub("-1","",barcodes)
colnames(M) <- barcodes
rownames(M) <- features$V2 

bar.table <- read.table(BARTABLE, sep = "\t", header = TRUE)
NMBC <- ncol(bar.table)-2
M_mbc_only <- t(bar.table[,1:NMBC])

# concatenate
M_with_mbc <- rbind(M, M_mbc_only)

# update features
mbc_names <- rownames(M_mbc_only)
df <- data.frame(V1 = mbc_names, V2 = mbc_names, V3 = rep("CRISPR Guide Capture",NMBC))
features <- rbind(features, df)

dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
writeMM(M_with_mbc, file = file.path(OUTDIR, "matrix.mtx"))
write.table(features, file = file.path(OUTDIR, "features.tsv"), quote = FALSE, col.names = FALSE, row.names = FALSE, sep = "\t")
write.table(barcodes, file = file.path(OUTDIR, "barcodes.tsv"), quote = FALSE, col.names = FALSE, row.names = FALSE)

q()

