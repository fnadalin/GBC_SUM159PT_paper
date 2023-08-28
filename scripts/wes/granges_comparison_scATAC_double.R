
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
    cat("\nUsage: <wes> <sc_name> <sc_dir> <cl_info>\n")
    q()
}

WES_EXP_FILE <- args[1]
SC_EXP_NAME <- args[2]
SC_EXP_DIR <- args[3]
CL_INFO <- args[4]

library("GenomicRanges")
library("Matrix")
library("pROC")

FEATURES <- file.path(SC_EXP_DIR, "features.tsv.gz")
BARCODES <- file.path(SC_EXP_DIR, "barcodes.tsv.gz")
MATRIX <- file.path(SC_EXP_DIR, "matrix.mtx.gz")

OUTDIR <- paste0("out_double_treatment_", SC_EXP_NAME, "-WES")
dir.create(OUTDIR, showWarnings = FALSE)

atac <- read.table(FEATURES, sep = "\t")
M <- readMM(MATRIX)
barcodes <- read.table(BARCODES)[,1]

idx <- which(atac[,3] == "Peaks")
n_peaks <- length(idx)

M <- M[idx,]
atac <- atac[idx,]
colnames(M) <- barcodes

gr_atac <- GRanges(seqnames = Rle(atac[,4]), 
                   ranges = IRanges(atac[,5], end = atac[,6], names = 1:n_peaks), 
                   strand = rep("*",n_peaks))

tot_count <- colSums(M)
M_cpm <- 1e+6*t(t(M)/tot_count)

cl_info <- read.table(CL_INFO, sep = ";", header = TRUE)

data <- read.table(WES_EXP_FILE, sep = "\t", header = TRUE)
n_cnv <- nrow(data)
gr <- GRanges(seqnames = data[,1], 
              ranges = IRanges(data[,2], end = data[,3]), 
              strand = rep("*",nrow(data)))

M_cnv_cov <- matrix(0, nrow = n_cnv, ncol = ncol(M))
colnames(M_cnv_cov) <- colnames(M)
coord <- apply(data.frame(start(gr@ranges), end(gr@ranges)), 1, function(x) paste(x, collapse = "-"))
rownames(M_cnv_cov) <- apply(data.frame(gr@seqnames, coord), 1, function(x) paste(x, collapse = ":"))
for (k in 1:n_cnv) {
    ov_loc <- subsetByOverlaps(gr_atac, gr[k])
    idx1 <- as.numeric(names(ov_loc))
    if (sum(width(ov_loc)) == 0) {
        norm_counts <- rep(0,ncol(M_cpm))
    } else {
        if (length(idx1) == 1) {
            norm_counts <- 1e+3*M_cpm[idx1,]/sum(width(ov_loc)) # compute sort of RPKM
        } else {
            norm_counts <- 1e+3*colSums(M_cpm[idx1,])/sum(width(ov_loc)) # compute sort of RPKM
        }
    }
    M_cnv_cov[k,] <- norm_counts
}

idx <- match(cl_info[,1], colnames(M_cnv_cov))
M_cnv_cov <- M_cnv_cov[,idx]
clusters <- cl_info[,2]

write.table(M_cnv_cov, file = file.path(OUTDIR,"matrix.tsv"), sep = "\t", quote = FALSE)
write.table(clusters, file = file.path(OUTDIR,"cl_annotation.txt"), row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(data[,4:ncol(data)], file = file.path(OUTDIR, "logfc_annotation.tsv"), row.names = FALSE, col.names = FALSE, quote = FALSE)


q()

