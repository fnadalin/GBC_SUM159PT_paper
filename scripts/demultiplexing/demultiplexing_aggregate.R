
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
    cat("\nUsage: <demux_info> <outdir>\n\n")
    q()
}

DEMUX_INFO <- args[1]
OUTDIR <- args[2]

library("Matrix")

demux_info <- read.table(DEMUX_INFO, sep = "\t")
colnames(demux_info) <- c("sampleID","MBC","matrix.dir")

SAMPLES <- unique(demux_info$sampleID)

for (s in SAMPLES) {
    
    outdir <- file.path(OUTDIR, s)
    dirs <- demux_info$matrix.dir[demux_info$sampleID == s]
    mbc <- demux_info$MBC[demux_info$sampleID == s]
    
    i <- 1
    M_aggr <- NULL
    for (d in dirs) {

        M <- readMM(file.path(d,"matrix.mtx"))
        barcodes <- drop(as.matrix(read.table(file.path(d,"barcodes.tsv"))))
        features <- read.table(file.path(d,"features.tsv"), sep = "\t")
        barcodes <- gsub("$",paste0("-",i),barcodes)
        colnames(M) <- barcodes
        rownames(M) <- features$V2 

        # update the aggregate matrix
        if (is.null(M_aggr)) {
           M_aggr <- M
        } else {
            M_aggr <- cbind(M_aggr, M)
        }
        
        i <- i + 1
    }
    
    dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
    writeMM(M_aggr, file = file.path(outdir, "matrix.mtx"))
    write.table(features, file = file.path(outdir, "features.tsv"), quote = FALSE, col.names = FALSE, row.names = FALSE, sep = "\t")
    write.table(colnames(M_aggr), file = file.path(outdir, "barcodes.tsv"), quote = FALSE, col.names = FALSE, row.names = FALSE)
    
    system(paste("gzip", file.path(outdir, "matrix.mtx")))
    system(paste("gzip", file.path(outdir, "features.tsv")))
    system(paste("gzip", file.path(outdir, "barcodes.tsv")))
}


q()
