
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
    cat("\nUsage: <CB_metadata> <barcode_list> <out>\n")
    q()
}

CB_META <- args[1]
BARCODES <- args[2]
OUT <- args[3]

cb_meta <- read.table(CB_META, sep = "\t", header = TRUE)
if (length(grep("gz$",BARCODES)) > 0) {
    barcodes <- read.table(gzfile(BARCODES), sep = "\t")[,1]
} else {
    barcodes <- read.table(BARCODES, sep = "\t")[,1]
}

idx <- which(rownames(cb_meta) %in% barcodes)
cb_meta <- cb_meta[idx,]

write.table(cb_meta, file = OUT, sep = "\t", quote = FALSE)

q()

