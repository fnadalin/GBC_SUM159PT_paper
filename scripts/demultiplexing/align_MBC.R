# conda activate /hpcnfs/data/FN/R403

# .libPaths("/hpcnfs/data/FN/R403/lib/R/library")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
    cat("\nUsage: Rscript call_MBC_myFunc.R <in_dir>\n")
    cat("\n<in_dir>       dir containing MBC.txt, barcodes.tsv.gz, lanes.csv\n\n")
    q()
}
      
IN_DIR <- args[1]

library("deMULTIplex")
# require("KernSmooth")
# library("ggplot2")

setwd(IN_DIR)

MBC_FILE <- "MBC.txt"
CB_FILE <- "barcodes.tsv.gz"
LANES <- "lanes.csv"

OUT_DIR <- "out"
dir.create(OUT_DIR, showWarnings = FALSE)

MBC <- read.table(MBC_FILE, sep = "\t")
bar.ref <- MBC$V2
names(bar.ref) <- MBC$V1

cell.id.vec <- drop(as.matrix(read.table(gzfile(CB_FILE))))
cell.id.vec <- gsub("-1$", "", cell.id.vec)

lanes <- read.table(LANES, header = TRUE, sep = ",")

for (i in 1:nrow(lanes)) {
    R1 <- lanes$R1[i]
    R2 <- lanes$R2[i]
    readTable <- MULTIseq.preProcess(R1 = R1, R2 = R2, cellIDs = cell.id.vec, cell=c(1,16), umi=c(17,28), tag=c(1,8))
    bar.table.l1 <- MULTIseq.align(readTable, cell.id.vec, bar.ref)
    if (i == 1) {
        bar.table <- bar.table.l1
    } else {
        bar.table <- bar.table + bar.table.l1
    }
}

colnames(bar.table)[1:length(bar.ref)] <- names(bar.ref)
outfile <- file.path(OUT_DIR, "bar.table.tsv")
write.table(bar.table, file = outfile, sep = "\t")


q()

