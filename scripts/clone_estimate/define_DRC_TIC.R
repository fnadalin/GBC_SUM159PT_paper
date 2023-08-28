
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
    cat("\nUsage: <indir> <outdir>\n")
    cat("\n<indir>    where cpm_count_stat.tsv and count_sel_stat.tsv are stored\n")
    cat("<outdir>  where gain values are stored\n\n")
    q()
}

INDIR <- args[1]
OUTDIR <- args[2]

COUNT_SEL <- file.path(INDIR, "count_sel.tsv")
GAIN <- file.path(OUTDIR, "gain.tsv")

# select the min count as the average minimum count across samples such that gain > 1
gain <- read.table(GAIN, sep = "\t", header = TRUE) 
m <- nrow(gain)
n <- ncol(gain)
idx <- apply(gain, 2, function(x) which(x > 1)[1])
sel_idx <- floor(mean(idx))

cat(paste0("Min count for DRC/TIC is ", sel_idx, "\n"))

count_sel <- read.table(COUNT_SEL, sep = "\t")
drc_tic <- count_sel[count_sel[,2] >= sel_idx,1]

LIST <- file.path(OUTDIR, "list.txt")

write.table(drc_tic, file = LIST, quote = FALSE, col.names = FALSE, row.names = FALSE)

q()

