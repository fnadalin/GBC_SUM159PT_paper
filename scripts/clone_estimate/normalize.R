
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
    cat("\nUsage: <in.tsv> <out.tsv>\n")
    cat("\n<in.tsv>   table containing the count for each GBC (as well as of their neighbors)\n")
    cat("<out.tsv>  table containing the library-normalized count (cpm)\n\n")
    q()
}

IN <- args[1]
OUT <- args[2]

count <- read.table(IN, header = TRUE, sep = "\t")
count <- as.matrix(count)

# compute cpm adding
lib_size <- colSums(count)
cpm <- t(1e+6*t(count)/lib_size)

# write to file
write.table(cpm, file = OUT, sep = "\t", quote = FALSE)

q()

