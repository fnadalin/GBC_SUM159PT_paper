
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
    cat("\nUsage: <clone_count> <out>\n")
    cat("\n<clone_count>  contains cell count for each clone\n")
    cat("<out>          table containing the library-normalized count (clones per thousand)\n\n")
    q()
}

IN <- args[1]
OUT <- args[2]

count <- read.table(IN, header = TRUE, sep = "\t")
count <- as.matrix(count)

# compute cpt
lib_size <- colSums(count)
cpt <- t(1e+3*t(count)/lib_size)

# write to file
write.table(cpt, file = OUT, sep = "\t", quote = FALSE)


q()
