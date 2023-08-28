
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
    cat("\nUsage: <in_meta.tsv> <cells> <out_meta.tsv>\n")
    cat("\nsubset metadata on <cells>\n\n")
    q()
}

IN_META <- args[1]
CELLS <- args[2]
OUT_META <- args[3]

in_meta <- read.table(IN_META, sep = "\t", header = TRUE)
cells <- read.table(CELLS)[,1]

out_meta <- in_meta[rownames(in_meta) %in% cells,]

write.table(out_meta, file = OUT_META, sep = "\t", quote = FALSE)


q()

