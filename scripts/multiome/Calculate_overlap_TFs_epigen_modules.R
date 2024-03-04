# Script that takes a list of TFs and return the matrix showing the overlap size between each module

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  cat("\nUsage: <max_pval> <out> <dea1> [<tfea2> <tfea3> ... <tfeaN>] \n")
  cat("<max_pval>   adj. p-value cutoff\n")
  cat("<out>        output file containing the matrix overlap\n")
  cat("<tfeaI>       input tsv file containing TF enrichments (TF and FDR fields at least)\n\n")
  q()
}

PVAL.THRESH <- as.numeric(args[1])
OUT.FILE <- as.character(args[2])

TF.lists <- lapply(X = args[3:length(args)], FUN=read.table, header=TRUE)
overlap.matrix.names <- sapply(args[3:length(args)], function(x) sub(".txt","",basename(x)))

filtered.TF.lists <- lapply(TF.lists, function(df) {df[df$FDR < PVAL.THRESH,]})

length.intersect <- function(x,y){length(intersect(x$TF, y$TF))}

overlap.matrix <- outer(filtered.TF.lists, filtered.TF.lists, Vectorize(length.intersect))

colnames(overlap.matrix) <- overlap.matrix.names
rownames(overlap.matrix) <- overlap.matrix.names

print(overlap.matrix)

write.table(overlap.matrix, file = OUT.FILE, sep="\t", quote = FALSE)
q()
