args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
    cat("\nUsage: <infile> <samples> <top_clones> <out>\n")
    cat("\n<infile>     data frame with stats (rows are clones and columns are samples)\n")
    cat("<samples>    columns to select\n")
    cat("<top_clones> how many top clones to select (max average stats across samples)\n\n")
    cat("<out>        output file with the list of selected clones, one per line")
    q()
}

CPT <- args[1]
SAMPLES <- args[2]
TOP_CLONES <- args[3]
OUT <- args[4]

samples <- unlist(strsplit(SAMPLES, split = ",")) 

# select the top cpm clones

cpt <- read.table(CPT, sep = "\t", header = TRUE)
cpt <- cpt[,colnames(cpt) %in% samples, drop = FALSE] 

idx <- order(rowSums(cpt), decreasing = TRUE)[1:TOP_CLONES]
CLONE_LIST <- rownames(cpt)[idx]

write.table(CLONE_LIST, file = OUT, quote = FALSE, col.names = FALSE, row.names = FALSE)

q()


