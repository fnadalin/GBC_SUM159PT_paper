args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
    cat("\nUsage: <file_list> <rows> <out>\n")
    cat("\n<file_list>  comma-separated list of files containing the (same) clone stat in different experiments\n")
    cat("<rows>       file containing the row names to be saved in the matrix (possibly intersected with the rownames of the files above)\n")
    cat("<out>        txt file containing the DRCs (union of the above)\n")
    q()
}

FILE_LIST <- args[1]
ROWS <- args[2]
OUT <- args[3]

file_list <- unlist(strsplit(FILE_LIST, split = ","))

ncols <- 0
columns <- c()
for (file in file_list) {
    df <- read.table(file, sep = "\t", header = TRUE)
    ncols <- ncols + ncol(df)
    columns <- c(columns, colnames(df))
}

union <- read.table(ROWS, sep = "\t")[,1]

M <- matrix(0, nrow = length(union), ncol = ncols)
rownames(M) <- union
colnames(M) <- columns

ncols <- 0
for (file in file_list) {

    df <- read.table(file, sep = "\t", header = TRUE)

    gbc <- union[union %in% rownames(df)]
    idx1 <- which(union %in% gbc)
    idx2 <- match(gbc, rownames(df))
    for (i in 1:ncol(df)) {
        M[idx1,ncols+i] <- df[idx2,i]
    }
    
    ncols <- ncols + ncol(df)
}

write.table(M, file = OUT, sep = "\t", quote = FALSE)

q()


