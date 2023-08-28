args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
    cat("\nUsage: <file_list> <out>\n")
    cat("\n<file_list>  comma-separated list of files containing the lineage information\n")
    cat("<out>        txt file containing the overall lineage information\n")
    q()
}

FILE_LIST <- args[1]
OUT <- args[2]

file_list <- unlist(strsplit(FILE_LIST, split = ","))

# check that rownames are all the same
names <- c()
for (f in file_list) {
    lin <- read.table(f, sep = "\t")
    if (length(names) > 0) {
        if (sum(lin[,1] == names) != nrow(lin)) {
            cat("different clone names across samples\n")
            q()
        } 
    } else {
        names <- lin[,1]
    }
}
n <- length(names)

# populate the data frame
df <- as.data.frame(matrix(NA, nrow = n, ncol = 0))
for (f in file_list) {
    lin <- read.table(f, sep = "\t")
    df <- cbind(df, lin[,2])
}
rownames(df) <- names

# count the number of lineages each clone is assigned to
count <- apply(df, 1, function(x) {
               v <- unique(drop(x))
               sum(!is.na(v)) })

idx <- which(count == 1)
robust_id <- names[count == 1] # only clones consistently assigned to the same lineage are kept

robust_lin <- apply(df[idx,], 1, function(x) {
               v <- unique(drop(x))
               v[!is.na(v)] })

df <- data.frame(robust_id, robust_lin)

write.table(df, file = OUT, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

q()

