args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
    cat("\nUsage: <file_list> <out>\n")
    cat("\n<file_list>  comma-separated list of files containing clones in different experiments\n")
    cat("<out_dir>    dir containing the union of the above: cell count and CPT\n\n")
    q()
}

FILE_LIST <- args[1]
OUT_DIR <- args[2]

file_list <- unlist(strsplit(FILE_LIST, split = ","))
clone_id <- c()
n_samples <- 0
for (f in file_list) {
    df <- read.table(f, sep = "\t", header = TRUE)
    clone_id <- c(clone_id, rownames(df))
    n_samples <- n_samples + ncol(df)
}
clone_id <- unique(clone_id)

df_all <- data.frame(matrix(0, nrow = length(clone_id), ncol = n_samples))
n <- 0
for (f in file_list) {
    df <- read.table(f, sep = "\t", header = TRUE)
    idx <- match(rownames(df), clone_id)
    df_all[idx,n+1:ncol(df)] <- df
    colnames(df_all)[n+1:ncol(df)] <- colnames(df)
    n <- n + ncol(df)
}
rownames(df_all) <- clone_id

M <- as.matrix(df_all)
cpm <- t(1000000*t(M)/colSums(M))
cpt <- t(1000*t(M)/colSums(M))

OUT <- file.path(OUT_DIR, "counts.tsv")
write.table(df_all, file = OUT, sep = "\t", quote = FALSE)
OUT <- file.path(OUT_DIR, "cpm.tsv")
write.table(cpm, file = OUT, sep = "\t", quote = FALSE)
OUT <- file.path(OUT_DIR, "cpt.tsv")
write.table(cpt, file = OUT, sep = "\t", quote = FALSE)

q()


