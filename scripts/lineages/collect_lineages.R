args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
    cat("\nUsage: <file_list> <exp_list> <out>\n")
    cat("\n<file_list>  comma-separated list of files containing the lineage information\n")
    cat("<exp_list>   comma-separated list of experiments referring to <file_list>\n")
    cat("<out>        txt file containing the overall lineage information\n\n")
    q()
}

FILE_LIST <- args[1]
EXP_LIST <- args[2]
OUT <- args[3]

file_list <- unlist(strsplit(FILE_LIST, split = ","))
exp_list <- unlist(strsplit(EXP_LIST, split = ","))

# collect all clone IDs 
names <- c()
for (f in file_list) {
    lin <- read.table(f, sep = "\t")
    names <- c(names,lin[,1])
}
names <- unique(names)
n <- length(names)

# populate the data frame
df <- as.data.frame(matrix(0, nrow = n, ncol = length(file_list)))
rownames(df) <- names
i <- 1
for (f in file_list) {
    lin <- read.table(f, sep = "\t")
    common_names <- names[names %in% lin[,1]]
    idx1 <- match(common_names,names)
    idx2 <- match(common_names,lin[,1])
    df[idx1,i] <- lin[idx2,2]
    i <- i + 1
}

# re-label the lineages according to the best concordance between samples
for (i in 2:ncol(df)) {
    l1 <- sort(unique(df[,1]))
    l2 <- sort(unique(df[,i]))
    l1 <- l1[l1 > 0] # 0 means no lineage assignment
    l2 <- l2[l2 > 0]
    M <- matrix(0, nrow = length(l1), ncol = length(l2))
    for (j in 1:nrow(M)) {
        for (k in 1:ncol(M)) {
            M[j,k] <- sum(df[,1] == l1[j] & df[,i] == l2[k])
        }
    }
    l2_renamed <- apply(M, 2, function(x) l1[which(x == max(x))])
    v <- df[,i]
    for (j in 1:length(l2_renamed)) {
        idx <- which(df[,i] == j)
        v[idx] <- rep(l2_renamed[j],length(idx))
    }
    df[,i] <- v
}

# count the number of lineages each clone is assigned to
count <- apply(df, 1, function(x) {
               v <- unique(drop(x))
               sum(v > 0) })

# count the number of samples a clone is assigned to a lineage for each experiment
df_exp <- as.data.frame(matrix(0, nrow = nrow(df), ncol = length(exp_list)))
i <- 1
for (exp in unique(exp_list)) {
    df_exp[,i] <- rowSums(df[,exp_list == exp] > 0)
}

idx <- which(count == 1 & rowSums(df_exp) == length(unique(exp_list)))
robust_id <- names[idx] # only clones consistently assigned to the same lineage and found at least one sample in all experiments are kept

robust_lin <- apply(df[idx,], 1, function(x) {
               v <- unique(drop(x))
               v[v > 0] })

df <- data.frame(robust_id, robust_lin)

write.table(df, file = OUT, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

q()

