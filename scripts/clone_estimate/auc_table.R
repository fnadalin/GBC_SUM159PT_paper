
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
    cat("\nUsage: <data.frame> <response> <out>\n")
    cat("\n<data.frame>    data frame containing a predictor in each column\n")
    cat("<response>      rownames subset, one name per line\n")
    cat("<out>           output named vector with AUC values\n\n")
    q()
}

DF <- args[1]
RES <- args[2]
OUT <- args[3]

library("pROC")

df <- read.table(DF, header = TRUE, sep = "\t")
res <- read.table(RES, sep = "\t")[,1]

res <- res[res %in% rownames(df)]
is.sel <- as.numeric(rownames(df) %in% res)

aucval <- c()
for (i in 1:ncol(df)) {
    df_tmp <- data.frame(a = is.sel, b = df[,i])
    rocobj <- roc(data = df_tmp, response = "a", predictor = "b")
    aucval <- c(aucval, rocobj$auc)
}
names(aucval) <- colnames(df)

write.table(aucval, file = OUT, col.names = FALSE, quote = FALSE, sep = "\t")

q()



