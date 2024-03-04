SUBPOP <- c("S1","S2","S3")
LABELS_HIGHLIGHT <- 1:24 # only highlight the programs common to more than one tumour

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
    cat("\nUsage: <dir1> <dir2> <outdir>\n")
    q()
}

DIR1 <- args[1]
DIR2 <- args[2]
OUT_DIR <- args[3]

library("ComplexHeatmap")

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

for (type in c("bycell","bygene")) {

    basename <- paste0("auc_",type,".tsv")
    df1 <- read.table(file.path(DIR1, basename), sep = "\t", header = TRUE, check.names = FALSE)
    df2 <- read.table(file.path(DIR2, basename), sep = "\t", header = TRUE, check.names = FALSE)
    
    df1 <- df1[match(SUBPOP, rownames(df1)),]
    df2 <- df2[match(SUBPOP, rownames(df2)),]
    rownames(df1) <- paste0(rownames(df1), rep("_1",length(SUBPOP)))
    rownames(df2) <- paste0(rownames(df2), rep("_2",length(SUBPOP)))
    df <- rbind(df1, df2)
    idx <- unlist(lapply(1:3, function(x) c(0,3)+x))
    df <- df[idx,]
    
    fontface <- rep("italic",ncol(df))
    fontface[LABELS_HIGHLIGHT] <- rep("bold",length(LABELS_HIGHLIGHT))

    basename <- paste0("auc_",type,".pdf")
    pdf(file.path(OUT_DIR, basename), width = 4.7, height = 7.5)
    draw(Heatmap(t(as.matrix(df)), name = paste0("auc_", type), cluster_columns = FALSE, column_split = c(1,1,2,2,3,3), 
    		row_names_gp = gpar(fontface = fontface)))
    dev.off()
    
    basename <- paste0("auc_",type,"_common.pdf")
    pdf(file.path(OUT_DIR, basename), width = 4.7, height = 5.5)
    draw(Heatmap(t(as.matrix(df[,LABELS_HIGHLIGHT])), name = paste0("auc_", type), cluster_columns = FALSE, column_split = c(1,1,2,2,3,3)))
    dev.off()
    
    df_mean <- sapply(1:length(SUBPOP), function(x) apply(df[c(2*x-1,2*x),], 2, mean))
    colnames(df_mean) <- SUBPOP
    rownames(df_mean) <- colnames(df)
    
    basename <- paste0("auc_",type,"_mean.pdf")
    pdf(file.path(OUT_DIR, basename), width = 4.7, height = 7.5)
    draw(Heatmap(as.matrix(df_mean), name = paste0("auc_", type), cluster_columns = FALSE, row_names_gp = gpar(fontface = fontface)))
    dev.off()
    
    basename <- paste0("auc_",type,"_common_mean.pdf")
    pdf(file.path(OUT_DIR, basename), width = 4.7, height = 5.5)
    draw(Heatmap(as.matrix(df_mean[LABELS_HIGHLIGHT,]), name = paste0("auc_", type), cluster_columns = FALSE))
    dev.off()
    
}

q()


