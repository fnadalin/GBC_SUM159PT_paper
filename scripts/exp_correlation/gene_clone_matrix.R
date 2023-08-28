
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
    cat("\nUsage: <object> <out_dir>\n")
    cat("\n<object>      Seurat object containing the <sampleID>\n")
    cat("<out_dir>     where output correlations are saved\n\n")
    q()
}

OBJECT <- args[1]
OUT_DIR <- args[2]

library("Seurat")
library("Matrix")

object <- readRDS(OBJECT)
samples <- unique(as.character(object@meta.data$sample.name))

for (s in samples) {

    # extract cells in selected sample only
    idx <- which(object@meta.data$sample.name == s)
    cells <- colnames(object)[idx]
    obj_sub <- subset(object, cells = cells)
    
    # subset on the clones that are present in the sample! 
    # otherwise, zero expression is going to be misleading...
    clones <- unique(as.character(obj_sub@meta.data$expr.GBC.list))

    # compute the UMI sum for each clone in the sample (raw)
    cloneUMI <- matrix(NA, nrow=nrow(obj_sub), ncol=length(clones))
    for (i in 1:length(clones)) {
        idx_c <- which(obj_sub@meta.data$expr.GBC.list == clones[i])
        if (length(idx_c) > 1) {
            cloneUMI[,i] <- rowSums(obj_sub@assays$RNA@counts[,idx_c])
        } else {
            cloneUMI[,i] <- obj_sub@assays$RNA@counts[,idx_c]
        }
    }

    # compute the average expression for each clone in the sample (norm)
    cloneExpr <- matrix(NA, nrow=nrow(obj_sub), ncol=length(clones))
    for (i in 1:length(clones)) {
        idx_c <- which(obj_sub@meta.data$expr.GBC.list == clones[i])
        if (length(idx_c) > 1) {
            cloneExpr[,i] <- apply(obj_sub@assays$RNA@data[,idx_c], 1, function(xx) log(mean(exp(xx))))
        } else {
            cloneExpr[,i] <- obj_sub@assays$RNA@data[,idx_c]
        }
    }
    
    # save as sparse matrix 
    M <- Matrix(cloneUMI)
    
    outdir <- file.path(OUT_DIR, s, "raw")
    dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
    
    mtx <- file.path(outdir, "matrix.mtx")
    rows <- file.path(outdir, "features.tsv")
    cols <- file.path(outdir, "clones.tsv")
    
    writeMM(M, file = mtx)
    write.table(rownames(object), file = rows, row.names = FALSE, col.names = FALSE, quote = FALSE)
    write.table(clones, file = cols, row.names = FALSE, col.names = FALSE, quote = FALSE)
    
    system(paste("gzip -f", mtx))
    system(paste("gzip -f", rows))
    system(paste("gzip -f", cols))
    
    # save as sparse matrix 
    M <- Matrix(cloneExpr)
    
    outdir <- file.path(OUT_DIR, s, "norm")
    dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
    
    mtx <- file.path(outdir, "matrix.mtx")
    rows <- file.path(outdir, "features.tsv")
    cols <- file.path(outdir, "clones.tsv")
    
    writeMM(M, file = mtx)
    write.table(rownames(object), file = rows, row.names = FALSE, col.names = FALSE, quote = FALSE)
    write.table(clones, file = cols, row.names = FALSE, col.names = FALSE, quote = FALSE)
    
    system(paste("gzip -f", mtx))
    system(paste("gzip -f", rows))
    system(paste("gzip -f", cols))
}

q()




q()
