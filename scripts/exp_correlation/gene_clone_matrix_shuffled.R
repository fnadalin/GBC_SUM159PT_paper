
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
    cat("\nUsage: <object> <out_dir> <nsamples> <seed> <samples>\n")
    cat("\n<object>      Seurat object containing the <sampleID>\n")
    cat("<out_dir>     where output correlations are saved\n")
    cat("<seeds>       file containing random seeds, one per line\n")
    cat("<samples>     comma-separated list of samples to consider\n")
    q()
}

OBJECT <- args[1]
OUT_DIR <- args[2]
SEEDS <- args[3]
SAMPLES <- args[4]

library("Seurat")
library("Matrix")

object <- readRDS(OBJECT)
samples <- unlist(strsplit(SAMPLES, split = ","))

seeds <- read.table(SEEDS)[,1]

for (s in samples) {

    # extract cells in selected sample only
    idx <- which(object@meta.data$sample.name == s)
    cells <- colnames(object)[idx]
    n_cells <- length(cells)
    obj_sub <- subset(object, cells = cells)
    
    # subset on the clones that are present in the sample! 
    # otherwise, zero expression is going to be misleading...
    clones <- unique(as.character(obj_sub@meta.data$expr.GBC.list))
    
    # print rows and columns 
    outdir <- file.path(OUT_DIR, s, "shuffled")
    dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

    rows <- file.path(outdir, "features.tsv")
    cols <- file.path(outdir, "clones.tsv")
    
    write.table(rownames(outdir), file = rows, row.names = FALSE, col.names = FALSE, quote = FALSE)
    write.table(clones, file = cols, row.names = FALSE, col.names = FALSE, quote = FALSE)

    system(paste("gzip -f", rows))
    system(paste("gzip -f", cols))
    
    for (seed in seeds) {

        set.seed(seed)
        x <- sample(1:n_cells)
        lab_shuff <- obj_sub@meta.data$expr.GBC.list[x]
        
        # compute the UMI sum for each clone in the sample (raw)
        cloneUMI <- matrix(NA, nrow=nrow(obj_sub), ncol=length(clones))
        for (i in 1:length(clones)) {
            idx_c <- which(lab_shuff == clones[i])
            if (length(idx_c) > 1) {
                cloneUMI[,i] <- rowSums(obj_sub@assays$RNA@counts[,idx_c])
            } else {
                cloneUMI[,i] <- obj_sub@assays$RNA@counts[,idx_c]
            }
        }
        
        # save as sparse matrix 
        M <- Matrix(cloneUMI)
        outdir_shuffled <- file.path(outdir, seed, "raw")
        dir.create(outdir_shuffled, recursive = TRUE, showWarnings = FALSE)
        
        mtx <- file.path(outdir_shuffled, "matrix.mtx")
        
        writeMM(M, file = mtx)
        
        system(paste("gzip -f", mtx))
        system(paste("cd", outdir_shuffled, "; ln -s ../../features.tsv.gz; ln -s ../../clones.tsv.gz"))
        
        # compute the average expression for each clone in the sample (norm)
        cloneExpr <- matrix(NA, nrow=nrow(obj_sub), ncol=length(clones))
        for (i in 1:length(clones)) {
            idx_c <- which(lab_shuff == clones[i])
            if (length(idx_c) > 1) {
                cloneExpr[,i] <- apply(obj_sub@assays$RNA@data[,idx_c], 1, function(xx) log(mean(exp(xx))))
            } else {
                cloneExpr[,i] <- obj_sub@assays$RNA@data[,idx_c]
            }
        }
 
        # save as sparse matrix 
        M <- Matrix(cloneExpr)
        outdir_shuffled <- file.path(outdir, seed, "norm")
        dir.create(outdir_shuffled, recursive = TRUE, showWarnings = FALSE)
        
        mtx <- file.path(outdir_shuffled, "matrix.mtx")
        
        writeMM(M, file = mtx)
        
        system(paste("gzip -f", mtx))
        system(paste("cd", outdir_shuffled, "; ln -s ../../features.tsv.gz; ln -s ../../clones.tsv.gz"))
    }   
}

q()




q()
