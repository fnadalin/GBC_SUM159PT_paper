
# shuffle clone labels to create null sets
# compute rank correlation coefficient between each gene's expression in sample S (averaged across cells) and clone expansion in sample S'
# only clones that are present in both S and S' are considered, otherwise a zero expression is confounded with no clone present

MIN_EXPR <- 20 # compute correlation only for genes detected in >= MIN_EXPR clones

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
    cat("\nUsage: <object> <sample> <cl_dir> <out_dir> <seeds>\n")
    cat("\n<object>      input Seurat object containing <sample> in meta.data$sample.name\n")
    cat("<sample>      sample ID to extract\n")
    cat("<cl_dir>      contains counts.tsv and log2fc.tsv\n")
    cat("<out_dir>     where output correlations are saved\n")
    cat("<seeds>       file containing seeds, one per line\n\n")
    q()
}

OBJECT <- args[1]
SAMPLE <- args[2]
CL_DIR <- args[3]
OUT_DIR <- args[4]
SEEDS <- args[5]

library("Seurat")
library("Matrix")
# library("ggplot2")
library("Matrix")
library("pcaPP") # for cor.fk
library("future.apply")  # for future_lapply

plan(multisession)

seeds <- read.table(SEEDS)[,1]

COUNTS <- file.path(CL_DIR, "counts.tsv")
LOG2FC <- file.path(CL_DIR, "log2fc.tsv")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

counts <- read.table(COUNTS, sep = "\t", header = TRUE)
log2fc <- read.table(LOG2FC, sep = "\t", header = TRUE)

object <- readRDS(OBJECT)

# extract cells in selected sample only
idx <- which(object@meta.data$sample.name == SAMPLE)
genes <- rownames(object)
cells <- colnames(object)[idx]
n_cells <- length(cells)
object <- subset(object, cells = cells)

# subset on the clones that are present in the sample! 
# otherwise, zero expression is going to be misleading...
clones <- unique(as.character(object@meta.data$expr.GBC.list))

# subset on clones that are present in the gene-clone matrix
clones_in_exp <- rownames(counts)
sel_clones <- clones_in_exp[clones_in_exp %in% clones]
log2fc_sel <- log2fc[match(sel_clones, clones_in_exp),]


################

samples_late <- colnames(log2fc_sel)
for (s in samples_late) {

    # expansion rate of selected clones in current (late) sample
    exp_rate <- log2fc_sel[[s]]

    # dataframe where rows are genes, columns are seeds and entries are correlation values 
    df <- data.frame()

    # compute correlation values
    for (seed in seeds) {

        cat(paste("seed=", seed, "\n"))

        set.seed(seed)
        x <- sample(1:n_cells)
        lab_shuff <- object@meta.data$expr.GBC.list[x]
        
        # compute the UMI sum for each clone in the sample (raw)
        cloneUMI <- matrix(NA, nrow=nrow(object), ncol=length(clones))
        for (i in 1:length(clones)) {
            idx_c <- which(lab_shuff == clones[i])
            if (length(idx_c) > 1) {
                cloneUMI[,i] <- rowSums(object@assays$RNA@counts[,idx_c])
            } else {
                cloneUMI[,i] <- object@assays$RNA@counts[,idx_c]
            }
        }
        colnames(cloneUMI) <- clones
        rownames(cloneUMI) <- genes
  
        # compute the average expression for each clone in the sample (norm)
        cloneExpr <- matrix(NA, nrow=nrow(object), ncol=length(clones))
        for (i in 1:length(clones)) {
            idx_c <- which(lab_shuff == clones[i])
            if (length(idx_c) > 1) {
                cloneExpr[,i] <- apply(object@assays$RNA@data[,idx_c], 1, function(xx) log(mean(exp(xx))))
            } else {
                cloneExpr[,i] <- object@assays$RNA@data[,idx_c]
            }
        }
        colnames(cloneExpr) <- clones
        rownames(cloneExpr) <- genes

        # UMI count in selected clones only
        clone_umi_sel <- cloneUMI[,match(sel_clones, clones)]

        # for each gene, number of clones where it is detected 
        num_non_zero_clones <- rowSums(clone_umi_sel > 0)
        sel_genes <- genes[num_non_zero_clones >= MIN_EXPR]
            
        # gene expression of selected genes and selected clones only
        idx <- match(sel_genes,genes)
        clone_expr_sel <- cloneExpr[match(sel_genes,genes),match(sel_clones, clones)]

        # compute the correlation between each gene's expression and the expansion rate across clones
        ng <- length(sel_genes)
        corr_clones <- unlist(future_lapply(1:ng, function(i) cor(exp_rate, clone_expr_sel[i,], method = "spearman"), future.packages = "Matrix"))
        
        v <- rep(NA, length(genes))
        if (ncol(df) == 0) {
            df <- data.frame(v)
        } else {
            df <- cbind(df, v)
        }
        df[idx,ncol(df)] <- corr_clones
    }
    
    rownames(df) <- genes
    colnames(df) <- seeds
    
    # print the dataframe
    filename <- file.path(OUT_DIR, paste0("rank_corr_clones_", s, "_permutest.tsv"))
    write.table(df, file = filename, quote = FALSE, sep = "\t")

}



q()
