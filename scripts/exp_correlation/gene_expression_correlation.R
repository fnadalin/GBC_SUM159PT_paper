
# compute 3 correlation coefficient between each gene's expression in sample S (averaged across cells) and clone expansion in sample S'
# only clones that are present in both S and S' are considered, otherwise a zero expression is confounded with no clone present

MIN_EXPR <- 20 # compute correlation only for genes detected in >= MIN_EXPR clones

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
    cat("\nUsage: <cl_dir> <mtx_dir> <out_dir>\n")
    cat("\n<cl_dir>      contains counts.tsv and log2fc.tsv\n")
    cat("<mtx_dir>     where the gene-clone matrices (raw and norm) are stored\n")
    cat("<out_dir>     where output correlations are saved\n\n")
    q()
}

CL_DIR <- args[1]
MTX_DIR <- args[2]
OUT_DIR <- args[3]

# library("ggplot2")
library("Matrix")
library("pcaPP") # for cor.fk
library("future.apply")  # for future_lapply

plan(multisession)

COUNTS <- file.path(CL_DIR, "counts.tsv")
LOG2FC <- file.path(CL_DIR, "log2fc.tsv")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

counts <- read.table(COUNTS, sep = "\t", header = TRUE)
log2fc <- read.table(LOG2FC, sep = "\t", header = TRUE)

cloneUMI <- readMM(gzfile(file.path(MTX_DIR, "raw/matrix.mtx.gz")))
clones <- colnames(cloneUMI) <- read.table(gzfile(file.path(MTX_DIR, "raw/clones.tsv.gz")))[,1]
genes <- rownames(cloneUMI) <- read.table(gzfile(file.path(MTX_DIR, "raw/features.tsv.gz")))[,1]

cloneExpr <- readMM(gzfile(file.path(MTX_DIR, "norm/matrix.mtx.gz")))
rownames(cloneExpr) <- read.table(gzfile(file.path(MTX_DIR, "norm/features.tsv.gz")))[,1]
colnames(cloneExpr) <- read.table(gzfile(file.path(MTX_DIR, "norm/clones.tsv.gz")))[,1]

# subset on clones that are present in the gene-clone matrix
clones_in_exp <- rownames(counts)
sel_clones <- clones_in_exp[clones_in_exp %in% clones]
log2fc_sel <- log2fc[match(sel_clones, clones_in_exp),]

# UMI count in selected clones only
clone_umi_sel <- cloneUMI[,match(sel_clones, clones)]

# for each gene, number of clones where it is detected 
num_non_zero_clones <- rowSums(clone_umi_sel > 0)
sel_genes <- genes[num_non_zero_clones >= MIN_EXPR]
    
# gene expression of selected genes and selected clones only
clone_expr_sel <- cloneExpr[match(sel_genes,genes),match(sel_clones, clones)]


################

samples_late <- colnames(log2fc_sel)
for (s in samples_late) {

    # expansion rate of selected clones in current (late) sample
    exp_rate <- log2fc_sel[[s]]

    # compute the correlation between each gene's expression and the expansion rate across clones
    ng <- length(sel_genes)
    lin_corr_clones <- unlist(future_lapply(1:ng, function(i) cor(exp_rate, clone_expr_sel[i,], method = "pearson"), future.packages = "Matrix"))
    rank_corr_clones <- unlist(future_lapply(1:ng, function(i) cor(exp_rate, clone_expr_sel[i,], method = "spearman"), future.packages = "Matrix"))
    rank_k_corr_clones <- unlist(future_lapply(1:ng, function(i) cor.fk(exp_rate, clone_expr_sel[i,]), future.packages = "Matrix"))
    names(rank_corr_clones) <- names(rank_k_corr_clones) <- names(lin_corr_clones) <- sel_genes
    
    filename <- file.path(OUT_DIR, paste0("rank_corr_clones_", s, ".tsv"))
    write.table(rank_corr_clones, filename, col.names = FALSE, sep = "\t", quote = FALSE)
    filename <- file.path(OUT_DIR, paste0("hist_rank_corr_clones_", s, ".pdf"))
    pdf(filename, width = 5, height = 5)
    hist(rank_corr_clones[!is.na(rank_corr_clones)], xlab = "rho", main = paste0("Gene expression VS clonal expansion [", s, "]"))
    dev.off()

    filename <- file.path(OUT_DIR, paste0("rank_k_corr_clones_", s, ".tsv"))
    write.table(rank_k_corr_clones, filename, col.names = FALSE, sep = "\t", quote = FALSE)
    filename <- file.path(OUT_DIR, paste0("hist_rank_k_corr_clones_", s, ".pdf"))
    pdf(filename, width = 5, height = 5)
    hist(rank_k_corr_clones[!is.na(rank_k_corr_clones)], xlab = "tau", main = paste0("Gene expression VS clonal expansion [", s, "]"))
    dev.off()

    filename <- file.path(OUT_DIR, paste0("lin_corr_clones_", s, ".tsv"))
    write.table(lin_corr_clones, filename, col.names = FALSE, sep = "\t", quote = FALSE)
    filename <- file.path(OUT_DIR, paste0("hist_lin_corr_clones_", s, ".pdf"))
    pdf(filename, width = 5, height = 5)
    hist(lin_corr_clones[!is.na(lin_corr_clones)], xlab = "R", main = paste0("Gene expression VS clonal expansion [", s, "]"))
    dev.off()

}


q()
