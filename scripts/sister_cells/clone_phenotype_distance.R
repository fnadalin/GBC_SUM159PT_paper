
non.rand.sd.frac <- 0.5
DISTANCE <- "euclidean"
SEED <- 1
# NUM_ITER <- 10
NUM_ITER <- 1

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) {
        cat("\nUsage: Rscript clone_phenotype_distance.R <obj> <slot> <dr> <nPC> <outdir>\n")
        q()
}

OBJ <- args[1]
SLOT <- args[2]
DR <- args[3]
NPC <- args[4]
OUTDIR <- args[5]

dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

library("Seurat")

# cells in objects have been filtered for doublets and low quality cells
# they also contain relevant metadata
obj <- readRDS(OBJ)
ncells <- ncol(obj)

# extract phenotype
phenotype <- obj@meta.data[[SLOT]]

# compute distances in the PCA space
data <- as.matrix(obj@reductions[[DR]]@cell.embeddings[,1:NPC])

intra_pheno_1_dist <- NULL
intra_pheno_0_dist <- NULL
inter_pheno_dist <- NULL

n <- min(sum(phenotype),sum(!phenotype))
idx_1 <- which(phenotype == 1)
idx_0 <- which(phenotype == 0)
n_1 <- length(idx_1)
n_0 <- length(idx_0)

set.seed(SEED)
for (iter in 1:NUM_ITER) {
   
    idx_a <- sample(n_1)[1:n]
    idx_b <- sample(n_1)[1:n]
    distances <- as.matrix(dist(data[c(idx_1[idx_a],idx_1[idx_b]),], method = DISTANCE))
    M <- distances[1:n,n+1:n]
    intra_pheno_1_dist <- c(intra_pheno_1_dist, drop(M))

    idx_a <- sample(n_0)[1:n]
    idx_b <- sample(n_0)[1:n]
    distances <- as.matrix(dist(data[c(idx_0[idx_a],idx_0[idx_b]),], method = DISTANCE))
    M <- distances[1:n,n+1:n]
    intra_pheno_0_dist <- c(intra_pheno_0_dist, drop(M))

    idx_a <- sample(n_1)[1:n]
    idx_b <- sample(n_0)[1:n]
    distances <- as.matrix(dist(data[c(idx_1[idx_a],idx_0[idx_b]),], method = DISTANCE))
    M <- distances[1:n,n+1:n]
    inter_pheno_dist <- c(inter_pheno_dist, drop(M))
    
}

max_dist <- max(intra_pheno_1_dist, intra_pheno_0_dist, inter_pheno_dist)

# plot density
plotfile <- file.path(OUTDIR, paste0("clone_pheno_dist_density.pdf"))
pdf(plotfile, width = 5, height = 5)
d_pheno_1 <- density(intra_pheno_1_dist)
d_pheno_0 <- density(intra_pheno_0_dist)
d_inter <- density(inter_pheno_dist)
max_density <- max(c(d_pheno_1$y,d_pheno_0$y,d_inter$y))
plot(d_inter, lty = 2, lwd = 2, xlim = c(0,max_dist), ylim = c(0,max_density), main = "Cell-cell distance", xlab = "euclidean distance", ylab = "density")
lines(d_pheno_0, col = "blue", lwd = 2)
lines(d_pheno_1, col = "red", lwd = 2) 
legend("topright", legend = c("positive pheno","negative pheno","different pheno"), lty = c(1,1,2), lwd = 2, col = c("red","blue","black"), bty = "n")
dev.off()

    


q()

