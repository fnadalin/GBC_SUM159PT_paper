
non.rand.sd.frac <- 0.5
DISTANCE <- "euclidean"
SEED <- 1

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 6) {
        cat("\nUsage: Rscript clone_distance_intra_sample.R <obj> <slot> <id> <dr> <nPC> <outdir>\n")
        q()
}

OBJ <- args[1]
SLOT <- args[2]
ID <- args[3]
DR <- args[4]
NPC <- args[5]
OUTDIR <- args[6]

dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

library("Seurat")

# cells in objects have been filtered for doublets and low quality cells
# they also contain relevant metadata
obj <- readRDS(OBJ)

# subset object
cells <- colnames(obj)[obj@meta.data[[SLOT]] == ID] 
obj <- subset(obj, cells = cells)
ncells <- ncol(obj)

# extract clones
expr.GBC.list <- obj@meta.data$expr.GBC.list
clones <- unique(expr.GBC.list)

# compute distances in the PCA space
data <- as.matrix(obj@reductions[[DR]]@cell.embeddings[,1:NPC])

# extract inter-sample distances
intra_clone_dist <- NULL
inter_clone_dist <- NULL
set.seed(SEED)
for (cl in clones) {
    idx_rand <- sample(ncells)
    idx1 <- which(expr.GBC.list == cl)
    if (length(idx1) > 1) {
        idx2 <- which(expr.GBC.list[idx_rand] != cl)[1:length(idx1)]
        dm <- as.matrix(dist(data[idx1,], method = DISTANCE))
        intra_clone_dist <- c(intra_clone_dist, drop(dm))
        distances <- as.matrix(dist(data[c(idx1,idx2),], method = DISTANCE))
        dm <- distances[1:length(idx1),length(idx1)+1:length(idx2)]
        inter_clone_dist <- c(inter_clone_dist, drop(dm))
    }
}
max_dist <- max(c(intra_clone_dist,inter_clone_dist))

# plot density
plotfile <- file.path(OUTDIR, paste0("clone_dist_density.pdf"))
pdf(plotfile, width = 5, height = 5)
d_intra <- density(intra_clone_dist)
d_inter <- density(inter_clone_dist)
max_density <- max(c(d_intra$y,d_inter$y))
plot(d_inter, lty = 2, lwd = 2, xlim = c(0,max_dist), ylim = c(0,max_density), main = "Cell-cell distance", xlab = "euclidean distance", ylab = "density")
lines(d_intra, col = "blue", lwd = 2) 
legend("topright", legend = c("sisters","random pairs"), lty = c(1,2), lwd = 2, col = c("blue","black"), bty = "n")
dev.off()

    


q()

