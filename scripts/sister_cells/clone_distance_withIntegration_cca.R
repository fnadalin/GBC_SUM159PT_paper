
non.rand.sd.frac <- 0.5
DISTANCE <- "euclidean"
SEED <- 1

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 7) {
        cat("\nUsage: Rscript clone_distance_withIntegration_cca.R <obj1> <obj2> <cl1> <cl2> <hvg1> <hvg2> <outdir>\n")
        q()
}

OBJ1 <- args[1]
OBJ2 <- args[2]
CL1 <- args[3]
CL2 <- args[4]
HVG1 <- args[5]
HVG2 <- args[6]
OUTDIR <- args[7]

dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

library("Seurat")

# cells in objects have been filtered for doublets and low quality cells
# they also contain relevant metadata
obj1 <- readRDS(OBJ1)
obj2 <- readRDS(OBJ2)
ncells1 <- ncol(obj1)
ncells2 <- ncol(obj2)
cl1 <- obj1@meta.data[[CL1]]
cl2 <- obj2@meta.data[[CL2]]

# subset on common genes
genes1 <- rownames(obj1)
genes2 <- rownames(obj2)
# extract common hvgs
hvg1 <- read.table(HVG1)[,1]
hvg2 <- read.table(HVG2)[,1]
hvg <- hvg1[hvg1 %in% hvg2]

obj1 <- subset(obj1, features = hvg)
obj2 <- subset(obj2, features = hvg)

# perform integration
object <- RunCCA(obj1, obj2, set.seed = SEED)

nCCs <- 20
object <- RunUMAP(object, dims = 1:nCCs, reduction = "cca")
object <- AddMetaData(object, metadata = c(cl1,cl2), col.name = "clusters")
object <- AddMetaData(object, metadata = c(rep(1,length(cl1)),rep(2,length(cl2))), col.name = "experiment")

Idents(object) <- object@meta.data$experiment
pdf(file.path(OUTDIR, "umap_CCAcorrected.pdf"))
print(DimPlot(object))
dev.off()

Idents(object) <- object@meta.data$sample.name
pdf(file.path(OUTDIR, "umap_CCAcorrected_bysample.pdf"), width = 14)
print(DimPlot(object, split.by = "experiment"))
dev.off()

Idents(object) <- object@meta.data$clusters
pdf(file.path(OUTDIR, "umap_CCAcorrected_bycluster.pdf"), width = 14)
print(DimPlot(object, split.by = "experiment"))
dev.off()

objfile <- file.path(OUTDIR, "object_cca.Rds")
saveRDS(object, file = objfile)

# extract common clones
expr.GBC.list1 <- as.character(obj1@meta.data$expr.GBC.list)
expr.GBC.list2 <- as.character(obj2@meta.data$expr.GBC.list)
clones <- unique(expr.GBC.list1[expr.GBC.list1 %in% expr.GBC.list2])

# compute distances in the CCA space
data <- as.matrix(object@reductions$cca@cell.embeddings[,1:nCCs])

# extract inter-sample distances
intra_clone_dist <- NULL
inter_clone_dist <- NULL
set.seed(SEED)
idx_rand <- sample(ncells2)
for (cl in clones) {
	idx1 <- which(expr.GBC.list1 == cl)
	x <- which(expr.GBC.list2 == cl)
	idx2 <- ncells1 + x
	idx2_rand <- ncells1 + idx_rand[x] 
	distances <- as.matrix(dist(data[c(idx1,idx2),], method = DISTANCE))
	dm <- distances[1:length(idx1),length(idx1)+1:length(idx2)]
	intra_clone_dist <- c(intra_clone_dist, drop(dm))
	distances <- as.matrix(dist(data[c(idx1,idx2_rand),], method = DISTANCE))
	dm <- distances[1:length(idx1),length(idx1)+1:length(idx2)]
	inter_clone_dist <- c(inter_clone_dist, drop(dm))
}
max_dist <- max(c(intra_clone_dist,inter_clone_dist))

# plot density
plotfile <- file.path(OUTDIR, paste0("clone_dist_density.pdf"))
pdf(plotfile, width = 5, height = 5)
d_intra <- density(intra_clone_dist)
d_inter <- density(inter_clone_dist)
plot(d_inter, lty = 2, lwd = 2, xlim = c(0,max_dist), main = "Cell-cell distance", xlab = "euclidean distance", ylab = "density")
lines(d_intra, col = "blue", lwd = 2) 
legend("topright", legend = c("sisters","random pairs"), lty = c(1,2), lwd = 2, col = c("blue","black"), bty = "n")
dev.off()

	


q()

