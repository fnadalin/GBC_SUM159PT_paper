
ComputeStats <- function(cl, expr.GBC.list) {

    df <- as.data.frame(matrix(NA, nrow = length(unique(cl)), ncol = 6))
    colnames(df) <-  c("cluster","clones","cells.in.cluster","cells.in.clones","clones.only.in.cluster","cells.in.clones.only.in.cluster")
    df[,1] <- unique(cl)

    ncl <- length(unique(cl))
    i <- 1
    for (c in unique(cl)) {
        idx1 <- which(cl == c)
        df[i,3] <- length(idx1)
        clus <- unique(expr.GBC.list[idx1])
        df[i,2] <- length(clus)
        idx2 <- which(expr.GBC.list %in% clus)
        df[i,4] <- length(idx2)
        idx3 <- idx2[!(idx2 %in% idx1)] # cells outside the cluster with the label of a clone inside the cluster
        clones_out <- unique(expr.GBC.list[idx3])
        clones_in <- clus[!(clus %in% clones_out)]
        df[i,5] <- length(clones_in)
        idx4 <- which(cl == clones_in)
        df[i,6] <- length(idx4)
        i <- i + 1
    }
    
    return(df)
}


args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
        cat("\nUsage: Rscript cluster_clone_pairs.R <obj1> <obj2> <cl1> <cl2> <out_dir>\n")
        q()
}

OBJ1 <- args[1]
OBJ2 <- args[2]
CL1 <- args[3]
CL2 <- args[4]
OUTDIR <- args[5]

dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

library("Seurat")

# cells in objects have been filtered for doublets and low quality cells
# they also contain relevant metadata
obj1 <- readRDS(OBJ1)
obj2 <- readRDS(OBJ2)

# extract common clones
expr.GBC.list1 <- as.character(obj1@meta.data$expr.GBC.list)
expr.GBC.list2 <- as.character(obj2@meta.data$expr.GBC.list)
clones <- unique(expr.GBC.list1[expr.GBC.list1 %in% expr.GBC.list2])

# subset to common clones
cells1 <- which(expr.GBC.list1 %in% clones)
cells2 <- which(expr.GBC.list2 %in% clones)
obj1 <- subset(obj1, cells = cells1)
obj2 <- subset(obj2, cells = cells2)
expr.GBC.list1 <- as.character(obj1@meta.data$expr.GBC.list)
expr.GBC.list2 <- as.character(obj2@meta.data$expr.GBC.list)
ncells1 <- ncol(obj1)
ncells2 <- ncol(obj2)

# extract cluster labels
cl1 <- obj1@meta.data[[CL1]]
cl2 <- obj2@meta.data[[CL2]]
ncl1 <- length(unique(cl1))
ncl2 <- length(unique(cl2))

# compute cluster assignment
M1 <- matrix(0, nrow = ncl1, ncol = ncl2)
M2 <- matrix(0, nrow = ncl1, ncol = ncl2)
x <- 1:nrow(M1)
y <- 1:ncol(M1)
rownames(M1) <- rownames(M2) <- x-1
colnames(M1) <- colnames(M2) <- y-1
for (cl in clones) {
        idx1 <- which(expr.GBC.list1 == cl)
        idx2 <- which(expr.GBC.list2 == cl)
        n1 <- as.numeric(table(cl1[idx1]))
        n2 <- as.numeric(table(cl2[idx2]))
        # count the cells in cl1 having at least one sister in cl2
        M1[n1>0,n2>0] <- M1[n1>0,n2>0] + n1[n1>0] 
        # count the cells in cl2 having at least one sister in cl1
        M2[n1>0,n2>0] <- t(t(M2[n1>0,n2>0]) + n2[n2>0])
}
M <- M1 * M2 # cells with sisters in obj1 * cells with sisters in obj2

# normalize
cl1_size <- as.numeric(table(cl1))
cl2_size <- as.numeric(table(cl2))
M_norm <- t(t(M/cl1_size)/cl2_size) # compute the fraction of cells with sisters in both objects

# print matrix
write.table(M, file = file.path(OUTDIR, "cluster_pairs.tsv"), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(M_norm, file = file.path(OUTDIR, "cluster_pairs_norm.tsv"), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

# print clone lists for each cluster pair
for (i in unique(cl1)) {
	for (j in unique(cl2)) {
		x1 <- expr.GBC.list1[cl1 == i]
		x2 <- expr.GBC.list2[cl2 == j]
		cl <- x1[x1 %in% x2] # clones in cluster i and j
		cl <- unique(cl)
		cl <- cl[cl!=""]
		write(cl, file = file.path(OUTDIR, paste0("GBC_list_",i,"-",j,".txt"))) 
	}
}

# compute stats:
# count the clones per cluster and the number of cells for those clones (both within clusters and in the whole sample)

df1 <- ComputeStats(cl1, expr.GBC.list1)
df2 <- ComputeStats(cl2, expr.GBC.list2)

# print
write.table(df1, file = file.path(OUTDIR, "sample1_stats.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(df2, file = file.path(OUTDIR, "sample2_stats.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)


q()


