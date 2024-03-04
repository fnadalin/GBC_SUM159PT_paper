# FRAC_DETECTED <- 0.01
LABELS_HIGHLIGHT <- 1:24 # only highlight the programs common to more than one tumour

# FIXME: only works for 3 cell subpopulations!!!

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 7) {
    cat("\nUsage: <table3CA> <object> <cl_mode> <subpop_id> <subpop_name> <out_dir> <sign_dir>\n")
    q()
}

TABLE_3CA <- args[1]
OBJECT <- args[2]
CL_MODE <- args[3]
SUBPOP_ID <- unlist(strsplit(args[4], split = ","))
SUBPOP_NAME <- unlist(strsplit(args[5], split = ","))
OUT_DIR <- args[6]
SIGN_DIR <- args[7]

library("Seurat")
library("ROCR")
library("ComplexHeatmap")
library("circlize")
library("scales")
library("RColorBrewer")

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

df_3CA <- read.table(TABLE_3CA, sep = "\t", header = TRUE, check.names = FALSE)
object <- readRDS(OBJECT)

sign1 <- read.table(file.path(SIGN_DIR, "signature1.txt"))[,1]
sign2 <- read.table(file.path(SIGN_DIR, "signature2.txt"))[,1]
sign3 <- read.table(file.path(SIGN_DIR, "signature3.txt"))[,1]

object <- subset(object, cells = colnames(object)[!is.na(object@meta.data[[CL_MODE]])])
CL_ID <- as.character(unique(object@meta.data[[CL_MODE]]))
ordering <- c(SUBPOP_ID, sort(setdiff(CL_ID,SUBPOP_ID)))
CL_ID <- CL_ID[match(ordering, CL_ID)]
CL_NAME <- CL_ID
CL_NAME[match(SUBPOP_ID, CL_ID)] <- SUBPOP_NAME

# select genes detected in at least 1% of the cells
# cc <- rowSums(object@assays$RNA@counts > 0)
# g <- rownames(object)[cc > round(FRAC_DETECTED*ncol(object))]
# object <- subset(object, features = g)
Idents(object) <- object@meta.data[[CL_MODE]]

# strategy 1: for each subpopulation S and for each meta-program M, use the log2FC(S,complement(S)) 
#             of gene g to predict whether g belongs to M

M <- matrix(NA, nrow = length(CL_ID), ncol = ncol(df_3CA))
colnames(M) <- colnames(df_3CA)
rownames(M) <- CL_NAME

for (i in 1:length(CL_ID)) {
    logFC <- FoldChange(object, ident.1 = CL_ID[i])
    M[i,] <- apply(df_3CA, 2, function(x) { labels = as.numeric(rownames(logFC) %in% x)
                                            pred <- prediction(predictions = logFC$avg_log2FC, labels = labels)
                                            auc <- performance(pred, measure = "auc")
                                            auc@y.values[[1]] })
}

write.table(M, file = file.path(OUT_DIR,"auc_bygene.tsv"), sep = "\t", quote = FALSE)

pdf(file.path(OUT_DIR,"auc_bygene.pdf"), width = 5, height = 7)
fontface <- rep("plain",ncol(M))
fontface[LABELS_HIGHLIGHT] <- rep("bold",length(LABELS_HIGHLIGHT))
Heatmap(t(M), name = "auc_byGene", cluster_columns = FALSE, row_names_gp = gpar(fontface = fontface))
dev.off()

# strategy 2: for each gene g in meta-program M and for each subpopulation S, use the normalised expression 
#             of g in cell c to predict whether c belong to S

name <- "MP"
object <- AddModuleScore(object, features = as.list(df_3CA), search = FALSE, name = "MP") # search = FALSE because the synomyms have been replaced already

M <- matrix(NA, nrow = length(CL_ID), ncol = ncol(df_3CA))
colnames(M) <- colnames(df_3CA)
rownames(M) <- CL_NAME

for (i in 1:length(CL_ID)) {
    M[i,] <- sapply(1:ncol(df_3CA), function(x) { labels = as.numeric(Idents(object) == CL_ID[i])
                                                  pred <- prediction(predictions = object@meta.data[[paste0(name,x)]], labels = labels)
                                                  auc <- performance(pred, measure = "auc")
                                                  auc@y.values[[1]] }) 
}

write.table(M, file = file.path(OUT_DIR,"auc_bycell.tsv"), sep = "\t", quote = FALSE)

pdf(file.path(OUT_DIR,"auc_bycell.pdf"), width = 5, height = 7)
fontface <- rep("plain",ncol(M))
fontface[LABELS_HIGHLIGHT] <- rep("bold",length(LABELS_HIGHLIGHT))
Heatmap(t(M), name = "auc_byCell", cluster_columns = FALSE, row_names_gp = gpar(fontface = fontface))
dev.off()

# plot heatmaps of meta-programs

cols <- c(brewer.pal(n = 3, name = "Set1"), "gray20", "white")
names(cols) <- c("S1","S2","S3","shared","not found")

for (j in 1:ncol(df_3CA)) {
    v <- df_3CA[,j]
    m <- object@assays$RNA@scale.data[v[!is.na(v) & v != ""],]
    
    bins <- seq(from = 0, to = ncol(object), length.out = 101)
    bins <- round(bins)
    
    idx_cells <- order(object@meta.data[[paste0(name,j)]])
    w <- as.numeric(object@meta.data[[CL_MODE]] == 6)[idx_cells]
    bin_count_S1 <- lapply(2:length(bins), function(x) rep(sum(w[bins[x-1]:bins[x]]),bins[x]-bins[x-1]))
    w <- as.numeric(object@meta.data[[CL_MODE]] == 2)[idx_cells]
    bin_count_S2 <- lapply(2:length(bins), function(x) rep(sum(w[bins[x-1]:bins[x]]),bins[x]-bins[x-1]))
    w <- as.numeric(object@meta.data[[CL_MODE]] == 5)[idx_cells]
    bin_count_S3 <- lapply(2:length(bins), function(x) rep(sum(w[bins[x-1]:bins[x]]),bins[x]-bins[x-1]))
    
    by <- length(w)/100
    S1_bar <- anno_barplot(unlist(bin_count_S1), height = unit(1, "cm"))
    S2_bar <- anno_barplot(unlist(bin_count_S2), height = unit(1, "cm"))
    S3_bar <- anno_barplot(unlist(bin_count_S3), height = unit(1, "cm"))

    # map the signatures

    ann <- rep("not found", nrow(m))
    idx3 <- which(rownames(m) %in% sign3)
    ann[idx3] <- rep("S3", length(idx3))
    idx2 <- which(rownames(m) %in% sign2)
    ann[idx2] <- rep("S2", length(idx2))
    idx1 <- which(rownames(m) %in% sign1)
    ann[idx1] <- rep("S1", length(idx1))
    idx_genes <- sort(unique(c(idx1,idx2,idx3)))

    v <- c(idx1,idx2,idx3)
    w <- match(idx_genes,v)
    shared <- v[-w]
    ann[shared] <- rep("shared", length(shared))
    
    gene_labels <- rownames(m)[idx_genes]
    
    ht <- columnAnnotation(S1 = S1_bar, S2 = S2_bar, S3 = S3_bar)
    hr <- rowAnnotation(signature = ann, col = list(signature = cols), show_legend = TRUE)
    anno <- rowAnnotation(highlight = anno_mark(at = idx_genes, labels = gene_labels))
    title <- colnames(df_3CA)[j]
#    h <- Heatmap(m[,idx_cells], name = "z-score", top_annotation = ht, right_annotation = hr, 
#                 row_names_gp = gpar(fontsize = 8), column_title = title, show_column_names = FALSE, cluster_columns = FALSE)
    h <- Heatmap(m[,idx_cells], name = "z-score", top_annotation = ht, left_annotation = hr, right_annotation = anno,
                 column_title = title, show_column_names = FALSE, show_row_names = FALSE, cluster_columns = FALSE)
    
    title <- gsub("/", "-", title)
    title <- gsub("\\(|\\)", "", title)
    title <- gsub(" ", "_", title)
    pdf(file.path(OUT_DIR, paste0(title, ".pdf")), width = 5, height = 3.5)
    draw(h)
    dev.off()
}



q()

