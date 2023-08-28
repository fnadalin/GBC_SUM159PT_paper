
PADDING <- 100000

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
    cat("\nUsage: <wes> <sc_name> <sc_obj> <sc_samples> <features> <lin1_marks> <lin2_marks>\n")
    q()
}

WES_EXP_FILE <- args[1]
SC_EXP_NAME <- args[2]
SC_EXP_OBJ <- args[3]
SC_SAMPLES <- args[4]          
FEATURES_INFO <- args[5]
LINEAGE_1_MARKERS <- args[6]
LINEAGE_2_MARKERS <- args[7]
   
library("Seurat")         
library("GenomicRanges")
library("Matrix")
library("pROC")

features <- read.table(FEATURES_INFO, sep = "\t")
idx <- which(features[,3] == "Gene Expression")
features <- features[idx,]

# rename duplicated features
uniq_features <- unique(features[,2])
uniq_features_freq <- unlist(lapply(uniq_features, function(x) sum(features[,2] == x)))
idx <- which(uniq_features_freq > 1)
uniq_features <- uniq_features[idx]
uniq_features_freq <- uniq_features_freq[idx]
for (x in 1:length(uniq_features)) {
    idx <- which(features[,2] == uniq_features[x])
    features[idx,2] <- paste(uniq_features[x],1:uniq_features_freq[x],sep=".")
}

lineage1_markers <- read.table(LINEAGE_1_MARKERS)[,1]
lineage2_markers <- read.table(LINEAGE_2_MARKERS)[,1]

OUTDIR <- paste0("out_", SC_EXP_NAME, "-WES_intersect")
dir.create(OUTDIR, showWarnings = FALSE)

object <- readRDS(SC_EXP_OBJ)
n_features <- nrow(object)

# map genes to their positions on the genome and create the granges object

all_genes <- rownames(object)

idx_list <- lapply(features[,2], function(x) grep(paste0("^",x,"$"),all_genes))
seq <- rep("",n_features)
start <- rep(-1,n_features)
end <- rep(-1,n_features)
for (x in 1:length(idx_list)) {
    for (el in idx_list[[x]]) {
        seq[el] <- features[x,4]
        start[el] <- features[x,5]
        end[el] <- features[x,6]
    }
} 

idx <- which(seq != "" & start != -1 & end != -1)
seq <- seq[idx]
start <- start[idx]
end <- end[idx]
genes <- all_genes[idx]
n_genes <- length(genes)

gr_rna <- GRanges(seqnames = Rle(seq), 
                  ranges = IRanges(start, end = end, names = 1:n_genes), 
                  strand = rep("*",n_genes))

# create the granges object for CNVs

# select only amplified / depleted CNVs
data <- read.table(WES_EXP_FILE, sep = "\t", header = TRUE)
avg_logfc <- apply(data[,4:ncol(data)], 1, mean)
idx <- which(abs(avg_logfc) > 0.1)
data <- data[idx,]
n_cnv <- nrow(data)

cnv_desc <- paste0(data[,1],":",data[,2],"-",data[,3])

# padding
data[,2] <- data[,2]-PADDING
data[data[,2] < 0,2] <- rep(0,sum(data[,2] < 0))
data[,3] <- data[,3]+PADDING

gr <- GRanges(seqnames = data[,1],
              ranges = IRanges(data[,2], end = data[,3]), 
              strand = rep("*",nrow(data)))
                  
write.table(data[,4:ncol(data)], file = file.path(OUTDIR,"logfc_annotation.tsv"), 
            row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
   
lineage <- object@meta.data$lineage
idx_row <- match(genes,rownames(object))
idx_col <- which(!is.na(lineage))
 
M <- object@assays$RNA@data[idx_row,idx_col]
lineage <- lineage[idx_col]
sample.name <- object@meta.data$sample.name[idx_col]

samples <- unlist(strsplit(SC_SAMPLES, split = ","))

for (s in samples) {

    idx_sample <- which(sample.name == s)
    
    lineage_sub <- lineage[idx_sample]
    M_sub <- M[,idx_sample]
    
    lin <- sort(unique(lineage_sub))
    if (length(lin) < 2) {
        next
    }
        
    for (k in 1:n_cnv) {
    
        dir.create(paste0(OUTDIR,"/",k), showWarnings = FALSE)

        write(cnv_desc[k], file = paste0(OUTDIR,"/",k,"/cnv_desc.txt"))
        
        # extract the normalized expression for each group of genes overlapping a CNV
       
        ov_loc <- subsetByOverlaps(gr_rna, gr[k])
        idx1 <- as.numeric(names(ov_loc))
        if (sum(width(ov_loc)) == 0) {
            next
        }
        norm_counts <- M_sub[idx1,]
        
        is_lineage_marker <- rep(0,length(idx1))
        idx <- which(rownames(norm_counts) %in% lineage1_markers)
        is_lineage_marker[idx] <- rep(1,length(idx))
        idx <- which(rownames(norm_counts) %in% lineage2_markers)
        is_lineage_marker[idx] <- rep(2,length(idx))
      
        write.table(norm_counts, file = paste0(OUTDIR,"/",k,"/matrix_",s,".tsv"), sep = "\t", quote = FALSE)
        write.table(lineage_sub, file = paste0(OUTDIR,"/",k,"/lineage_annotation_",s,".txt"), 
                    row.names = FALSE, col.names = FALSE, quote = FALSE)
        write.table(is_lineage_marker, file = file.path(OUTDIR,"/",k,"is_lineage_marker.txt"), 
                    row.names = FALSE, col.names = FALSE, quote = FALSE)
        
        # compute AUC for each CNV and for each lineage
        auc <- matrix(0.5, nrow = nrow(norm_counts), ncol = length(lin))
        rownames(auc) <- rownames(norm_counts)
        colnames(auc) <- paste0("auc_lineage", lin)
        for (h in 1:nrow(norm_counts)) {
            for (l in 1:length(lin)) {
                df <- data.frame(a = norm_counts[h,], b = as.numeric(lineage_sub == lin[l]))
                if (min(df$a) < max(df$a)) {
                    rocres <- roc(data = df, response = "b", predictor = "a", direction = "<")
                    auc[h,l] <- rocres$auc
                }
            }
        }
        df_auc <- cbind(geneID = rownames(norm_counts), as.data.frame(gr_rna[idx1]), auc)
    
        write.table(df_auc, file = paste0(OUTDIR,"/",k,"/auc_",s,".txt"), row.names = FALSE, sep = "\t", quote = FALSE)
        
    } 

}


q()

