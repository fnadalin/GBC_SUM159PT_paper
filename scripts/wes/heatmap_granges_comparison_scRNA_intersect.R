
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
    cat("\nUsage: <sc_name> <sc_samples>\n")
    q()
}

SC_EXP_NAME <- args[1]
SC_SAMPLES <- args[2]
 
library("ComplexHeatmap")
library("circlize")
library("scales")
library("ggrepel")
library("RColorBrewer")

OUTDIR <- paste0("out_", SC_EXP_NAME, "-WES_intersect")

cnv_ann <- read.table(file.path(OUTDIR,"logfc_annotation.tsv"), sep = "\t")
n_cnv <- nrow(cnv_ann)
samples <- unlist(strsplit(SC_SAMPLES, split = ","))

for (s in samples) {
        
    for (k in 1:n_cnv) {
    
        cnv_desc <- read.table(paste0(OUTDIR,"/",k,"/cnv_desc.txt"), sep = "\t")[,1]
        norm_counts <- read.table(paste0(OUTDIR,"/",k,"/matrix_",s,".tsv"), sep = "\t", header = TRUE)
        lineage <- read.table(paste0(OUTDIR,"/",k,"/lineage_annotation_",s,".txt"))[,1]
        is_lineage_marker <- read.table(file.path(OUTDIR,"/",k,"is_lineage_marker.txt"))[,1]
        
        NORM_COUNT_HEATMAP <- paste0(OUTDIR,"/",k,"/heatmap_matrix_",s,".pdf")
        
        avg_exp1 <- apply(norm_counts[,lineage == 1], 1, function(x) mean(x = expm1(x = x)))
        avg_exp2 <- apply(norm_counts[,lineage == 2], 1, function(x) mean(x = expm1(x = x)))
        log2fc <- log2(avg_exp1) - log2(avg_exp2)
        M <- t(scale(t(norm_counts)))
        
        idx <- which(rowSums(is.na(M))==0)
        M <- M[idx,]
        log2fc <- log2fc[idx]
        avg_exp1 <- avg_exp1[idx]
        is_lineage_marker <- is_lineage_marker[idx]
        
        idx2 <- order(avg_exp1, decreasing = TRUE) 
        M <- M[idx2,]
        is_lineage_marker <- is_lineage_marker[idx2]
        
        idx3 <- order(lineage)
        M <- M[,idx3]
        lineage <- lineage[idx3]

        # select palette
     
        col_fun_lin <- hue_pal()(2)
        names(col_fun_lin) <- sort(unique(lineage))
        col_fun <- colorRamp2(c(-1,0,1), c("white","gray90","black"))
        
        gene_labels_id <- which(is_lineage_marker == "1")
        genes <- rownames(M)[gene_labels_id]
        
        gene_labels <- rowAnnotation(highlight = anno_mark(at = gene_labels_id, labels = genes, 
                         which = "row", side = 'right')) # labels_gp = gpar(fontsize = 9)))
        hb <- HeatmapAnnotation(lineage = lineage, which = "column", col = list(lineage = col_fun_lin),
                        show_legend = TRUE, show_annotation_name = TRUE)
        h <- Heatmap(M, name = "Z-score", col = col_fun, top_annotation = hb, 
                     cluster_rows = FALSE, show_column_names = FALSE, show_row_names = FALSE,
                     column_title = cnv_desc, column_split = lineage, right_annotation = gene_labels, 
                     show_row_dend = FALSE, show_column_dend = FALSE, row_title = "gene loci spanned by CNV")
                     
        pdf(NORM_COUNT_HEATMAP, width = 5, height = 4.5)
        draw(h)
        dev.off()
                    
    } 

}

q()

