
EXP <- c("1C_ARC","1E_ARC")

library("Seurat")

for (exp in EXP) {
    
    OBJECT <- paste0("../../analysis/multiome/", exp, "/", exp, "_seurat.Rds")

    object <- readRDS(OBJECT)
    cistopic <- object[["cisTopic"]]
    M <- cistopic@misc$cisTopicObject@count.matrix
    bed <- gsub("-","	",gsub(":","	",rownames(M)))
    
    write(bed, file = paste0("cistopic_input_peaks_", exp, ".bed"))
}

q()


