
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
    cat("\nUsage: <obj> <sc_obj> <sc_out_dir> <pops>\n")
    q()
}

OBJECT <- args[1]
SC_OBJECT <- args[2]
DIR <- args[3]
POP <- unlist(strsplit(args[4], split = ","))

library("Scissor")
library("Seurat")

object <- readRDS(OBJECT)
sc_dataset <- readRDS(SC_OBJECT)

for (s in POP) {

    SCISSOR <- file.path(DIR,s,"/scissor_out.RData")
    META <- paste0("scissor.",s)
    
    load(SCISSOR)
    
    Scissor_select <- rep(0, ncol(sc_dataset))
    names(Scissor_select) <- colnames(sc_dataset)

    Scissor_select[infos$Scissor_pos] <- 1
    Scissor_select[infos$Scissor_neg] <- 2
    sc_dataset <- AddMetaData(sc_dataset, metadata = Scissor_select, col.name = META)

}

sample_name <- object@meta.data$sample.name
df <- data.frame(sample.name = sample_name)
rownames(df) <- colnames(sc_dataset)
sc_dataset <- AddMetaData(sc_dataset, metadata = df)

saveRDS(sc_dataset, SC_OBJECT)

q()

