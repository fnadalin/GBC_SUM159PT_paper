
ALPHA <- c(0.7, 0.8, 0.9)
CUTOFF <- 0.15 # select at most 15% cells

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
    cat("\nUsage: <sc_dataset> <id> <out_prefix>\n")
    q()
}

OBJECT <- args[1]
IN_PREFIX <- args[2]
ID <- args[3]
OUT_PREFIX <- args[4]

library("Scissor")
library("Seurat")

PHENOTYPE <- paste0(IN_PREFIX, "_phenotype.tsv")
BULK_DATASET <- paste0(IN_PREFIX, "_geneExprTable.tsv")

sc_dataset <- readRDS(OBJECT)

phenotype <- read.table(PHENOTYPE, sep = "\t", header = TRUE)
bulk_dataset <- read.table(BULK_DATASET, sep = "\t", header = TRUE)

tag <- c(paste0("in-",ID), paste0("outside-",ID))

infos <- Scissor(as.matrix(bulk_dataset), sc_dataset, as.matrix(phenotype), tag = tag, 
                  alpha = ALPHA, family = "binomial", cutoff = CUTOFF,
                  Save_file = file.path(OUT_PREFIX,"scissor_in.RData"))   
save(infos, file = file.path(OUT_PREFIX,"scissor_out.RData"))

q()

