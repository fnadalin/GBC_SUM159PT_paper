# FRAC_DETECTED <- 0.01

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
    cat("\nUsage: <table3CA> <object> <out_dir>\n")
    q()
}

TABLE_3CA <- args[1]
OBJECT <- args[2]
OUT_DIR <- args[3]

library("Seurat")

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

df_3CA <- read.table(TABLE_3CA, sep = "\t")
colnames(df_3CA) <- df_3CA[1,]
df_3CA <- df_3CA[2:nrow(df_3CA),]

object <- readRDS(OBJECT)

# select genes detected in at least 1% of the cells
# cc <- rowSums(object@assays$RNA@counts > 0)
# g <- rownames(object)[cc > round(FRAC_DETECTED*ncol(object))]
# object <- subset(object, features = g)

# find current gene names by synonims
for (i in 1:ncol(df_3CA)) {
    missing_genes <- setdiff(df_3CA[,i], rownames(object))
    if (length(missing_genes) > 0) {
        mp_genes <- intersect(df_3CA[,i], rownames(object))
        missing_genes <- UpdateSymbolList(missing_genes)
        to_add <- intersect(missing_genes, rownames(object))
        if (length(to_add) > 0) {
            mp_genes <- c(mp_genes, to_add)
        }
        if (length(mp_genes) < nrow(df_3CA)) { 
            mp_genes <- c(mp_genes, rep("", nrow(df_3CA)-length(mp_genes)))
        }
        df_3CA[,i] <- mp_genes
    }
}

write.table(df_3CA, file = file.path(OUT_DIR, "MP_table.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

q()

