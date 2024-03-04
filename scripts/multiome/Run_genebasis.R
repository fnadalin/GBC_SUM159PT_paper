
# Takes as input a Seurat object and outputs a geneBasis object
############################## OPTIONS MENU ###################################

library("optparse")

# specify our desired options in a list
# by default OptionParser will add an help option equivalent to 
# make_option(c("-h", "--help"), action="store_true", default=FALSE, 
#               help="Show this help message and exit")

cat("\n")
option_list <- list(
  make_option("--Seurat", type="character", 
              help="[REQUIRED] Seurat Object path"),
  make_option("--out_RDS", type="character",
              help="[OPTIONAL] Output file name, otherwise input file name will be used"),
  make_option("--n_genes", type="integer", default=100,
              help="[OPTIONAL] Maximum number of genes to search [default=%default]"))

############################# PARSE OPTIONS ###################################

opt <- parse_args(OptionParser(option_list=option_list))

if (is.null(opt$Seurat)){
  write("Option --Seurat is required\nTry --help for help", stderr())
  q()
} else {
  SEURAT <- opt$Seurat
}

if (is.null(opt$out_RDS)) {
  OUT_FILE <- opt$Seurat
}

if (!is.null(opt$n_genes)) {
  N_GENES - opt$n_genes
}

################################# EXECUTION ###################################
library(geneBasisR)
library(Seurat)
library(Signac)

object <- readRDS(SEURAT)

sce <- retain_informative_genes(as.SingleCellExperiment(object, assay = "RNA"))
genes <- gene_search(sce,n_genes_total = N_GENES)

Misc(object, slot="geneBasis") <- genes

saveRDS(object, file = OUT_FILE)
