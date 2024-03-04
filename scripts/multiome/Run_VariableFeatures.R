
# Takes a Seurat object as input, normalizes the data and calculates the variable features for each modality


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
  make_option("--RNA", type="integer", default=3000,
              help="[OPTIONAL] Number of RNA features to include in the set [default=%default]"),
  make_option("--ATAC", type="integer", default=5000,
              help="[OPTIONAL] Number of ATAC features to include in the set [default=%default]"))


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

if (!is.null(opt$RNA)) {
  RNA <- opt$RNA
}

if (!is.null(opt$ATAC)) {
  ATAC <- opt$ATAC
}

################################# EXECUTION ###################################

library(Seurat)
library(Signac)

object <- readRDS(SEURAT)

object <- NormalizeData(object, assay = "RNA")
object <- FindVariableFeatures(object = object, 
                               assay = "RNA", 
                               nfeatures = RNA)

cat("\nCalculated Variable Features for RNA\n")

object <- RunTFIDF(object, assay = "ATAC")
object <- FindVariableFeatures(object = object, 
                               assay = "ATAC",
                               selection.method = "disp",
                               nfeatures = ATAC)

saveRDS(object, file = OUT_FILE)
sessionInfo()
q()