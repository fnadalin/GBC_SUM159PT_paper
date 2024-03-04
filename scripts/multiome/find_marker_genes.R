set.seed(1234)
ALGORITHM <- 3 # SLM algorithm is set as default
RESOLUTION <- 0.5 # clustering resolution is set to 0.5
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
  make_option("--clustering_name", type="character",
              help="[REQUIRED] name of the metadata column of the clusters"),
  make_option("--assay", type="character", default="RNA",
              help="[OPTIONAL] Specifiy assay to select features (genes or regions)"),
  make_option("--test_use", default="wilcox",
              help="[OPTIONAL] name of the statistical test to be used"))

############################# PARSE OPTIONS ###################################

opt <- parse_args(OptionParser(option_list=option_list))

if (is.null(opt$Seurat)){
  write("Option --Seurat is required\nTry --help for help", stderr())
  q()
} else {
  SEURAT <- opt$Seurat
}

if (is.null(opt$clustering_name)){
  write("Option --clustering_name is required\nTry --help for help", stderr())
  q()
} else {
  CLUSTERING_NAME <- opt$clustering_name
}

ASSAY <- opt$assay

TEST_USE <- opt$test_use

############################### EXECUTION #####################################
library(Seurat)
library(Signac)
library(ggplot2)

object <- readRDS(SEURAT)

Idents(object) <- object@meta.data[[CLUSTERING_NAME]]

if (TEST_USE == "LR") {
  df <- FindAllMarkers(object, 
                       assay = ASSAY, 
                       test.use = TEST_USE, 
                       latent.vars = "frequency_count",
                       min.pct = 0.05)
} else {
  df <- FindAllMarkers(object, assay = ASSAY, test.use = TEST_USE)
}


filename <- paste(sub(".Rds","",basename(SEURAT), ignore.case = T),CLUSTERING_NAME,TEST_USE,ASSAY,"marker_features.csv", sep='_')

write.csv2(df, file = file.path(dirname(SEURAT), filename))

q()