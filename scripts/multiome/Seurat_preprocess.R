
# Take a Seurat object as input and performs subsetting, normalization, scaling

SUBSET <- NULL
set.seed(1234)
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
  make_option("--subset", type="character", default=NULL,
              help="[OPTIONAL] logical expression with the subset of cells 
              to filter")
)

############################## PARSE OPTIONS ##################################

# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults
opt <- parse_args(OptionParser(option_list=option_list))

if (is.null(opt$Seurat)) {
  write("Option --Seurat is required\nTry --help for help", stderr())
  q()
} else {
  SEURAT <- opt$Seurat
}

if (!is.null(opt$subset)) {
  SUBSET <- opt$subset
}


############################## EXECUTION ######################################

library(Seurat)
library(Signac)
library(ggplot2)

# we load the object
object <- readRDS(SEURAT)
EXP_NAME <- sub(".rds", "", basename(SEURAT), ignore.case = T)
# we check if the subset logical expression is empty, then proceed with the 
# filtering
if (!is.null(SUBSET)) {
  object <- subset(
    x = object,
    subset = str2lang(SUBSET)
  )
}

## RNA preprocessing

DefaultAssay(object) <- "RNA"

object <- ScaleData(object)
object <- RunPCA(object, npcs = 30, reduction.name="pca_rna")
object <- RunUMAP(object, reduction="pca_rna", dims = 1:30, reduction.name = "umap_rna")

## ATAC preprocessing

DefaultAssay(object) <- "ATAC"

object <- FindTopFeatures(object, min.cutoff = 10)
object <- RunSVD(object, reduction.name = "lsi_atac")
object <- RunUMAP(object, reduction = 'lsi_atac', dims = 2:50, reduction.name = 'umap_atac')

dir.create("plots", recursive = T)

DepthCor(object, reduction="lsi_atac")
ggsave(paste(paste(dirname(SEURAT),"/plots/", sep=''), EXP_NAME,"_LSI_depth_correlation.png", sep=''), width = 10, height = 5)

saveRDS(object, file = SEURAT)

sessionInfo()
q()




