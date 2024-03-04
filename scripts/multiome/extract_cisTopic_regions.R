# Take a Seurat object and perform cisTopic analysis on it, 
# extract the cell-topics matrix and attaches it to the original seurat object as a reduction
# the full cistopic object is stored in the misc slot of the reduction of the seurat object

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
              help="[REQUIRED] Seurat Object path"))

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

################################# EXECUTION ###################################

library(Seurat)
library(Signac)
library(ggplot2)
library(cisTopic)
library(EnsDb.Hsapiens.v86)

# we load the Seurat object of our experiment
object <- readRDS(SEURAT)

cisTopicObject <- object@reductions$cisTopic@misc$cisTopicObject

cisTopicObject <- getRegionsScores(cisTopicObject, method = "NormTop", scaled = TRUE)

cisTopicObject <- binarizecisTopics(cisTopicObject)

bed_path <- file.path(dirname(SEURAT), paste(sub(".Rds","",basename(SEURAT)), "cisTopic", "bedfiles", sep="_"))

getBedFiles(object = cisTopicObject, path = bed_path)

object@reductions$cisTopic@misc$cisTopicObject <- cisTopicObject

saveRDS(object, file = SEURAT)
q()
