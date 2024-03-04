
# Take a Seurat object and perform cisTopic analysis on it, 
# extract the cell-topics matrix and attaches it to the original seurat object as a reduction
# the full cistopic object is stored in the misc slot of the reduction of the seurat object

set.seed(1234)
read.as.input <- function(input_text) {
  input_text = paste0("c(", input_text, ")")
  result = eval(parse(text = input_text))
  return(result)
}

localMaxima <- function(x) {
  # Use -Inf instead if x is numeric (non-integer)
  y <- diff(c(-.Machine$integer.max, x)) > 0L
  rle(y)$lengths
  y <- cumsum(rle(y)$lengths)
  y <- y[seq.int(1L, length(y), 2L)]
  if (x[[1]] == x[[2]]) {
    y <- y[-1]
  }
  y
}
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
  make_option("--Seurat2", type="character", 
              help="[REQUIRED] Seurat Object path with which select overlapping regions only"),
  make_option("--out_RDS", type="character",
              help="[OPTIONAL] Output file name, otherwise input file name will be used"),
  make_option("--topics", type="character", default="20:50",
              help="[OPTIONAL] Range and/or comma-separated of topics number to test [start:end=%default]"),
  make_option("--n_cpus", type="integer", default=1,
              help="[OPTIONAL] Number of CPUs to be utilized during model training"),
  make_option("--min_cutoff", type="character", default=-1,
              help="[OPTIONAL] Minimum cutoff for regions to be included in the cisTopic analysis (available as variable features, input 'none' if Variable Features)"),
  make_option("--n_iter", type="integer", default=500,
              help="[OPTIONAL] Total number of iterations for the algorithm"),
  make_option("--n_burnin", type="integer", default=250,
              help="[OPTIONAL] Number of burn-in iterations"),
  make_option("--var_feat", type="character", default="all",
              help="[OPTIONAL] Number of variable features to be included in the cisTopic analysis (available as variable features)"),
  make_option("--clean_regions", type="logical", default=TRUE,
              help="[OPTIONAL] Removes the lowly accessible regions included when running feature selection")
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

if (is.null(opt$Seurat2)) {
  write("Option --Seurat2 is required\nTry --help for help", stderr())
  q()
} else {
  SEURAT2 <- opt$Seurat2
}

if (is.null(opt$out_RDS)) {
  OUT_FILE <- opt$out_RDS
} else {
  OUT_FILE <- opt$Seurat
}

if (!is.null(opt$topics)) {
  TOPICS <- read.as.input(opt$topics)
} else {
  TOPICS <- 20:50
}

NCORES <- opt$n_cpus
if (!is.na(as.numeric(opt$min_cutoff))) {
  # if it is a percentage or a minimum cutoff of cells
  MIN.CUTOFF <- as.numeric(opt$min_cutoff)
} else {
  # if it is a quantile
  MIN.CUTOFF <- opt$min_cutoff
}

ITERATIONS <- opt$n_iter
BURNIN <- opt$n_burnin
VARFEAT <- opt$var_feat
CLEAN_REGIONS <- opt$clean_regions
################################# EXECUTION ###################################

suppressPackageStartupMessages({library(Seurat)
library(Signac)
library(ggplot2)
library(cisTopic)
library(EnsDb.Hsapiens.v86)})
cat("##--------------------------- cisTopic training ---------------------------------\n##\n")
cat("## FILE READING:\n##\n")

plots_dir <- paste(paste(dirname(SEURAT),"plots","",sep='/'))

if (!dir.exists(plots_dir)){
  dir.create(plots_dir, recursive = T)
}

# we load the Seurat object of our experiment
object <- readRDS(SEURAT)
object2 <- readRDS(SEURAT2)

atac.ranges.1 <- GetAssayData(object, assay = "ATAC", slot="ranges")

# we don't need to select variable features in Object 2 because we use it only to 
# match regions, independently of their importance
atac.ranges.2 <- GetAssayData(object2, assay = "ATAC", slot="ranges")

overlaps <- findOverlaps(query=atac.ranges.1, subject=atac.ranges.2)

object.filt <- subset(object, 
                      features=rownames(object[["ATAC"]])[unique(queryHits(overlaps))])

object.filt.2 <- subset(object2, 
                      features=rownames(object2[["ATAC"]])[unique(subjectHits(overlaps))])

# we clean the dataset from all possible completely inactive regions
object.filt <- FindTopFeatures(object.filt, assay = "ATAC", min.cutoff = MIN.CUTOFF)
object.filt.2 <- FindTopFeatures(object.filt.2, assay = "ATAC", min.cutoff = MIN.CUTOFF)

if (VARFEAT == "all") {
  VariableFeatures(object.filt) <- rownames(object.filt@assays$ATAC)
} else{
  par(mfrow=c(1,2))
  var.feat <- as.integer(VARFEAT)
  object.filt <- FindVariableFeatures(object.filt,
                                 assay="ATAC",
                                 selection.method="disp",
                                 nfeatures = var.feat)
  if (CLEAN_REGIONS) {
    active.features <- rowSums(object.filt[["ATAC"]]@counts[VariableFeatures(object.filt[["ATAC"]]),] > 0)
    active.features.hist <- hist(active.features, breaks=100)
    active.features.maxima <- localMaxima(active.features.hist$counts)
    active.features.minimum <- active.features.hist$breaks[active.features.maxima[1] + which.min(active.features.hist$counts[active.features.maxima[1]:active.features.maxima[2]]) - 1]
    abline(v=active.features.minimum, col="red", lty=2, lwd=2)
    VariableFeatures(object.filt) <- VariableFeatures(object.filt)[which(active.features > active.features.minimum)]
  }
  
  object.filt.2 <- FindVariableFeatures(object.filt.2,
                                      assay="ATAC",
                                      selection.method="disp",
                                      nfeatures = var.feat)
  if (CLEAN_REGIONS) {
    active.features <- rowSums(object.filt.2[["ATAC"]]@counts[VariableFeatures(object.filt.2[["ATAC"]]),] > 0)
    active.features.hist <- hist(active.features, breaks=100)
    active.features.maxima <- localMaxima(active.features.hist$counts)
    active.features.minimum <- active.features.hist$breaks[active.features.maxima[1] + which.min(active.features.hist$counts[active.features.maxima[1]:active.features.maxima[2]]) - 1]
    abline(v=active.features.minimum, col="red", lty=2, lwd=2)
    VariableFeatures(object.filt.2) <- VariableFeatures(object.filt.2)[which(active.features > active.features.minimum)]
  }
  # we get a vector of length Variablefeatures(object.filt) with 1s where the regions overlap and 0 otherwise
  overlapping.var.feats <- findOverlaps(StringToGRanges(VariableFeatures(object.filt)), 
                                         StringToGRanges(VariableFeatures(object.filt.2)))
  
  VariableFeatures(object.filt) <- VariableFeatures(object.filt)[unique(queryHits(overlapping.var.feats))]
}
cat("## cisTopic matrix preparation\n##\n")
# retaining, among the overlapping regions, the unique occurrence of it
cistopic.count.matrix <- object.filt[["ATAC"]]@counts[VariableFeatures(object.filt),]
cat(paste("Matrix dimension:", dim(cistopic.count.matrix), collapse = " "))
# we reformat annotations for a bug in cisTopic
rownames(cistopic.count.matrix) <- sub('-',':', 
                                       rownames(cistopic.count.matrix))

# we create the cisTopic object
cisTopicObject <- createcisTopicObject(
  count.matrix = cistopic.count.matrix
)
cat()
cat("## cisTopic run:\n##\n")
# Collapsed Gibbs Sampler model training

cisTopicObject <- runCGSModels(
  cisTopicObject,
  topic = TOPICS,
  nCores = NCORES,
  iterations = ITERATIONS, 
  burnin = BURNIN
)

saveRDS(cisTopicObject, file.path(dirname(SEURAT), paste0(sub(".Rds","",basename(SEURAT)),"_cisTopicObject.Rds")))

png(filename=file.path(dirname(SEURAT), 
                       "plots", 
                       paste(sub(".Rds","",basename(SEURAT)), "cisTopic_loglikelihood_plot.png", sep="_")),
    height=16, width = 20, res=300, units="cm")
logLikelihoodByIter(cisTopicObject)
dev.off()

cat("## Topic model selection\n##\n")
png(filename=file.path(dirname(SEURAT), 
                       "plots", 
                       paste(sub(".Rds","",basename(SEURAT)), "cisTopic_model_selection.png", sep="_")),
    height=8, width = 20, res=300, units="cm")
par(mfrow=c(1,2))
# Model selection
cisTopicObject <- selectModel(cisTopicObject, type="maximum")
dev.off()
cat("## Cell Embedding preparation\n##\n")
embeddings  <- modelMatSelection(
  cisTopicObject, 
  target = "cell", 
  method = "Z-score"
)

# Sort samples to match same order as in the RNA
embeddings <- t(embeddings[,colnames(object)])

# We will set it up as a DimRed object to attach to our Seurat Object
object[["cisTopic"]] <- CreateDimReducObject(embeddings = embeddings, 
                                             key = "Topic_",
                                             assay = "ATAC",
                                             misc = list("cisTopicObject"=cisTopicObject))

cat("## Saving Seurat object\n##\n")
saveRDS(object, file=SEURAT)
q()