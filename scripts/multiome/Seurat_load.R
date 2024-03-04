
# Take a CellRanger output matrix folder as input and generate a Seurat
# output as object

MIN_GENES <- 0
MIN_CELLS <- 0
VARS_TO_REGRESS <- NULL
set.seed(1234)
############################# OPTIONS MENU ####################################

library("optparse")

# specify our desired options in a list
# by default OptionParser will add an help option equivalent to 
# make_option(c("-h", "--help"), action="store_true", default=FALSE, 
#               help="Show this help message and exit")
cat("\n")
option_list <- list(
    make_option("--matrix_dir", type = "character",
        help="[REQUIRED] path to the input matrix directory"),
    make_option("--GBC", type="character", 
                help="[REQUIRED] path to the input genomic barcodes (GBC) file"),
    make_option("--out_RDS", type = "character",
        help="[REQUIRED] output seurat object containing the cell barcodes from the matrices"),
    make_option("--frag_path", type="character", 
                help="[REQUIRED] fragment file path"),
    make_option("--min_genes", type = "integer", default = MIN_GENES,
        help="[OPTIONAL] minimum number of detected genes to keep a cell [default=%default]"),
    make_option("--min_cells", type = "integer", default = MIN_CELLS,
        help="[OPTIONAL] minimum number of cells where a gene is detected to keep it [default=%default]")
)

############################# PARSE OPTIONS ###################################

# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults
opt <- parse_args(OptionParser(option_list=option_list))

if (is.null(opt$matrix_dir)) {
  write("Option --matrix_dir is required\nTry --help for help", stderr()) 
  q()
} else {
  MATRIX_DIR <- opt$matrix_dir
}

if (is.null(opt$out_RDS)) {
  write("Option --out_RDS is required\nTry --help for help", stderr()) 
  q()
} else {
  OUT_FILE <- opt$out_RDS
}

if (is.null(opt$GBC)) {
  write("Option --GBC is required\nTry --help for help", stderr())
  q()
} else {
  GBC <- opt$GBC
}

if (is.null(opt$frag_path)) {
  write("Option --frag_path is required\nTry --help for help", stderr())
  q()
} else {
  FRAGPATH <- opt$frag_path
}

if (!is.null(opt$min_genes))
  MIN_GENES <- as.numeric(opt$min_genes)

if (!is.null(opt$min_cells))
  MIN_CELLS <- as.numeric(opt$min_cells)

if (!is.null(opt$vars_to_regress))
  VARS_TO_REGRESS <- unlist(strsplit(opt$vars_to_regress, sep = ","))

############################ EXECUTION ########################################

# Read10X, create seurat object, load GBC dataset and filter cells accordingly
# run QC over cells for all parameters and return plots and datasets (for later
# filtering)

# import relevant libraries
library(Seurat)
library(Signac)
library(ggplot2)
library(EnsDb.Hsapiens.v86)

# Read the CellRanger output with Seurat function Read10X
data <- Read10X(data.dir=MATRIX_DIR)

# This function will generate a list of matrices for multimodal data 
# (Gene Expression and Peaks)

object <- CreateSeuratObject(
  count=data$`Gene Expression`, 
  project="P0_Multiome", 
  assay="RNA")

object[["ATAC"]] <- CreateChromatinAssay(
  count=data$Peaks,
  sep = c(":", "-"),
  genome = "hg38",
 fragments = FRAGPATH
)

annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotations) <- 'UCSC' # convert chromosome names

Annotation(object[["ATAC"]]) <- annotations

# lineage information barcode is processed into a dataframe and used both as 
# metadata and to filter cells in the Seurat experiment
clonal.meta <- read.csv2(GBC, sep="\t")

# we add a columns of barcode names because in the file are fused with the 
# experiment name

barcode.names <- gsub("^[A-Za-z0-9\\_]*-", "", rownames(clonal.meta))
clonal.meta$barcode.names <- barcode.names
rownames(clonal.meta) <- barcode.names
# we keep only cells that have a GBC barcode

object <- subset(x = object, 
                 cells = clonal.meta$barcode.names)

object <- AddMetaData(object, clonal.meta)
# We now return information about number of features, TSS Enrichment and 
# nucleosome signal

DefaultAssay(object) <- "ATAC"
object <- TSSEnrichment(object)
object <- NucleosomeSignal(object)
object$blacklist_fraction <- FractionCountsInRegion(
  object = object,
  assay = 'ATAC',
  regions = blacklist_hg38
)

cat("\nCalculating FRiP...\n")
# calculate FRiP
counts.df <- CountFragments(fragments = FRAGPATH)
rownames(counts.df) <- counts.df$CB
object <- AddMetaData(object, counts.df[Cells(object),])

object <- FRiP(object, assay = "ATAC", total.fragments = "frequency_count")
cat("\nFRiP calculated!\n")
# removing blacklist regions from the analysis

peaks <- GetAssayData(object, assay = "ATAC")
peaks <- peaks[is.na(findOverlaps(object[["ATAC"]], blacklist_hg38_unified, select = "first")),]

object[["ATAC"]] <- subset(x = object[["ATAC"]],
                 features = rownames(peaks))

# plotting all features is a single image to be saved 
VlnPlot(object=object, 
        features = c("nCount_RNA", "nCount_ATAC", "TSS.enrichment",
                     "nucleosome_signal", "blacklist_fraction", "FRiP"),
        pt.size = 0.1,
        ncol = 6
        )

plot.name <- paste("QC_metrics", 
                   sub(".RDS", ".png", basename(OUT_FILE), ignore.case = T), 
                   sep = "_")

plot.path <- paste(dirname(OUT_FILE),"plots",plot.name, sep="/")

dir.create(paste(dirname(OUT_FILE), "plots", sep = "/"), recursive = T)

ggsave(plot.path, height = 7, width = 18)

# Now we save the obtained Seurat object we can save it and close the session

saveRDS(object, file = OUT_FILE)

sessionInfo()
q()
