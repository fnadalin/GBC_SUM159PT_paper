#' adds marker distances to bed file with all 
#' 
#' 
#' 

############################## OPTIONS MENU ###################################

library("optparse")

# specify our desired options in a list
# by default OptionParser will add an help option equivalent to 
# make_option(c("-h", "--help"), action="store_true", default=FALSE, 
#               help="Show this help message and exit")

cat("\n")
option_list <- list(
  make_option("--idr", type="character", 
              help="[REQUIRED] BED file containing regions that passed IDR filtering"),
  make_option("--signatures", type="character",
              help="[REQUIRED] Signature file obtained by calculate_signatures.R"),
  make_option("--out", type="character",
              help="[REQUIRED] Output file name")
)

############################## PARSE OPTIONS ##################################


# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults
opt <- parse_args(OptionParser(option_list=option_list))

if (is.null(opt$idr)) {
  write("Option --idr is required\nTry --help for help", stderr())
  q()
} else {
  IDR <- opt$idr
}

if (is.null(opt$signatures)) {
  write("Option -$signatures is required\nTry --help for help", stderr())
  q()
} else {
  SIGNATURES <- opt$signatures
}

if (is.null(opt$out)) {
  write("Option --out is required\nTry --help for help", stderr())
  q()
} else {
  OUT <- opt$out
}
################################# EXECUTION ###################################

suppressPackageStartupMessages({
  library(Signac)
  library(GenomicRanges)})


get.closest.marker.features <- function(regions, signatures, gapsize=2.5e+5, ignore.strand=TRUE) {
  # for each gene, measure overlap with an arbitrary gap distance
  # return gene names and whether they fall within one of the two clusters
  regions <- makeGRangesFromDataFrame(regions)
  markers.coords <- StringToGRanges(signatures$coords, sep=c(":","-"))
  elementMetadata(markers.coords)$gene_name <- signatures$gene
  regions.markers.overlap <- findOverlaps(query=makeGRangesFromDataFrame(regions), 
                                          subject=markers.coords,
                                          maxgap=gapsize)
  
  if (length(regions.markers.overlap) == 0) {
    empty.df <- data.frame(list(NA,NA,NA,NA,NA))
    return(empty.df)
  }
  
  hits.regions <- regions[queryHits(regions.markers.overlap)]
  hits.regions$gene_name <- markers.coords[subjectHits(regions.markers.overlap)]$gene_name
  
  hits.regions.df <- cbind(as.data.frame(hits.regions), marker=hits.regions$gene_name)
  hits.regions.df$chr <- hits.regions.df$seqnames
  hits.regions.df$seqnames <- NULL
  hits.regions.df[,c("width", "strand")] <- NULL
  results <- aggregate(marker~chr + start + end, data=hits.regions.df, paste, collapse=",")
  # now we need to find the cluster combination for each region associated with it
  
  results$cluster <- sapply(results$marker, 
                            function(x) paste0(subset(signatures, gene %in% strsplit(x, split=",")[[1]])$cluster, collapse=","))
  
  return(results)
}


signatures <- read.table(SIGNATURES, header=T)
IDR.bed <- read.table(IDR)

IDR.bed[,c(7,9,10,12,13)] <- NULL
colnames(IDR.bed) <- c("chr", "start", "end", "name", "IDR", "strand", "gIDR", "1C_probability", "1E_probability")

marker.features.500k <- get.closest.marker.features(regions = IDR.bed,
                                                    signatures = signatures,
                                                    gapsize = 5e+5)
colnames(marker.features.500k) <- c("chr","start","end", "marker_500kb", "cluster_500kb")

marker.features.50k <- get.closest.marker.features(regions = IDR.bed,
                                                   signatures = signatures,
                                                   gapsize = 5e+4)
colnames(marker.features.50k) <- c("chr","start","end", "marker_50kb", "cluster_50kb")

marker.features.5k <- get.closest.marker.features(regions = IDR.bed,
                                                  signatures = signatures,
                                                  gapsize = 5e+3)
colnames(marker.features.5k) <- c("chr","start","end", "marker_5kb", "cluster_5kb")

merged.regions.marker.df <- merge(x=IDR.bed,
                                  y=marker.features.500k,
                                  all.x=TRUE)
merged.regions.marker.df <- merge(x=merged.regions.marker.df,
                                  y=marker.features.50k,
                                  all.x=TRUE)
merged.regions.marker.df <- merge(x=merged.regions.marker.df,
                                  y=marker.features.5k,
                                  all.x=TRUE)
write.csv2(subset(merged.regions.marker.df, IDR > 540), OUT,
           row.names = F, quote = F)

