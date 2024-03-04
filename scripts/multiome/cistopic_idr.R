#'
#' Calculate IDR scores between two cisTopic experiments.
#'
#' Takes as input two cisTopic objects and return for each pair of topics a bed file containing IDR values
#' 
#' @param cisTopic1 name of the first cisTopic experiment. It must link to an RDS file obtained by cistopic_training.R
#' @param cisTopic2 same as cisTopic1
#' @param idr_threshold threshold of significance for the IDR calculated for each region. It defaults to 0.05
#' @param ncpus Number of cores for the analysis
#' @param skip_bed_writing a flag to avoid rewriting bed files for each topic (time intensive procedure). Defaults to FALSE
#' 
#' @return a folder nested in the folder of cisTopic1 named IDR_files containing all bed files written \
#' and a file containing information about the number of significantly reproducible among each topic couple \
#' and the respective quantile of the pvalues (called size_list)
#' 

############################## OPTIONS MENU ###################################

library("optparse")

# specify our desired options in a list
# by default OptionParser will add an help option equivalent to 
# make_option(c("-h", "--help"), action="store_true", default=FALSE, 
#               help="Show this help message and exit")

cat("\n")
option_list <- list(
  make_option("--cisTopic1", type="character", 
              help="[REQUIRED] cisTopic1 Object path"),
  make_option("--cisTopic2", type="character", 
              help="[REQUIRED] cisTopic1 Object path"),
  make_option("--idr_threshold", type="numeric", default="0.05",
              help="[OPTIONAL] IDR threshold for lowly reproducible regions"),
  make_option("--ncpus", type="numeric", default=4,
              help="[OPTIONAL] Number of cores to use in the IDR calculation [default=%default]"),
  make_option("--skip_bed_writing", type="logical", default=FALSE,
              help="[OPTIONAL] Logical to indicate whether to write topic-specific bed files with associated region probabilities"),
  make_option("--skip_IDR_calculations", type="logical", default=FALSE,
              help="[OPTIONAL] Logical to indicate whether to skip all calculations and prepare only the sizelist")
)

############################## PARSE OPTIONS ##################################


# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults
opt <- parse_args(OptionParser(option_list=option_list))

if (is.null(opt$cisTopic1)) {
  write("Option --cisTopic1 is required\nTry --help for help", stderr())
  q()
} else {
  CISTOPIC1 <- opt$cisTopic1
}

if (is.null(opt$cisTopic2)) {
  write("Option --cisTopic2 is required\nTry --help for help", stderr())
  q()
} else {
  CISTOPIC2 <- opt$cisTopic2
}

ncpus <- opt$ncpus

IDR_THR <- opt$idr_threshold
SKIP_BED <- opt$skip_bed_writing
SKIP_CALC <- opt$skip_IDR_calculations
################################# EXECUTION ###################################

suppressPackageStartupMessages({
  library(cisTopic)
  library(GenomicRanges)
  library(parallel)
  })

cat("##-------------------------------- cisTopic analysis --------------------------------------\n##\n")
cat("## DISCLAIMER: the reproducibility score calculated within this script is deprecated.\n") 
cat("## Please run calculate_IDR_repr_score.R after this.\n##\n")
cat("## FILE READING:...")

cisTopic1 <- readRDS(CISTOPIC1)
cisTopic2 <- readRDS(CISTOPIC2)

cisTopic1 <- selectModel(cisTopic1, type="maximum")
cisTopic2 <- selectModel(cisTopic2, type="maximum")

cisTopic1.name <- sub(".Rds","", CISTOPIC1)
cisTopic2.name <- sub(".Rds","", CISTOPIC2)

cat("done!\n##\n")
if(!SKIP_CALC) {
cat("## Extracting probability matrices from cisTopics...")
object1.reg.emb.p <- modelMatSelection(cisTopic1,
                                        target = "region",
                                        method = "Probability",
                                        all.regions = TRUE)

object2.reg.emb.p <- modelMatSelection(cisTopic2,
                                       target = "region",
                                       method = "Probability",
                                       all.regions = TRUE)
cat("done!\n##\n")

if (!SKIP_BED) {
cat("## Saving .bed files for each topic with region probabilities...\n## ")
for (topic in 1:40) {
  cat(topic)
  cat(" ")
  bed.file.1 <- data.frame(chrom=seqnames(cisTopic1@region.ranges),
                           chromStart=start(cisTopic1@region.ranges)-1,
                           chromEnd=end(cisTopic1@region.ranges))
  
  bed.file.2 <- data.frame(chrom=seqnames(cisTopic2@region.ranges),
                           chromStart=start(cisTopic2@region.ranges)-1,
                           chromEnd=end(cisTopic2@region.ranges))
  
  bed.file.1 <- cbind(bed.file.1,
                     object1.reg.emb.p[topic,],
                     rep(".", length(cisTopic1@region.ranges)),
                     strand(cisTopic1@region.ranges),
                     rep(0.0, length(cisTopic1@region.ranges)),
                     rep(0.0, length(cisTopic1@region.ranges)),
                     rep(0.0, length(cisTopic1@region.ranges))
                     )
  
  bed.file.2 <- cbind(bed.file.2,
                     object2.reg.emb.p[topic,],
                     rep(".", length(cisTopic2@region.ranges)),
                     strand(cisTopic2@region.ranges),
                     rep(0.0, length(cisTopic2@region.ranges)),
                     rep(0.0, length(cisTopic2@region.ranges)),
                     rep(0.0, length(cisTopic2@region.ranges))
  )
  
  bed.file.1.name <- paste0(cisTopic1.name, "_topic_",topic,"_prob",".bed")
  bed.file.2.name <- paste0(cisTopic2.name, "_topic_",topic,"_prob", ".bed")
  
  write.table(bed.file.1, file=bed.file.1.name, row.names=F, col.names = F, quote = F, sep="\t")
  write.table(bed.file.2, file=bed.file.2.name, row.names=F, col.names = F, quote = F, sep="\t")
}
cat("\n##\n")
}
run.idr <- function(topics) {
  topic1 <- topics[1]
  topic2 <- topics[2]
  bed.file.1.name <- paste0(cisTopic1.name, "_topic_",topic1,"_prob",".bed")
  bed.file.2.name <- paste0(cisTopic2.name, "_topic_",topic2,"_prob", ".bed")
  cisTopic1.path <- dirname(cisTopic1.name)
  dir.create(file.path(cisTopic1.path, "IDR_files"), showWarnings = FALSE)
  out.file.name <- file.path(cisTopic1.path, "IDR_files", 
                             paste(basename(cisTopic1.name),
                             basename(cisTopic2.name), 
                             topic1, topic2,
                             "IDR",
                             paste0(IDR_THR,".bed"), sep="_"))
  
  command <- "idr"
  args <- c("--samples", bed.file.1.name, bed.file.2.name,
                "--input-file-type bed",
                "--rank", 4,
                "--output-file", out.file.name,
                "--peak-merge-method", "max",
                "--idr-threshold", IDR_THR,
                "--max-iter", 100,
                "--verbose")
  system2(
    command = command,
    args = args
  )
}

cat("## Running IDR calculations (this may take some time)...")
combinations.topics <- gtools::permutations(40,2,repeats.allowed = TRUE)
mcmapply(run.idr, 
         split(combinations.topics, f=1:nrow(combinations.topics)),
         mc.cores = ncpus)
cat("done!\n##\n")


}
cat("## Measuring set size of significantly reproducible peaks...")
count.lines <- function(file) {
  total_records <- as.integer(system2("wc",
                                      args = c("-l",
                                               file,
                                               " | awk '{print $1}'"),
                                      stdout = TRUE))
  return(total_records)
}


combinations.topics <- gtools::permutations(40,2,repeats.allowed = TRUE)
cisTopic1.path <- dirname(cisTopic1.name)
measure.idr.overlap <- function(topic1, topic2) {
  out.file.name <- file.path(cisTopic1.path, "IDR_files", 
                             paste(basename(cisTopic1.name),
                                   basename(cisTopic2.name), 
                                   topic1, topic2,
                                   "IDR",
                                   paste0(IDR_THR,".bed"), sep="_"))
  if (count.lines(out.file.name) == 0) {
    return(0)
  } else {
    topic.regions <- read.table(out.file.name)
    topic.regions[,9:14] <- NULL
    colnames(topic.regions) <- c("chr", "start", "end", "name", "IDR", "strand", "lIDR", "gIDR")
    topic.region.filtered <- subset(topic.regions, IDR > 0)
    return(dim(topic.region.filtered)[1])
  }
}


measure.idr.overlap.q75 <- function(topic1, topic2) {
  out.file.name <- file.path(cisTopic1.path, "IDR_files", 
                             paste(basename(cisTopic1.name),
                                   basename(cisTopic2.name), 
                                   topic1, topic2,
                                   "IDR",
                                   paste0(IDR_THR,".bed"), sep="_"))
  if (count.lines(out.file.name) == 0) {
    return(0)
  } else {
    topic.regions <- read.table(out.file.name)
    topic.regions[,9:14] <- NULL
    colnames(topic.regions) <- c("chr", "start", "end", "name", "IDR", "strand", "lIDR", "gIDR")
    topic.region.filtered <- subset(topic.regions, IDR > 0)
    return(quantile(topic.region.filtered$IDR, 0.75))
  }
}

overlap <- cbind(combinations.topics, 
                 apply(combinations.topics,1, function(topics) measure.idr.overlap(topics[1], topics[2])))

overlap <- cbind(overlap, 
                 apply(combinations.topics,1, function(topics) measure.idr.overlap.q75(topics[1], topics[2])))
# 
# 
# ## This reproducibility score does not work. To have a 
# overlap$repr.score <- rowMeans(overlap[,c(3,4)])
# 
# # before normalization, we apply a clipping to contain outliers effect
# fence <- function(vec, UB=1, LB=0) pmax( LB, pmin( vec, UB))
# overlap$repr.score <- with(overlap, fence(repr.score, UB =quantile(repr.score, 0.95)))
# overlap$repr.score <- with(overlap, repr.score/max(repr.score))

cat("done!\n##\n")
colnames(overlap) <- c("cisTopic1", "cisTopic2", "IDR5", "IDR.q75")
cat("## Writing overlap list to file...")
write.csv2(overlap, paste0(cisTopic1.name, "_", basename(cisTopic2.name), "_IDR_",IDR_THR,"_size_list.csv"), 
           row.names = FALSE, quote = FALSE)
cat("done!\n")
