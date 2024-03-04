
options(warn=-1)

# Takes as input a pair of cisTopic objects and return a file with a column named topic and other n columns named "Cluster.\d.cisTopic{1,2}"

############################## OPTIONS MENU ###################################

library("optparse")

# specify our desired options in a list
# by default OptionParser will add an help option equivalent to 
# make_option(c("-h", "--help"), action="store_true", default=FALSE, 
#               help="Show this help message and exit")

cat("\n")
option_list <- list(
  make_option("--cisTopic1", type="character", 
              help="[REQUIRED] cisTopic Object 1 path (must have same topic search space size as cisTopic2)"),
  make_option("--cisTopic2", type="character", 
              help="[REQUIRED] cisTopic Object 2 path (must have same topic search space size as cisTopic1)"),
  make_option("--metadata1", type="character",
              help="[REQUIRED] File containing Seurat metadata information for object cisTopic1"),
  make_option("--metadata2", type="character",
              help="[REQUIRED] File containing Seurat metadata information for object cisTopic2"),
  make_option("--clusters1", type="character",
              help="[REQUIRED] Column name of Seurat clustering for object cisTopic1"),
  make_option("--clusters2", type="character",
              help="[REQUIRED] Column name of Seurat clustering for object cisTopic2"),
  make_option("--out", type="character",
              help="[REQUIRED] Output file name")
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

if (is.null(opt$metadata1)) {
  write("Option --metadata1 is required\nTry --help for help", stderr())
  q()
} else {
  METADATA1 <- opt$metadata1
}

if (is.null(opt$metadata2)) {
  write("Option --metadata2 is required\nTry --help for help", stderr())
  q()
} else {
  METADATA2 <- opt$metadata2
}

if (is.null(opt$clusters1)) {
  write("Option --clusters1 is required\nTry --help for help", stderr())
  q()
} else {
  CLUSTERS1 <- opt$clusters1
}

if (is.null(opt$clusters2)) {
  write("Option --clusters2 is required\nTry --help for help", stderr())
  q()
} else {
  CLUSTERS2 <- opt$clusters2
}

if (is.null(opt$out)) {
  write("Option --out is required\nTry --help for help", stderr())
  q()
} else {
  OUT <- opt$out
}


################################# EXECUTION ###################################
suppressPackageStartupMessages({
  library(cisTopic)
  library(pROC)
  })

topic.auc <- function(topic,labels) {
  return(sapply(
    X=levels(labels),
    FUN=function(x, topic, labels) auc(labels == x, topic, direction='<'),
    topic=topic,
    labels=labels
  ))
}

calculate.auc <- function(df, labels) {
  return(apply(
    X=df,
    MARGIN=1,
    FUN = function(x, labels) topic.auc(x, labels),
    labels=labels
  ))
}

# reading input files

cisTopic1 <- readRDS(CISTOPIC1)
cisTopic2 <- readRDS(CISTOPIC2)
meta1 <- read.csv2(METADATA1, row.names = 1)
meta2 <- read.csv2(METADATA2, row.names = 1)

cisTopic1 <- selectModel(cisTopic1, type = "maximum")
cisTopic2 <- selectModel(cisTopic2, type = "maximum")


# calculate cell embeddings
cell.embeddings.1 <- modelMatSelection(cisTopic1,
                                       target = "cell",
                                       method = "Z-score")
cell.embeddings.2 <- modelMatSelection(cisTopic2,
                                       target = "cell",
                                       method = "Z-score")

# calculating AUC datasets
cisTopic1.AUC <- as.data.frame(t(calculate.auc(cell.embeddings.1, as.factor(meta1[[CLUSTERS1]]))))
cat(paste(dim(cisTopic1.AUC), collapse = " by "))
colnames(cisTopic1.AUC) <- paste("Cluster", seq(0,ncol(cisTopic1.AUC)-1), "cisTopic1", sep=".")
cisTopic1.AUC$topic <- seq(1,nrow(cisTopic1.AUC))

cisTopic2.AUC <- as.data.frame(t(calculate.auc(cell.embeddings.2, as.factor(meta2[[CLUSTERS2]]))))
colnames(cisTopic2.AUC) <- paste("Cluster", seq(0,ncol(cisTopic2.AUC)-1), "cisTopic2", sep=".")
cisTopic2.AUC$topic <- seq(1,nrow(cisTopic2.AUC))

result <- merge(cisTopic1.AUC, cisTopic2.AUC, by="topic")

write.csv2(result, OUT, row.names = F, quote = F)

q()