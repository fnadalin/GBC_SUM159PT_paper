#' Merge IDR sizelist with cisTopic AUC
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
  make_option("--sizelist", type="character", 
              help="[REQUIRED] Output file of cistopic_idr.R"),
  make_option("--auc", type="character",
              help="[REQUIRED] Output file of cistopic_calculate_auc_IDR.R"),
  make_option("--out", type="character",
              help="[REQUIRED] Output file name")
)

############################## PARSE OPTIONS ##################################


# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults
opt <- parse_args(OptionParser(option_list=option_list))

if (is.null(opt$sizelist)) {
  write("Option --sizelist is required\nTry --help for help", stderr())
  q()
} else {
  SIZELIST <- opt$sizelist
}

if (is.null(opt$auc)) {
  write("Option --auc is required\nTry --help for help", stderr())
  q()
} else {
  AUC <- opt$auc
}

if (is.null(opt$out)) {
  write("Option --out is required\nTry --help for help", stderr())
  q()
} else {
  OUT <- opt$out
}
################################# EXECUTION ###################################

idr.sizelist <- read.csv2(SIZELIST, header = T)
cistopic.auc <- read.csv2(AUC, header = T)

# get colnames of columns containing cisTopic1
cisTopic1.idx <- grep(pattern = "cisTopic1", colnames(cistopic.auc))
cisTopic2.idx <- grep(pattern = "cisTopic2", colnames(cistopic.auc))

# merge sizelist with AUC data by cisTopic replicate
output.1 <- merge(x=idr.sizelist, y=cistopic.auc[,c(1,cisTopic1.idx)], by.x="cisTopic1", by.y="topic", all.x=TRUE)
output.2 <- merge(x=output.1, y=cistopic.auc[,c(1,cisTopic2.idx)], by.x="cisTopic2", by.y="topic", all.x=TRUE)

write.csv2(output.2, file = OUT, row.names = F, quote = F)
q()
