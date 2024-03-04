#' clean sizelist from topics that correlate with coverage over a specific threshold
#'
#' @return a sizelist-like table with only topics that do not correlate with coverage
#' 


############################## OPTIONS MENU ###################################

library("optparse")

# specify our desired options in a list
# by default OptionParser will add an help option equivalent to 
# make_option(c("-h", "--help"), action="store_true", default=FALSE, 
#               help="Show this help message and exit")

cat("\n")
option_list <- list(
  make_option("--idr_file", type="character", 
              help="[REQUIRED] Path to the result of cistopic_idr.R"),
  make_option("--cisTopic1", type="character", 
              help="[REQUIRED] Topics numbers from replicate 1 that do correlate with coverage"),
  make_option("--cisTopic2", type="character", 
              help="[REQUIRED] Topics numbers from replicate 1 that do correlate with coverage"),
  make_option("--out", type="character",
              help="[REQUIRED] Output file name")
)

############################## PARSE OPTIONS ##################################


# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults
opt <- parse_args(OptionParser(option_list=option_list))

if (is.null(opt$idr_file)) {
  write("Option --idr_file is required\nTry --help for help", stderr())
  q()
} else {
  IDR.FILE <- opt$idr_file
}

if (is.null(opt$cisTopic1)) {
  write("Option --cisTopic1 is required\nTry --help for help", stderr())
  q()
} else {
  CISTOPIC1 <- as.numeric(unlist(strsplit(opt$cisTopic1, ",")))
}

if (is.null(opt$cisTopic2)) {
  write("Option --cisTopic2 is required\nTry --help for help", stderr())
  q()
} else {
  CISTOPIC2 <- as.numeric(unlist(strsplit(opt$cisTopic2, ",")))
}

if (is.null(opt$out)) {
  write("Option --out is required\nTry --help for help", stderr())
  q()
} else {
  OUT <- opt$out
}
############################## EXECUTION #####################################

IDR <- read.csv2(IDR.FILE, header=T)

IDR.clean <- subset(IDR, !(cisTopic1 %in% CISTOPIC1) & !(cisTopic2 %in% CISTOPIC2))

write.csv2(IDR.clean, file = OUT, row.names = F, quote = F)
q()
