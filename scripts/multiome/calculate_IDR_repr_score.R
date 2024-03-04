#'
#' Script to fix \code{calculate_idr.R} reproducibility score.
#'
#' It takes as input the result file of cistopic_idr.R (size_list) 
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
              help="[REQUIRED] Output file of cistopic_idr.R")
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


################################# EXECUTION ###################################

sizelist <- read.csv2(SIZELIST, header = T)

min.max <- function(x) {
  (x-min(x))/(max(x)-min(x))
}

sizelist$repr.score <- rowMeans(apply(X = sizelist[,c(3,4)], MARGIN = 2, FUN = min.max))

write.csv2(sizelist, file = SIZELIST)
q()