# Obtain a set of files with signatures from a single cluster


############################## OPTIONS MENU ###################################

library("optparse")

# specify our desired options in a list
# by default OptionParser will add an help option equivalent to 
# make_option(c("-h", "--help"), action="store_true", default=FALSE, 
#               help="Show this help message and exit")

cat("\n")
option_list <- list(
  make_option("--signatures_file", type="character", 
              help="[REQUIRED] Path to the result of FindMarkerFeatures function, saved as a .csv file")
)

############################## PARSE OPTIONS ##################################


# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults
opt <- parse_args(OptionParser(option_list=option_list))

if (is.null(opt$signatures_file)) {
  write("Option --signatures_file is required\nTry --help for help", stderr())
  q()
} else {
  SIG.FILE <- opt$signatures_file
}

############################## EXECUTION #####################################

signatures <- read.csv2(SIG.FILE, header=T, row.names=1)

for (c in unique(signatures$cluster)) {
  df <- subset(signatures, cluster == c)
  filename <- sub(".csv", paste0("_cluster", c, ".tsv"), SIG.FILE)
  colnames(df) <- c("p_val", "avg_log2FC", "pct.1", "pct.2", "p_val_adj", "cluster", "geneID")
  write.table(df, file = filename, row.names = FALSE, sep = "\t")
}

