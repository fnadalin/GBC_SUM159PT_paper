#'
#' Script to fuse all signatures and get their coordinates into a single csv file
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
  make_option("--signatures", type="character", 
              help="[REQUIRED] signature file names as comma-separated list"),
  make_option("--names", type="character", 
              help="[REQUIRED] Signatures names as comma-separated list"),
  make_option("--out", type="character",
              help="[OPTIONAL] Output filename")
)

############################## PARSE OPTIONS ##################################


# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults
opt <- parse_args(OptionParser(option_list=option_list))

if (is.null(opt$signatures)) {
  write("Option --signatures is required\nTry --help for help", stderr())
  q()
} else {
  SIGNATURES <- unlist(strsplit(opt$signatures,","))
}

if (is.null(opt$names)) {
  write("Option --names is required\nTry --help for help", stderr())
  q()
} else {
  NAMES <- unlist(strsplit(opt$names,","))
}

if (is.null(opt$out)) {
  write("Option --out is required\nTry --help for help", stderr())
  q()
} else {
  OUT <- opt$out
}
################################# EXECUTION ###################################

library(biomaRt)

paste5 <- function(..., sep = " ", collapse = NULL, na.rm = T) {
  if (na.rm == F)
    paste(..., sep = sep, collapse = collapse)
  else
    if (na.rm == T) {
      paste.na <- function(x, sep) {
        x <- gsub("^\\s+|\\s+$", "", x)
        ret <- paste(na.omit(x), collapse = sep)
        is.na(ret) <- ret == ""
        return(ret)
      }
      df <- data.frame(..., stringsAsFactors = F)
      ret <- apply(df, 1, FUN = function(x) paste.na(x, sep))
      
      if (is.null(collapse))
        ret
      else {
        paste.na(ret, sep = collapse)
      }
    }
}

SIG.TABLES <- lapply(SIGNATURES, read.csv2, col.names = c("gene"))
names(SIG.TABLES) <- NAMES

SIG.TABLES <- lapply(1:length(SIG.TABLES),
                     function(i){
                       SIG.TABLES[[i]]$cluster <- rep(NAMES[i], nrow(SIG.TABLES[[i]]))
                       return(SIG.TABLES[[i]])
                     })

signatures <- Reduce(x=SIG.TABLES, f=function(df1,df2) merge(df1,df2, by="gene", all=TRUE))

signatures$cluster <- apply(signatures[,grep("cluster", colnames(signatures))], 1, paste5, collapse="/")

signatures[,grep('cluster.[xy]',  colnames(signatures))] <- NULL

ensembl <- useMart("ensembl")
ensembl <- useDataset("hsapiens_gene_ensembl",mart=ensembl)
coords <- getBM(attributes=c('hgnc_symbol', 'chromosome_name', 'start_position', 'end_position'),
      filters=c('hgnc_symbol', 'chromosome_name'),
      values=list(signatures$gene, c(1:22, "X","Y")),
      mart=ensembl)

coords$coords <- with(coords, paste0("chr",chromosome_name,":",start_position,"-", end_position))

output <- merge(signatures, coords[,c("hgnc_symbol","coords")], by.x="gene", by.y="hgnc_symbol")

write.table(output, OUT, row.names = F, quote=F)
q()

