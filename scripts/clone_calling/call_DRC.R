args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
    cat("\nUsage: <file_list> <out>\n")
    cat("\n<file_list>  comma-separated list of files containing clones in different experiments\n")
    cat("<out>        txt file containing the DRCs (intersection of the above)\n")
    q()
}

FILE_LIST <- args[1]
OUT <- args[2]

file_list <- unlist(strsplit(FILE_LIST, split=","))

if (length(file_list) == 0) {
    write("No list found\n", stderr()) 
    q()
}

DRC <- read.table(file_list[1])[,1]
if (length(file_list) > 1) {
    for (i in 2:length(file_list)) {
        x <- read.table(file_list[i])[,1]
        DRC <- DRC[DRC %in% x]
    }
}

write.table(DRC, file = OUT, quote = FALSE, row.names = FALSE, col.names = FALSE)

q()

