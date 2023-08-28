
# label the cell as "true" if ALL GBCs are expressed

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
    cat("\nUsage: <gbc_list> <metadata.tsv> <label>\n")
    cat("\n<gbc_list>      where GBC IDs are stored\n")
    cat("<metadata.tsv>  where cell metadata are stored; must contain expr.GBC.list\n")
    cat("<label>         name of the metadata field\n\n")
    q()
}


GBC_LIST <- args[1]
METADATA <- args[2]
LABEL <- args[3]

gbc_list <- read.table(GBC_LIST)[,1]
df <- read.table(METADATA, sep = "\t", header = TRUE)

# all GBC expressed by the cell 
is_found <- unlist(lapply(df$expr.GBC.list, function(x) {
                                            v <- unlist(strsplit(x, split=","))
                                            sum(v %in% gbc_list) == length(v)} ))

meta <- rep(0, length(is_found))
meta[is_found] <- rep(1, sum(is_found))

names <- colnames(df)
df <- cbind(df, meta)
colnames(df) <- c(names, LABEL)

write.table(df, file = METADATA, sep = "\t", quote = FALSE)


q()

