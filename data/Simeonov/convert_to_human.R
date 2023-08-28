IN <- "mmc6.tsv"
OUT <- "hgc6.tsv"

# library("org.Hs.eg.db")
library("org.Mm.eg.db")
library("limma")

df <- read.table(IN, sep = "\t", header = TRUE)

df[,1] <- alias2SymbolTable(df[,1], species = "Mm") # convert to official gene names
df[,1] <- toupper(df[,1])

write.table(df, file = OUT, row.names = FALSE, sep = "\t", quote = FALSE)

q()
 

