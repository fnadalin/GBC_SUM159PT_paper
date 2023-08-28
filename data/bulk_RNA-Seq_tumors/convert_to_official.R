IN <- "SUM159_in_vivo_Vs_in_vitro"
OUT <- "genes_official.tsv"

library("org.Hs.eg.db")
library("limma")

df <- read.table(IN, sep = "\t", header = TRUE)

df[,10] <- alias2SymbolTable(df[,10], species = "Hs") # convert to official gene names

write.table(df, file = OUT, row.names = FALSE, sep = "\t", quote = FALSE)

q()
 

