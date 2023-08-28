
# extract the relationship between GBC abundance (in bulk) and cell count (that express those GBCs)

STEP <- 20

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) {
    cat("\nUsage: <cpm_bulk.tsv> <CB_metadata_sc.tsv> <sample_IDs_bulk> <sample_IDs_sc> <outdir>\n\n")
    q()
}

CPM_BULK <- args[1]
CB_META <- args[2]
SAMPLES_BULK <- args[3]
SAMPLES_SC <- args[4]
OUTDIR <- args[5]

library("ggplot2")
library("scales") # for alpha

cpm_bulk <- read.table(CPM_BULK, sep = "\t", header = TRUE)
cb_meta <- read.table(CB_META, sep = "\t", header = TRUE)

samples_bulk <- unlist(strsplit(SAMPLES_BULK, split=","))
samples_sc <- unlist(strsplit(SAMPLES_SC, split=","))

df_cpm_cellcount <- data.frame(cutoff = NULL, frac_cells = NULL, sc_sample = NULL, bulk_sample = NULL)
for (s1 in samples_sc) {

    gbc_list <- cb_meta$expr.GBC.list[cb_meta$sample.name == s1]
    
    for (s2 in samples_bulk) {
    
        cpm <- cpm_bulk[[s2]]
        names(cpm) <- rownames(cpm_bulk)
        cpm <- cpm[order(cpm, decreasing = TRUE)]
        
        # count the number of cells where each GBC is found (co-infections are counted n times, where n is the number of GBCs)
        gbc_list_split <- unlist(lapply(gbc_list, function(x) strsplit(x, split=",")))
        gbc_uniq <- unique(gbc_list_split)
        gbc_split_count <- unlist(lapply(gbc_uniq, function(x) sum(gbc_list_split %in% gbc_uniq)))
        gbc_split_count <- c(gbc_split_count, rep(0, sum(!(names(cpm) %in% gbc_uniq))))
        names(gbc_split_count) <- c(gbc_uniq, names(cpm)[!(names(cpm) %in% gbc_uniq)])

        cell_count <- gbc_split_count[match(names(cpm), names(gbc_split_count))]
        
        n <- length(cpm)
        df <- data.frame(GBC = names(cpm), cpm = cpm, ncells = cell_count, sc_sample = rep(s1,n), bulk_sample = rep(s2,n))
  
        # compute the fraction of GBCs that are found in at least one cell, for some CPM cut-off
        # MAX <- ceiling(max(df$cpm)/20)
        # cutoffs <- (1:MAX)*STEP
        cutoffs <- cpm[1:floor(n/STEP)*STEP]
        frac_cells_with_gbc <- unlist(lapply(cutoffs, function(x) {
                        # a <- df$cpm >= x-STEP & df$cpm < x
                        a <- df$cpm > x
                        b <- df$ncells > 0
                        sum(a&b) / sum(a)
                        }
                    ))
        
        m <- length(cutoffs)
        df1 <- data.frame(cutoff = cutoffs, frac_cells = frac_cells_with_gbc, sc_sample = rep(s1,m), bulk_sample = rep(s2,m))
        
        df_cpm_cellcount <- rbind(df_cpm_cellcount, df1)
        
    }
}

TABLE <- file.path(OUTDIR, "cpm_cellcount.tsv")
write.table(df_cpm_cellcount, file = TABLE, sep = "\t", quote = FALSE, row.names = FALSE)

df_cpm_cellcount$bulk_sample <- factor(df_cpm_cellcount$bulk_sample, levels = samples_bulk)
df_cpm_cellcount$sc_sample <- factor(df_cpm_cellcount$sc_sample, levels = samples_sc)

PLOT <- file.path(OUTDIR, "cpm_cellcount.pdf")
g <- ggplot(data = df_cpm_cellcount, aes(x = cutoff, y = frac_cells)) + geom_point(colour = "blue", alpha = 0.2) + ylim(0,1)
g <- g + facet_grid(sc_sample ~ bulk_sample)
g <- g + xlab("CPM cutoff") + ylab("fraction of GBC in cells")
pdf(PLOT, width = 1.5*length(samples_bulk), height = 1.5*length(samples_sc))
g
dev.off()

q()

