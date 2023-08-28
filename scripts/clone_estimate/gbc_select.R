
MIN_CPM <- 1
MIN_LOG2FC <- -10

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
    cat("\nUsage: <cpm> <outdir>\n")
    cat("\n<cpm>     file containing cpm\n")
    cat("<outdir>  where log2fc, ma stats are found & where to save plots and stats\n\n")
    q()
}

CPM <- args[1]
OUT <- args[2]

LOG2FC <- file.path(OUT, "log2fc.tsv")
MA <- file.path(OUT, "ma.tsv")

cpm <- read.table(CPM, sep = "\t", header = TRUE)
log2fc <- read.table(LOG2FC, sep = "\t", header = TRUE)
ma <- read.table(MA, sep = "\t", header = TRUE)

gbc <- rownames(cpm)

library("ggplot2")

samples <- colnames(log2fc)
gbc_count <- rep(0,nrow(log2fc))

dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

# Compute stats on all GBCs

df <- data.frame(sample = NULL, log2FC = NULL, cpm = NULL)

stat <- NULL
for (s in samples) { 

    df_tmp <- data.frame(log2fc = log2fc[[s]], cpm = cpm[[s]])
    df_tmp <- cbind(sample = rep(s, nrow(df_tmp)), df_tmp)
    df <- rbind(df, df_tmp)
    
    idx <- which(df_tmp$log2fc > MIN_LOG2FC & df_tmp$cpm > MIN_CPM)
    stat <- c(stat, length(idx))
    gbc_count[idx] <- gbc_count[idx] + 1

    write.table(gbc[idx], paste0(OUT, "/", s,"_sel.txt"), quote = FALSE, col.names = FALSE, row.names = FALSE)
}

names(stat) <- samples
write.table(stat, file = paste0(OUT, "/stat.tsv"), sep = "\t", col.names = FALSE, quote = FALSE)

names(gbc_count) <- gbc
write.table(gbc_count, file = paste0(OUT, "/count_sel.tsv"), sep = "\t", col.names = FALSE, quote = FALSE)

gbc_count_stat <- unlist(lapply(0:max(gbc_count), function(x) sum(gbc_count >= x)))
names(gbc_count_stat) <- 0:max(gbc_count)
write.table(gbc_count_stat, file = paste0(OUT, "/count_sel_stat.tsv"), sep = "\t", col.names = FALSE, quote = FALSE)

# plot across samples and across counts

pdf(paste0(OUT, "/all.pdf"), width = 8, height = 5)
g <- ggplot(data = df, aes(x = log2fc, y = cpm)) + geom_point() 
g <- g + facet_wrap( ~ sample, nrow = 2)
print(g)
dev.off()

pdf(paste0(OUT, "/all_log.pdf"), width = 8, height = 5)
g <- ggplot(data = df, aes(x = log2fc, y = cpm)) + geom_point() + scale_y_continuous(trans='log10')
g <- g + geom_vline(xintercept=MIN_LOG2FC, linetype="dashed", color = "red")
g <- g + geom_hline(yintercept=MIN_CPM, linetype="dashed", color = "red")
g <- g + facet_wrap( ~ sample, nrow = 2)
print(g)
dev.off()

pdf(paste0(OUT, "/ma.pdf"), width = 5, height = 5)
g <- ggplot(data = ma, aes(x = A, y = M)) + geom_point(aes(color = gbc_count))
print(g)
dev.off()

df <- cbind(ma, count = gbc_count)
pdf(paste0(OUT, "/ma_by_counts.pdf"), width = 7, height = 5)
g <- ggplot(data = df, aes(x = A, y = M)) + geom_point() + facet_wrap( ~ count, nrow = 2)
print(g)
dev.off()

# Compute stats on selected GBC and link to the number of evidences (counts)

cum_count_stat <- as.data.frame(matrix(NA, nrow = length(samples), ncol = length(samples)))
colnames(cum_count_stat) <- samples
rownames(cum_count_stat) <- 1:length(samples)

df <- data.frame(sample = NULL, log2fc = NULL, cpm = NULL, count = NULL, avgrank = NULL, norm.coverage = NULL)
df1 <- data.frame(sample = NULL, count = NULL, cum.frac = NULL)

for (s in samples) { 

    log2fcRank <- match(1:nrow(log2fc),order(log2fc[[s]]))
    cpmRank <- match(1:nrow(cpm),order(cpm[[s]]))
    avgrank <- apply(data.frame(log2fcRank,cpmRank), 1, mean)
    
    # subset on selected GBCs
    gbc_sel <- read.table(paste0(OUT, "/", s,"_sel.txt"), sep = "\t", header = TRUE)[,1]
    df_tmp <- data.frame(log2fc = log2fc[[s]], cpm = cpm[[s]], count = gbc_count, avgrank = avgrank)
    df_tmp <- df_tmp[match(gbc_sel,gbc),]
    
    nclones <- length(gbc_sel)
    avg_cpm_per_clone <- sum(df_tmp$cpm)/nclones
    norm.coverage <- df_tmp$cpm/avg_cpm_per_clone
    df_tmp <- cbind(sample = rep(s, length(gbc_sel)), df_tmp, norm.coverage)
    df <- rbind(df, df_tmp)
    
    cpm_by_count <- unlist(lapply(1:length(samples), function(x) sum(df_tmp$cpm[df_tmp$count >= x])))
    cpm_by_count <- cpm_by_count / sum(df_tmp$cpm)
    exp_cpm_by_count <- unlist(lapply(1:length(samples), function(x) avg_cpm_per_clone*sum(df_tmp$count >= x)))
    exp_cpm_by_count <- exp_cpm_by_count / sum(df_tmp$cpm)
    df1 <- rbind(df1, data.frame(sample = rep(s,length(samples)), count = 1:length(samples), cum.frac = cpm_by_count, exp.cum.frac = exp_cpm_by_count))
    
    cum_count_stat[[s]] <- cpm_by_count
}

write.table(df, file = paste0(OUT, "/sel.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(df1, file = paste0(OUT, "/sel_cpmByCount.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(cum_count_stat, file = paste0(OUT, "/cpm_count_stat.tsv"), sep = "\t", quote = FALSE)

# plot stats on selected GBCs across samples

pdf(paste0(OUT, "/sel.pdf"), width=8, height = 5)
g <- ggplot(data = df, aes(x = log2fc, y = cpm)) + geom_point(aes(color = count)) + scale_y_continuous(trans='log10') 
g <- g + facet_wrap( ~ sample, nrow = 2)
print(g)
dev.off()

df$count <- as.factor(df$count)

pdf(paste0(OUT, "/sel_log2fc.pdf"), width=5, height = 5)
g <- ggplot(data = df, aes(x = count, y = log2fc)) + geom_boxplot() 
g <- g + facet_wrap( ~ sample, nrow = 2)
print(g)
dev.off()

pdf(paste0(OUT, "/sel_cpm.pdf"), width=5, height = 5)
g <- ggplot(data = df, aes(x = count, y = cpm)) + geom_boxplot() + scale_y_continuous(trans='log10') 
g <- g + facet_wrap( ~ sample, nrow = 2)
print(g)
dev.off()

pdf(paste0(OUT, "/sel_avgrank.pdf"), width=5, height = 5)
g <- ggplot(data = df, aes(x = count, y = avgrank)) + geom_boxplot()
g <- g + facet_wrap( ~ sample, nrow = 2) 
print(g)
dev.off()

df1$count <- as.factor(df1$count)

pdf(paste0(OUT, "/sel_cpmByCount.pdf"), width=5, height = 5)
g <- ggplot(data = df1, aes(x = count, y = cum.frac)) + geom_point() + ylim(0,1)
g <- g + facet_wrap( ~ sample, nrow = 2)
print(g)
dev.off()

q()


