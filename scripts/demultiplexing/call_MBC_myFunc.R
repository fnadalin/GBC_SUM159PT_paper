# conda activate /hpcnfs/data/FN/R403

.libPaths("/hpcnfs/data/FN/R403/lib/R/library")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
    cat("\nUsage: Rscript call_MBC_myFunc.R <in_dir> <script_dir>\n")
    cat("\n<in_dir>       dir containing MBC.txt, barcodes.tsv.gz, lanes.csv\n")
    cat("<script_dir>   dir containing myFunc.R\n\n")
    q()
}
      
IN_DIR <- args[1]
SCRIPT_DIR <- args[2]

source(file.path(SCRIPT_DIR, "myFunc.R"))

library("ggplot2")

setwd(IN_DIR)

MBC_FILE <- "MBC.txt"
CB_FILE <- "barcodes.tsv.gz"
LANES <- "lanes.csv"

OUT_DIR <- "out"
dir.create(OUT_DIR, showWarnings = FALSE)

MBC <- read.table(MBC_FILE, sep = "\t")
bar.ref <- MBC$V2
names(bar.ref) <- MBC$V1

cell.id.vec <- drop(as.matrix(read.table(gzfile(CB_FILE))))
cell.id.vec <- gsub("-1$", "", cell.id.vec)

lanes <- read.table(LANES, header = TRUE, sep = ",")

for (i in 1:nrow(lanes)) {
    R1 <- lanes$R1[i]
    R2 <- lanes$R2[i]
    readTable <- MULTIseq.preProcess(R1 = R1, R2 = R2, cellIDs = cell.id.vec, cell=c(1,16), umi=c(17,28), tag=c(1,8))
    bar.table.l1 <- MULTIseq.align(readTable, cell.id.vec, bar.ref)
    if (i == 1) {
        bar.table <- bar.table.l1
    } else {
        bar.table <- bar.table + bar.table.l1
    }
}

colnames(bar.table)[1:length(bar.ref)] <- names(bar.ref)
outfile <- file.path(OUT_DIR, "bar.table.tsv")
write.table(bar.table, file = outfile, sep = "\t")

bar.table <- bar.table.full <- bar.table[,names(bar.ref)]
bar.tsne <- barTSNE(bar.table) 
outfile <- file.path(OUT_DIR, "bc.check.pdf")
pdf(outfile)
for (i in 3:ncol(bar.tsne)) {
    g <- ggplot(bar.tsne, aes(x = TSNE1, y = TSNE2, color = bar.tsne[,i])) +
    geom_point() +
    scale_color_gradient(low = "black", high = "red") +
    ggtitle(colnames(bar.tsne)[i]) +
    theme(legend.position = "none") 
    print(g)
}
dev.off()

bar.table_sweep.list <- list()
n <- 0
for (q in seq(0.01, 0.99, by=0.02)) {
  print(q)
  n <- n + 1
  bar.table_sweep.list[[n]] <- classifyCellsFN(bar.table, q=q)
  names(bar.table_sweep.list)[n] <- paste("q=",q,sep="")
}

#### round 1 ####

threshold.results1 <- findThresh(call.list=bar.table_sweep.list)
outfile <- file.path(OUT_DIR, "threshold.myFunc.pdf")
pdf(outfile)
ggplot(data=threshold.results1$res, aes(x=q, y=Proportion, color=Subset)) + geom_line() + theme(legend.position = "none") +
  geom_vline(xintercept=threshold.results1$extrema, lty=2) + scale_color_manual(values=c("red","black","blue"))
dev.off()

thresh <- findQ(threshold.results1$res, threshold.results1$extrema)
round1.calls <- classifyCellsFN(bar.table, q=thresh)
outfile <- file.path(OUT_DIR, "thresh.myFunc.txt")
write(thresh, file = outfile)
outfile <- file.path(OUT_DIR, "round1.calls.myFunc.tsv")
write.table(round1.calls, file = outfile, sep = "\t", quote = FALSE, col.names = FALSE)

#### round 2 ####

neg.cells <- names(round1.calls)[which(round1.calls == "Negative")]
bar.table <- bar.table[-which(rownames(bar.table) %in% neg.cells), ]

bar.table_sweep.list <- list()
n <- 0
for (q in seq(0.01, 0.99, by=0.02)) {
  print(q)
  n <- n + 1
  bar.table_sweep.list[[n]] <- classifyCellsFN(bar.table, q=q)
  names(bar.table_sweep.list)[n] <- paste("q=",q,sep="")
}

threshold.results2 <- findThresh(call.list=bar.table_sweep.list)
outfile <- file.path(OUT_DIR, "threshold2.myFunc.pdf")
pdf(outfile)
ggplot(data=threshold.results2$res, aes(x=q, y=Proportion, color=Subset)) + geom_line() + theme(legend.position = "none") +
  geom_vline(xintercept=threshold.results2$extrema, lty=2) + scale_color_manual(values=c("red","black","blue"))
dev.off()

thresh2 <- findQ(threshold.results2$res, threshold.results2$extrema)
round2.calls <- classifyCellsFN(bar.table, q=thresh2)
outfile <- file.path(OUT_DIR, "thresh2.myFunc.txt")
write(thresh2, file = outfile)
outfile <- file.path(OUT_DIR, "round2.calls.myFunc.tsv")
write.table(round2.calls, file = outfile, sep = "\t", quote = FALSE, col.names = FALSE)

#### round 3 ####

neg.cells <- c(neg.cells, names(round2.calls)[which(round2.calls == "Negative")])
bar.table <- bar.table[-which(rownames(bar.table) %in% neg.cells), ]

bar.table_sweep.list <- list()
n <- 0
for (q in seq(0.01, 0.99, by=0.02)) {
  print(q)
  n <- n + 1
  bar.table_sweep.list[[n]] <- classifyCellsFN(bar.table, q=q)
  names(bar.table_sweep.list)[n] <- paste("q=",q,sep="")
}

threshold.results3 <- findThresh(call.list=bar.table_sweep.list)
outfile <- file.path(OUT_DIR, "threshold3.myFunc.pdf")
pdf(outfile)
ggplot(data=threshold.results3$res, aes(x=q, y=Proportion, color=Subset)) + geom_line() + theme(legend.position = "none") +
  geom_vline(xintercept=threshold.results3$extrema, lty=2) + scale_color_manual(values=c("red","black","blue"))
dev.off()

thresh3 <- findQ(threshold.results3$res, threshold.results3$extrema)
round3.calls <- classifyCellsFN(bar.table, q=thresh3)
outfile <- file.path(OUT_DIR, "thresh3.myFunc.txt")
write(thresh3, file = outfile)
outfile <- file.path(OUT_DIR, "round3.calls.myFunc.tsv")
write.table(round3.calls, file = outfile, sep = "\t", quote = FALSE, col.names = FALSE)

#### end rounds ####

final.calls <- c(round3.calls, rep("Negative",length(neg.cells)))
names(final.calls) <- c(names(round3.calls),neg.cells)
outfile <- file.path(OUT_DIR, "final.calls.myFunc.txt")
write.table(final.calls, file = outfile, sep = "\t", col.names = FALSE, quote = FALSE)

#### reclassify ####

reclass.cells <- findReclassCellsFN(bar.table.full, names(final.calls)[which(final.calls=="Negative")])
reclass.res <- rescueCells(bar.table.full, final.calls, reclass.cells)

outfile <- file.path(OUT_DIR, "rescued_cells.myFunc.pdf")
pdf(outfile)
ggplot(reclass.res[-1, ], aes(x=ClassStability, y=MatchRate_mean)) + 
    geom_point() + xlim(c(nrow(reclass.res)-1,1)) + 
    ylim(c(0,1.05)) +
    geom_errorbar(aes(ymin=MatchRate_mean-MatchRate_sd, ymax=MatchRate_mean+MatchRate_sd), width=.1) +
    geom_hline(yintercept = reclass.res$MatchRate_mean[1], color="red") +
    geom_hline(yintercept = reclass.res$MatchRate_mean[1]+3*reclass.res$MatchRate_sd[1], color="red",lty=2) +
    geom_hline(yintercept = reclass.res$MatchRate_mean[1]-3*reclass.res$MatchRate_sd[1], color="red",lty=2)
dev.off()

final.calls.rescued <- final.calls
rescue.ind <- which(reclass.cells$ClassStability >= 20) ## Note: Value will be dataset-specific
final.calls.rescued[rownames(reclass.cells)[rescue.ind]] <- reclass.cells$Reclassification[rescue.ind]
outfile <- file.path(OUT_DIR, "final.calls.rescued.myFunc.txt")
write.table(final.calls.rescued, file = outfile, sep = "\t", col.names = FALSE, quote = FALSE)


q()

