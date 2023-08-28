# conda activate /hpcnfs/data/FN/R403

.libPaths("/hpcnfs/data/FN/R403/lib/R/library")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
    cat("\nUsage: Rscript stat_MBC_myFunc.R <in_dir>\n")
    cat("\n<out_dir>       dir containing MBC.tsv\n\n")
    q()
}

IN_DIR <- args[1]

library("deMULTIplex")
library("KernSmooth")
library("ggplot2")

MBC <- read.table(file.path(IN_DIR, "MBC.txt"), sep = "\t")
bar.ref <- MBC$V2
names(bar.ref) <- MBC$V1

OUT_DIR <- file.path(IN_DIR, "out")
setwd(OUT_DIR)

bar.table <- read.table("bar.table.tsv", sep = "\t")
q <- drop(as.matrix(read.table("thresh3.myFunc.txt")))

bar.table <- bar.table[,names(bar.ref)]

# Normalize Data: Log2 Transform, mean-center
bar.table.n <- as.data.frame(log2(bar.table))
for (i in 1:ncol(bar.table)) {
    ind <- which(is.finite(bar.table.n[,i]) == FALSE)
    bar.table.n[ind,i] <- 0
    bar.table.n[,i] <- bar.table.n[,i]-mean(bar.table.n[,i])
}

n_BC <- ncol(bar.table.n) # Number of barcodes
pdf("density_BC_call.MBC.myFunc.pdf", width = min(3,n_BC)*2.4, height = 2.5*ceiling(n_BC/3))
par(mfrow=c(ceiling(n_BC/3),min(3,n_BC)))
for (i in 1:n_BC) {

    ## Step 1: GKDE with bad barcode detection, outlier trimming
    model <- tryCatch( { approxfun(bkde(bar.table.n[,i], kernel="normal")) },
          error=function(e) { print(paste0("No threshold found for ", colnames(bar.table.n)[i],"...")) } )
    if (class(model) == "character") { next }
    x <-  seq(from=quantile(bar.table.n[,i],0.001), to=quantile(bar.table.n[,i],0.999), length.out=100)
    x_quant <- seq(from=0.001, to=0.999, length.out=100)

    plot(x = x, y = model(x), main = colnames(bar.table)[i], xlab = "normalized MBC count", ylab = "density", type = "l", lwd = 2, col = "gray")

    ## Step 2: Local maxima definition
    extrema <- localMaxima(model(x))
    if (length(extrema) <= 1) {
      print(paste0("No threshold found for ", colnames(barTable.n)[i],"..."))
      next 
    }

    ## Step 3: Select maxima
    ## Assumes negative cells are largest mode
    ## Assumes positive cells are highest extreme -- favors negatives over doublets
    ## -> NEW: assumes positive cells are the second largest mode
    low.extreme <- extrema[which.max(model(x)[extrema])]
    # high.extreme <- max(extrema)
    extrema <- extrema[extrema != low.extreme]
    high.extreme <- extrema[which.max(model(x)[extrema])]
    
    abline(v = x[low.extreme])
    abline(v = x[high.extreme])
    thresh <- quantile(c(x[high.extreme], x[low.extreme]), q)
    abline(v = thresh, col = "blue", lwd = 2)
}
dev.off()


q()

