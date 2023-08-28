library("deMULTIplex")
require("KernSmooth")

classifyCellsFN <- function(barTable, q) {

  ## Normalize Data: Log2 Transform, mean-center
  barTable.n <- as.data.frame(log2(barTable))
  for (i in 1:ncol(barTable)) {
    ind <- which(is.finite(barTable.n[,i]) == FALSE)
    barTable.n[ind,i] <- 0
    barTable.n[,i] <- barTable.n[,i]-mean(barTable.n[,i])
  }

  ## Pre-allocate memory to save processing time/memory
  n_BC <- ncol(barTable.n) # Number of barcodes
  n_cells <- nrow(barTable.n) # Numer of cells
  bc_calls <- rep("Negative",n_cells) # Barcode classification for each cell
  names(bc_calls) <- rownames(barTable.n)

  ## Generate thresholds for each barcode:
  ## 1. Gaussian KDE with bad barcode detection, outlier trimming
  ## 2. Define local maxima for GKDE
  ## 3. Split maxima into low and high subsets, adjust low if necessary
  ## 4. Threshold and classify cells according to user-specified inter-maxima quantile

  for (i in 1:n_BC) {
    ## Step 1: GKDE
    model <- tryCatch( { approxfun(bkde(barTable.n[,i], kernel="normal")) },
                       error=function(e) { print(paste0("No threshold found for ", colnames(barTable.n)[i],"...")) } )
    if (class(model) == "character") { next }
    x <-  seq(from=quantile(barTable.n[,i],0.001), to=quantile(barTable.n[,i],0.999), length.out=100)

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
    if (low.extreme == high.extreme) {print(paste0("No threshold found for ", colnames(barTable.n)[i],"...")) }

    ## Step 4: Threshold and classify cells
    thresh <- quantile(c(x[high.extreme], x[low.extreme]), q)
    cell_i <- which(barTable.n[,i] >= thresh)
    n <- length(cell_i)
    if (n == 0) { next } # Skip to next barcode if no cells classified
    bc_calls[cell_i] <- sapply(bc_calls[cell_i],
                               FUN = function(x) {
                                 if (x == "Negative") {
                                   return(colnames(barTable.n)[i])
                                 } else {
                                   return("Doublet")
                                 } })
  }

  return(bc_calls)

}

findReclassCellsFN <- function(barTable, neg.cells) {

  ## Normalize Data: Log2 Transform, mean-center
  print("Normalizing barcode data...")
  barTable.n <- as.data.frame(log2(barTable))
  for (i in 1:ncol(barTable)) {
    ind <- which(is.finite(barTable.n[,i]) == FALSE)
    barTable.n[ind,i] <- 0
    barTable.n[,i] <- barTable.n[,i]-mean(barTable.n[,i])
  }
  barTable.n_neg.cells <- barTable.n[neg.cells, ]

  ## Pre-allocate memory to save processing time/memory
  print("Pre-allocating data structures...")
  n_BC <- ncol(barTable.n)
  bcs <- colnames(barTable.n)
  n_neg <- length(neg.cells)
  bc.mats <- list()
  q.range <- seq(0.01, 0.99, by=0.02)
  ref.mat <- as.data.frame(matrix(0L, nrow=n_neg, ncol=length(q.range)))
  rownames(ref.mat) <- neg.cells
  colnames(ref.mat) <- paste("q", q.range, sep="_")
  for (i in 1:n_BC) {
    bc.mats[[i]] <- ref.mat
  }
  names(bc.mats) <- bcs

  ## Iterate through BCs, modeling distribution using GKDE
  ## Define low and high maxima before testing all negative cells across all q
  print("Determining classifications for negative cells across all BCs, all q...")
  for (i in 1:n_BC) {
    ## Step 1: GKDE
    model <- approxfun(bkde(barTable.n[,i], kernel="normal"))
    x <- seq(from=quantile(barTable.n[,i],0.001), to=quantile(barTable.n[,i],0.999), length.out=100)

    ## Step 2: Local maxima definition
    extrema <- localMaxima(model(x))

    ## Step 3: Select maxima, avoiding noisy extremes (if possible)
    ## Assumes negative cells are largest mode
    ## -> NEW: assumes positive cells are the second largest mode
    low.extreme <- extrema[which.max(model(x)[extrema])]
    # high.extreme <- max(extrema)
    extrema <- extrema[extrema != low.extreme]
    high.extreme <- extrema[which.max(model(x)[extrema])]

    ## Iterate through q, tracking classiication status for negative cells
    col.ind <- 0
    for (q in q.range) {
      col.ind <- col.ind + 1
      thresh <- quantile(c(x[high.extreme], x[low.extreme]), q)
      cell_i <- which(barTable.n_neg.cells[,i] >= thresh)
      bc.mats[[i]][cell_i, col.ind] <- 1
    }
  }

  ## Determine classification status across q.range for negative cells across all barcodes
  print("Computing classification stability...")
  bc.mat.final <- Reduce('+',bc.mats)

  ## Compute classification stability ~ i.e., number of q for which cell has 1 classification
  res <- as.data.frame(matrix(0L, nrow=n_neg, ncol=2))
  colnames(res) <- c("ClassStability","Reclassification")
  rownames(res) <- neg.cells
  res$Reclassification <- "Negative"
  for (cell in neg.cells) {
    res[cell,"ClassStability"] <- length(which(bc.mat.final[cell, ] == 1))
  }

  ## Extract classifications for cells above cluster.stability reclassification threshold
  # Define reclassifiable cells
  print(paste("Extracting rescued classifications...",sep=""))
  cells2reclass <- neg.cells[which(res$ClassStability >= 1)]

  # Find q at which reclassifiable cells have single classification, extract corresponding bc
  for (cell in cells2reclass) {
    q.ind <- max(which(bc.mat.final[cell, ] == 1))

    for (bar in 1:n_BC) {
      temp <- bc.mats[[bar]][cell, q.ind]
      if (temp == 1) { res[cell, "Reclassification"] <- bcs[bar] }
    }
  }

  # Return only cells with an available classification
  res <- res[-which(res$Reclassification == "Negative"), ]
  return(res)

}

