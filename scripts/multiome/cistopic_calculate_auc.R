
options(warn=-1)


############################## OPTIONS MENU ###################################

library("optparse")

# specify our desired options in a list
# by default OptionParser will add an help option equivalent to 
# make_option(c("-h", "--help"), action="store_true", default=FALSE, 
#               help="Show this help message and exit")

cat("\n")
option_list <- list(
  make_option("--cisTopic1", type="character", 
              help="[REQUIRED] cisTopic Object 1 path (must have same topic search space size as cisTopic2)"),
  make_option("--cisTopic2", type="character", 
              help="[REQUIRED] cisTopic Object 2 path (must have same topic search space size as cisTopic1)"),
  make_option("--metadata1", type="character",
              help="[REQUIRED] File containing Seurat metadata information for object cisTopic1"),
  make_option("--metadata2", type="character",
              help="[REQUIRED] File containing Seurat metadata information for object cisTopic2"),
  make_option("--clusters1", type="character",
              help="[REQUIRED] Column name of Seurat clustering for object cisTopic1"),
  make_option("--clusters2", type="character",
              help="[REQUIRED] Column name of Seurat clustering for object cisTopic2")
)

############################## PARSE OPTIONS ##################################


# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults
opt <- parse_args(OptionParser(option_list=option_list))

if (is.null(opt$cisTopic1)) {
  write("Option --cisTopic1 is required\nTry --help for help", stderr())
  q()
} else {
  CISTOPIC1 <- opt$cisTopic1
}


if (is.null(opt$cisTopic2)) {
  write("Option --cisTopic2 is required\nTry --help for help", stderr())
  q()
} else {
  CISTOPIC2 <- opt$cisTopic2
}

if (is.null(opt$metadata1)) {
  write("Option --metadata1 is required\nTry --help for help", stderr())
  q()
} else {
  METADATA1 <- opt$metadata1
}

if (is.null(opt$metadata2)) {
  write("Option --metadata2 is required\nTry --help for help", stderr())
  q()
} else {
  METADATA2 <- opt$metadata2
}

if (is.null(opt$clusters1)) {
  write("Option --clusters1 is required\nTry --help for help", stderr())
  q()
} else {
  CLUSTERS1 <- opt$clusters1
}

if (is.null(opt$clusters2)) {
  write("Option --clusters2 is required\nTry --help for help", stderr())
  q()
} else {
  CLUSTERS2 <- opt$clusters2
}

################################# EXECUTION ###################################
suppressPackageStartupMessages({library(Signac)
  library(cisTopic)
  library(ggplot2)
  library(pROC)
  library(ComplexHeatmap)})

topic.auc <- function(topic,labels) {
  return(sapply(
    X=levels(labels),
    FUN=function(x, topic, labels) auc(labels == x, topic, direction='<'),
    topic=topic,
    labels=labels
  ))
}

calculate.auc <- function(df, labels) {
  return(apply(
    X=df,
    MARGIN=1,
    FUN = function(x, labels) topic.auc(x, labels),
    labels=labels
  ))
}


CISTOPIC1.NAME <- tools::file_path_sans_ext(basename(CISTOPIC1))
CISTOPIC2.NAME <- tools::file_path_sans_ext(basename(CISTOPIC2))

## Load the two cisTopic objects
cat("##-------------------------------- cisTopic analysis --------------------------------------\n##\n")
cat("## FILE READING:...")
cisTopicObject1 <- readRDS(CISTOPIC1)
cisTopicObject2 <- readRDS(CISTOPIC2)
meta.1 <- read.csv2(METADATA1, row.names = 1)
meta.2 <- read.csv2(METADATA2, row.names = 1)
cat("done! \n")

for (topics in cisTopicObject1@calc.params$runCGSModels$topic) {
  cat(paste("## CISTOPIC RUN:", topics, "Topics\n"))
  cisTopicObject1 <- selectModel(cisTopicObject1, select = topics)

  # Extract Z-score cell embeddings

  cat("## AUC CALCULATION\n##\n")
  cat("## \tCalculating cells embeddings...")
  cell.embeddings.1 <- modelMatSelection(cisTopicObject1,
                                         target = "cell",
                                         method = "Z-score")
  rownames(cell.embeddings.1) <- paste("Topic", seq(1,topics))
  
  cat("done!\n##\n")
  
  cat(paste("## \tCalculating AUC on", CLUSTERS1, "for", CISTOPIC1.NAME,"\n##\n"))
  
  cistopic1.clustering.AUC <- calculate.auc(cell.embeddings.1, as.factor(meta.1[[CLUSTERS1]]))
  rownames(cistopic1.clustering.AUC) <- paste("Cluster", seq(0,nrow(cistopic1.clustering.AUC)-1))
  
  
  cat("## Writing files...\n##\n")
  write.csv2(cistopic1.clustering.AUC, 
             file=file.path(dirname(CISTOPIC1), paste(CISTOPIC1.NAME, topics, CLUSTERS1, "AUC.csv", sep="_")),
             quote = F, dec = ".")
  
  cat("## Etching plots...\n##\n")
  png(file.path(dirname(CISTOPIC1), paste(CISTOPIC1.NAME, topics,CLUSTERS1,"AUC_heatmap.png", sep="_")), 
      width=15, 
      height=20,
      units="cm",
      res=300)
  cistopic1.auc.heatmap <- Heatmap(
    t(cistopic1.clustering.AUC),
    name = "AUC",
    column_title = CISTOPIC1.NAME,
    column_names_gp =gpar(fontsize=8)
  )
  draw(cistopic1.auc.heatmap)
  dev.off()
  
  ###################### LOOP CISTOPIC2 ###########################
  for (topics2 in cisTopicObject2@calc.params$runCGSModels$topic) {
    cisTopicObject2 <- selectModel(cisTopicObject2, select = topics2)
    
    cell.embeddings.2 <- modelMatSelection(cisTopicObject2,
                                           target = "cell",
                                           method = "Z-score")
    rownames(cell.embeddings.2) <- paste("Topic", seq(1,topics2))
    
    cat(paste("## \tCalculating AUC on", CLUSTERS2, "for", CISTOPIC2.NAME,"\n##\n"))
    
    cistopic2.clustering.AUC <- calculate.auc(cell.embeddings.2, as.factor(meta.2[[CLUSTERS2]]))
    rownames(cistopic2.clustering.AUC) <- paste("Cluster", seq(0,nrow(cistopic2.clustering.AUC)-1))
    
    
    write.csv2(cistopic2.clustering.AUC, 
               file=file.path(dirname(CISTOPIC2), paste(CISTOPIC2.NAME, topics2, CLUSTERS2, "AUC.csv", sep="_")),
               quote = F, dec = ".")
    
    png(file.path(dirname(CISTOPIC2), paste(CISTOPIC2.NAME, topics2,CLUSTERS2,"AUC_heatmap.png", sep="_")), 
        width=15, 
        height=20,
        units="cm",
        res=300)
    cistopic2.auc.heatmap <- Heatmap(
      t(cistopic2.clustering.AUC),
      name = "AUC",
      column_title = CISTOPIC2.NAME,
      column_names_gp =gpar(fontsize=8)
    )
    draw(cistopic2.auc.heatmap)
    dev.off()
    
    # Unify the matching topic information with AUC information
    cat("## Merging AUC information with overlap dataset...\n##\n")
    overlapping.experiments <- read.csv2(file.path(dirname(CISTOPIC1), 
                                                   paste(CISTOPIC1.NAME, topics,
                                                         CISTOPIC2.NAME, topics2, 
                                                         "cross_regions_overlap_size_list.csv", sep="_")), header = TRUE)
    
    deltas.AUC <- abs(cistopic1.clustering.AUC[c(4,5,7),overlapping.experiments$cisTopic1] - 
                        cistopic2.clustering.AUC[c(4,5,6),overlapping.experiments$cisTopic2])
    cat("checkpoint\n\n")
    
    rownames(deltas.AUC) <- c("delta.S2", "delta.S3", "delta.S1")
    
    cat("checkpoint2\n\n")
    
    cat(paste(rownames(cistopic1.clustering.AUC), collapse = " "))
    cat("\n\n")
    cat(paste(rownames(cistopic1.clustering.AUC),"cisTopic1", collapse = "-"))
    cat("\n\n")
    cat(rownames(cistopic1.clustering.AUC)[1] == "Cluster 0")
    cat("\n\n")
    
    rownames(cistopic1.clustering.AUC) <- ifelse(rep(rownames(cistopic1.clustering.AUC)[1] == "Cluster 0", dim(cistopic1.clustering.AUC)[1]), 
                                                 paste(rownames(cistopic1.clustering.AUC),"cisTopic1"), 
                                                 rownames(cistopic1.clustering.AUC))
    rownames(cistopic2.clustering.AUC) <- ifelse(rep(rownames(cistopic2.clustering.AUC)[1] == "Cluster 0", dim(cistopic2.clustering.AUC)[1]), 
                                                 paste(rownames(cistopic2.clustering.AUC),"cisTopic2"), 
                                                 rownames(cistopic2.clustering.AUC))
    
    cat("checkpoint3\n\n")
    
    overlapping.experiments <- cbind(overlapping.experiments, 
                                     t(cistopic1.clustering.AUC)[overlapping.experiments$cisTopic1,], 
                                     t(cistopic2.clustering.AUC)[overlapping.experiments$cisTopic2,],
                                     t(deltas.AUC))
    overlapping.experiments <- overlapping.experiments[order(overlapping.experiments$overlap, decreasing = TRUE),]
    
    write.csv2(overlapping.experiments, file = file.path(dirname(CISTOPIC1), 
                                                         paste(CISTOPIC1.NAME, topics, 
                                                               CISTOPIC2.NAME, topics2, 
                                                               "cross_regions_sorted_overlap_size_AUC_list.csv", sep="_")))

  }
}















































