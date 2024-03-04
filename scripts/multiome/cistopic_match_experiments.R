#' This script takes as input a cisTopic object and returns a set of plots and dataframes
#'Topic matching
#' Overlap plots 
#' Self and cross experiments
#' Dataframe with overlap sizes 
#' Dataframe for each topic couple
#' Select topic couple?
#' dataframe with regions, marker proximity

# shopping list
# cisTopicObject1
# cisTopicObject2
# marker gene list with coords in bed format containing 5 columns
# chrName chrStart chrEnd name cluster
#options(warn = -1)

paste5 <- function(..., sep = " ", collapse = NULL, na.rm = F) {
  if (na.rm == F)
    paste(..., sep = sep, collapse = collapse)
  else
    if (na.rm == T) {
      paste.na <- function(x, sep) {
        x <- gsub("^\\s+|\\s+$", "", x)
        ret <- paste(na.omit(x), collapse = sep)
        is.na(ret) <- ret == ""
        return(ret)
      }
      df <- data.frame(..., stringsAsFactors = F)
      ret <- apply(df, 1, FUN = function(x) paste.na(x, sep))
      
      if (is.null(collapse))
        ret
      else {
        paste.na(ret, sep = collapse)
      }
    }
}

add.ranking.column<- function(df) {
  df$rank <- order(df[[1]], decreasing = TRUE)
  return(df)
}

get.topics.mixed.ranks <- function(
    cisTopic1,
    cisTopic2,
    topic1,
    topic2) {
  topic1.binarized.ranks <- lapply(cisTopic1@binarized.cisTopics, add.ranking.column)
  topic2.binarized.ranks <- lapply(cisTopic2@binarized.cisTopics, add.ranking.column)
  topic1.regions <- StringToGRanges(rownames(topic1.binarized.ranks[[topic1]]), sep=c(":","-"))
  topic2.regions <- StringToGRanges(rownames(topic2.binarized.ranks[[topic2]]), sep=c(":","-"))
  overlapping.indexes <- findOverlaps(query = topic1.regions, subject = topic2.regions)
  
  topic1.ranks <- data.frame(scores1=cisTopic1@binarized.cisTopics[[topic1]][queryHits(overlapping.indexes),])
  topic1.ranks$regions1 <- rownames(cisTopic1@binarized.cisTopics[[topic1]])[queryHits(overlapping.indexes)]
  topic1.ranks$rank1 <- topic1.binarized.ranks[[topic1]][queryHits(overlapping.indexes), "rank"]
  topic2.ranks <- data.frame(scores2=cisTopic2@binarized.cisTopics[[topic2]][subjectHits(overlapping.indexes),])
  topic2.ranks$regions2 <- rownames(cisTopic2@binarized.cisTopics[[topic2]])[subjectHits(overlapping.indexes)]
  topic2.ranks$rank2 <- topic2.binarized.ranks[[topic2]][subjectHits(overlapping.indexes),"rank"]
  mean.scores <- apply(FUN=mean,MARGIN = 1, X = cbind(topic1.ranks$scores1,topic2.ranks$scores2))
  mean.rank <- apply(FUN=mean, MARGIN = 1, X = cbind(topic1.ranks$rank1, topic2.ranks$rank2))
  res <- cbind(topic1.ranks,topic2.ranks, mean.scores, mean.rank)
  
  return(res)
}

get.closest.marker.features <- function(regions, signatures, gapsize=2.5e+5, ignore.strand=TRUE) {
  # for each gene, measure overlap with an arbitrary gap distance
  # return gene names and whether they fall within one of the two clusters
  regions <- StringToGRanges(regions, sep=c(":","-"))
  markers.coords <- StringToGRanges(signatures$coords, sep=c(":","-"))
  elementMetadata(markers.coords)$gene_name <- signatures$genes
  regions.markers.overlap <- findOverlaps(query=regions, 
                                          subject=markers.coords,
                                          maxgap=gapsize)
  
  if (length(regions.markers.overlap) == 0) {
    empty.df <- data.frame(list(NA,NA,NA))
    return(empty.df)
  }
  
  hits.regions <- regions[queryHits(regions.markers.overlap)]
  hits.regions$gene_name <- markers.coords[subjectHits(regions.markers.overlap)]$gene_name
  
  hits.regions.df <- cbind(GRangesToString(hits.regions, sep=c(":","-")), hits.regions$gene_name)
  colnames(hits.regions.df) <- c("coords", "marker")
  results <- aggregate(marker~coords,data=hits.regions.df, paste, collapse=",")
  
  # now we need to find the cluster combination for each region associated with it
  
  results$cluster <- sapply(results$marker, 
                            function(x) paste0(subset(signatures, genes %in% strsplit(x, split=",")[[1]])$cluster, collapse=","))
  
  return(results)
}


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
  make_option("--markers_file", type="character",
              help="[REQUIRED] csv file containing marker genes genomic coordinates"),
  make_option("--thrP", type="numeric", default = 0.95,
              help="[OPTIONAL] Probability threshold for the Gamma Fit of regions for each topic in the binarization step [default=%default]"),
  make_option("--metadata1", type="character",
              help="[REQUIRED] File containing Seurat metadata information for object cisTopic1"),
  make_option("--metadata2", type="character",
              help="[REQUIRED] File containing Seurat metadata information for object cisTopic2")
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

if (is.null(opt$markers_file)) {
  write("Option --markers_file is required\nTry --help for help", stderr())
  q()
} else {
  MARKERS.FILE <- opt$markers_file
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

thrP <- opt$thrP

################################# EXECUTION ###################################
library(Signac)
library(cisTopic)
library(ComplexHeatmap)
library(GenomicRanges)
library(reshape2)
library(gtools)

CISTOPIC1.NAME <- tools::file_path_sans_ext(basename(CISTOPIC1))
CISTOPIC2.NAME <- tools::file_path_sans_ext(basename(CISTOPIC2))

## Load the two cisTopic objects
cat("##-------------------------------- cisTopic analysis --------------------------------------\n##\n")
cat("## FILE READING:...")
cisTopicObject1 <- readRDS(CISTOPIC1)
cisTopicObject2 <- readRDS(CISTOPIC2)
cat("done! \n")
# for each run (with different topic number)
# ASSUMING BOTH EXPERIMENT HAVE THE SAME TOPIC SEARCH SPACE SIZE
for (topics in cisTopicObject1@calc.params$runCGSModels$topic) {
  cat(paste("## CISTOPIC RUN:", topics, "Topics for cisTopic1\n"))
  
  png(file.path(dirname(CISTOPIC1), paste(CISTOPIC1.NAME, topics,"model_selection.png", sep="_")), 
      width=20, 
      height=10,
      units="cm",
      res=300)
  par(mfrow=c(1,2))
  cisTopicObject1 <- selectModel(cisTopicObject1, select = topics)
  dev.off()
  
  png(file.path(dirname(CISTOPIC2), paste(CISTOPIC2.NAME, topics,"model_selection.png", sep="_")), 
      width=20, 
      height=10,
      units="cm",
      res=300)
  par(mfrow=c(1,2))
  cisTopicObject2 <- selectModel(cisTopicObject2, select = topics)
  dev.off()
  cat("## \tSelection models plots saved!\n##\n")
  # Extract Z-score cell embeddings
  par(mfrow=c(1,1))
  cat("## OVERLAP ANALYSIS:\n##\n")
  cat("## \tComputing cell-embedding cased correlations...")
  cell.embeddings.1 <- modelMatSelection(cisTopicObject1,
                                        target = "cell",
                                        method = "Z-score")
  rownames(cell.embeddings.1) <- paste("Topic", seq(1,topics))
  
  cell.embeddings.2 <- modelMatSelection(cisTopicObject2,
                                        target = "cell",
                                        method = "Z-score")
  rownames(cell.embeddings.2) <- paste("Topic", seq(1,topics))
  
  
  # Calculate cell-embeddings based (self) topic Pearson's correlation coefficient
  
  cor.cells.1 <- cor(t(cell.embeddings.1))
  
  png(file.path(dirname(CISTOPIC1), paste(CISTOPIC1.NAME, topics,"topic_correlation_matrix_cellbased.png", sep="_")), 
      width=16, 
      height=15,
      units="cm",
      res=300)
  Heatmap(cor.cells.1,
          column_title=paste(CISTOPIC1.NAME, "Topic Correlation Matrix (by cell-embedding)"),
          row_names_gp = gpar(fontsize=10),
          column_names_gp = gpar(fontsize=10),heatmap_legend_param = list("title"="Pearson's R")
  )
  dev.off()
  
  cor.cells.2 <- cor(t(cell.embeddings.2))
  cat("...done...")
  png(file.path(dirname(CISTOPIC2), paste(CISTOPIC2.NAME, topics,"topic_correlation_matrix_cellbased.png", sep="_")), 
      width=16, 
      height=15,
      units="cm",
      res=300)
  Heatmap(cor.cells.2,
          column_title=paste(CISTOPIC1.NAME, "Topic Correlation Matrix (by cell-embedding)"),
          row_names_gp = gpar(fontsize=10),
          column_names_gp = gpar(fontsize=10),heatmap_legend_param = list("title"="Pearson's R")
  )
  dev.off()
  cat("complete!\n##\n")
  
  # Calculate regions scores for each topic
  cat("## \tCalculating regions scores and selecting top regions...")
  cisTopicObject1 <- getRegionsScores(cisTopicObject1, method = "NormTop")
  
  pdf(file=file.path(dirname(CISTOPIC1), paste(CISTOPIC1.NAME, topics,"thrP",thrP, "topics_selected_regions_fit_histogram.pdf", sep="_")), 
      width=10, 
      height=8,
      onefile=TRUE)
  cisTopicObject1 <- binarizecisTopics(cisTopicObject1, method = "GammaFit", thrP = thrP, plot=TRUE)
  dev.off()
  
  cat("...done...")
  
  regions.names.1 <- lapply(cisTopicObject1@binarized.cisTopics, 
                            FUN = function(x) StringToGRanges(regions=rownames(x), sep=c(":","-")))
  
  ## Compare topics in pairs of two across the experiments and count the overlaps for each topic
  combinations.topics <- gtools::permutations(topics,2,repeats.allowed = TRUE)
  cat("cisTopic1...")
  ## Self-comparisons
  common.regions.1 <- apply(X=combinations.topics, 
                            MARGIN = 1, 
                            FUN = function(ind) sum(countOverlaps(regions.names.1[[ind[1]]],
                                                                  regions.names.1[[ind[2]]])))
  
  df.1 <- data.frame(combinations.topics, common.regions.1)
  colnames(df.1) <- c("Topics1", "Topics2", "overlap")
  
  write.csv2(df.1, file = file.path(dirname(CISTOPIC1), paste(CISTOPIC1.NAME, topics,"self_regions_overlap_size_list.csv", sep="_")),
             quote = F, row.names = F)
  ## We also produce the same dataframes in matrix form, as we need them for plotting
  df.1.cast <- reshape2::dcast(data=df.1, 
                               formula = Topics1 ~ Topics2, 
                               value.var = "overlap")
  rownames(df.1.cast) <- paste("Topic", df.1.cast$Topics1)
  df.1.cast$Topics1 <- NULL
  colnames(df.1.cast) <- paste("Topic", colnames(df.1.cast))
  diag(df.1.cast) <- NA
  
  write.csv2(df.1.cast, file = file.path(dirname(CISTOPIC1), paste(CISTOPIC1.NAME, topics,"self_regions_overlap_size_matrix.csv", sep="_")),
             quote = F)
  
  png(filename=file.path(dirname(CISTOPIC1), paste(CISTOPIC1.NAME, topics,"self_topics_overlap_heatmaps.png", sep="_")), 
      width=18, 
      height=15,
      units="cm",
      res=300)
  # cisTopicObject1 topics self-overlap heatmap
  df.1.heatmap <- Heatmap(as.matrix(df.1.cast),
                          row_names_gp = gpar(fontsize=10),
                          column_names_gp = gpar(fontsize=10),
                          heatmap_legend_param = list("title"="Overlapping size"))
  draw(df.1.heatmap)
  dev.off()
  
  #################### LOOP CISTOPIC2 ##########################################
  for (topics2 in cisTopicObject1@calc.params$runCGSModels$topic) {
    # we reset the number of topics and do the compared analysis so to match all combinations of topics between experiments
    cisTopicObject2 <- selectModel(cisTopicObject2, select = topics2)
    cell.embeddings.2 <- modelMatSelection(cisTopicObject2,
                                           target = "cell",
                                           method = "Z-score")
    rownames(cell.embeddings.2) <- paste("Topic", seq(1,topics2))
    
    cisTopicObject2 <- getRegionsScores(cisTopicObject2, method = "NormTop")
    
    pdf(file=file.path(dirname(CISTOPIC2), paste(CISTOPIC2.NAME, topics2,"thrP",thrP,"topics_selected_regions_fit_histogram.pdf", sep="_")), 
        width=10, 
        height=8,
        onefile=TRUE)
    cisTopicObject2 <- binarizecisTopics(cisTopicObject2, method = "GammaFit", thrP = thrP, plot=TRUE)
    dev.off()
    cat("complete!\n##\n")
    # Calculate overlap sizes across all topic pairs
    
    ## Extract region names and transform them into genomic ranges
    cat("## \tCalculating overlap size between topics...")
    regions.names.2 <- lapply(cisTopicObject2@binarized.cisTopics, 
                               FUN = function(x) StringToGRanges(regions=rownames(x), sep=c(":","-")))
    
    ## Compare topics in pairs of two across the two experiments and count the overlaps for each topic
    combinations.topics.2 <- gtools::permutations(topics2,2,repeats.allowed = TRUE)
    combinations.topics.cross <- expand.grid(seq(topics),seq(topics2))
    ## Self-comparisons
    cat("cisTopic2...")
    common.regions.2 <- apply(X=combinations.topics.2, 
                               MARGIN = 1, 
                               FUN = function(ind) sum(countOverlaps(regions.names.2[[ind[1]]],
                                                                     regions.names.2[[ind[2]]])))
    
    # Cross-comparisons
    cat("cross-experiment.\n##\n")
    common.regions.cross <- apply(X=combinations.topics.cross, 
                                  MARGIN = 1, 
                                  FUN = function(ind) sum(countOverlaps(regions.names.1[[ind[1]]], 
                                                                        regions.names.2[[ind[2]]])))
    
    ## Pack informations into dataframes
    cat("## \tPacking overlap information into dataframes (list and matrix)...")
    df.2 <- data.frame(combinations.topics.2, common.regions.2)
    colnames(df.2) <- c("Topics1", "Topics2", "overlap")
    df.cross <- data.frame(combinations.topics.cross, common.regions.cross)
    colnames(df.cross) <- c("cisTopic1", "cisTopic2", "overlap")
    
    ## We also produce the same dataframes in matrix form, as we need them for plotting
    df.2.cast <- reshape2::dcast(data=df.2, 
                                 formula = Topics1 ~ Topics2, 
                                 value.var = "overlap")
    rownames(df.2.cast) <- paste("Topic", df.2.cast$Topics1)
    df.2.cast$Topics1 <- NULL
    colnames(df.2.cast) <- paste("Topic", colnames(df.2.cast))
    diag(df.2.cast) <- NA
    
    df.cross.cast <- reshape2::dcast(data=df.cross, 
                                     formula = cisTopic1 ~ cisTopic2, 
                                      value.var = "overlap")
    rownames(df.cross.cast) <- paste("Topic", df.cross.cast$cisTopic1)
    df.cross.cast$cisTopic1 <- NULL
    colnames(df.cross.cast) <- paste("Topic", colnames(df.cross.cast))
    
    
    write.csv2(df.2, file = file.path(dirname(CISTOPIC2), paste(CISTOPIC2.NAME, topics2,"self_regions_overlap_size_list.csv", sep="_")),
               quote = F, row.names = F)
    write.csv2(df.2.cast, file = file.path(dirname(CISTOPIC2), paste(CISTOPIC2.NAME, topics2,"self_regions_overlap_size_matrix.csv", sep="_")),
               quote = F)
    
    write.csv2(df.cross, file = file.path(dirname(CISTOPIC1), paste(CISTOPIC1.NAME, topics, CISTOPIC2.NAME, topics2,"cross_regions_overlap_size_list.csv", sep="_")),
               quote = F, row.names = F)
    write.csv2(df.cross.cast, file = file.path(dirname(CISTOPIC1), paste(CISTOPIC1.NAME, topics, CISTOPIC2.NAME, topics2,"cross_regions_overlap_size_matrix.csv", sep="_")),
               quote = F)
    
    cat("done!\n##\n")
    cat("## \tProducing overlap heatmaps (self and cross experiment)...")
    
    png(filename=file.path(dirname(CISTOPIC2), paste(CISTOPIC2.NAME, topics2,"self_topics_overlap_heatmaps.png", sep="_")), 
        width=18, 
        height=15,
        units="cm",
        res=300)
    # cisTopicObject2 topics self-overlap heatmap
    df.2.heatmap <- Heatmap(as.matrix(df.2.cast),
            row_names_gp = gpar(fontsize=10),
            column_names_gp = gpar(fontsize=10),
            heatmap_legend_param = list("title"="Overlapping size"))
    draw(df.2.heatmap)
    dev.off()
    
    # Add annotations for nCount_ATAC and nFeature_RNA
    
    meta.1 <- read.csv2(METADATA1)
    meta.2 <- read.csv2(METADATA2)
    
    corr_abs_ATAC_1 <- abs(cor(t(cell.embeddings.1), meta.1$nCount_ATAC)) # topics
    corr_abs_ATAC_2 <- abs(cor(t(cell.embeddings.2), meta.2$nCount_ATAC)) # topics2
    corr_abs_RNA_1 <- abs(cor(t(cell.embeddings.1), meta.1$nFeature_RNA)) # topics
    corr_abs_RNA_2 <- abs(cor(t(cell.embeddings.2), meta.2$nFeature_RNA)) # topics2
    
    cat(paste(length(corr_abs_ATAC_1),length(corr_abs_ATAC_2)))
    
    pch_1_ATAC <- rep("*",topics)
    pch_1_ATAC[corr_abs_ATAC_1 < 0.5] <- NA
    pch_2_ATAC <- rep("*",topics2)
    pch_2_ATAC[corr_abs_ATAC_2 < 0.5] <- NA
    
    pch_1_RNA <- rep("*",topics)
    pch_1_RNA[corr_abs_RNA_1 < 0.5] <- NA
    pch_2_RNA <- rep("*",topics2)
    pch_2_RNA[corr_abs_RNA_2 < 0.5] <- NA
    
    col_fun_RNA <- circlize::colorRamp2(c(0, 0.5, 1), c("green", "white", "red"))
    col_fun_ATAC <- circlize::colorRamp2(c(0, 0.5, 1), c("yellow", "white", "red"))
    
    anno_1 <- rowAnnotation(ATAC_Counts=anno_simple(corr_abs_ATAC_1, col=col_fun_ATAC, pch = pch_1_ATAC),
                            RNA_Features=anno_simple(corr_abs_RNA_1, col=col_fun_RNA, pch= pch_1_RNA),
                            annotation_name_side="bottom")
    
    anno_2 <- HeatmapAnnotation(ATAC_Counts=anno_simple(corr_abs_ATAC_2, col=col_fun_ATAC, pch = pch_2_ATAC),
                            RNA_Features=anno_simple(corr_abs_RNA_2, col=col_fun_RNA, pch= pch_2_RNA),
                            annotation_name_side="left")
    
    
    # cisTopicObject1 and cisTopicObject2 topics cross-overlap heatmap
    heatmap.cross <- Heatmap(as.matrix(df.cross.cast),column_title = CISTOPIC2.NAME,
                             row_title = CISTOPIC1.NAME,
                             heatmap_legend_param = list("title"="Overlap"),
                             row_km = 2,
                             column_km = 2,
                             bottom_annotation = anno_2,
                             right_annotation= anno_1,
                             column_names_rot = 45,
                             row_names_gp = gpar(fontsize = 10),
                             column_names_gp = gpar(fontsize = 9))
    
    lgd_corr_ATAC <- Legend(title = "Abs. correlation \nnCount_ATAC", col_fun = col_fun_ATAC, at=c(0,0.5,1))
    lgd_corr_RNA <- Legend(title = "Abs. correlation \nnFeatures_RNA", col_fun = col_fun_RNA, at=c(0,0.5,1))
    lgd_thresh <- Legend(pch = "*", type = "points", labels = "> 0.5")
    
    png(filename=file.path(dirname(CISTOPIC1), paste(CISTOPIC1.NAME,topics, CISTOPIC2.NAME, topics2,"cross_topics_overlap_heatmaps.png", sep="_")), 
        width=18, 
        height=15,
        units="cm",
        res=300)
    draw(heatmap.cross, annotation_legend_list=list(lgd_corr_ATAC, lgd_corr_RNA, lgd_thresh))
    dev.off()
    cat("done!\n##\n")
    cat("## \tPlotting correlation between Topic cell-embedding and nCount_ATAC...")
    corr_ATAC_1 <- cor(t(cell.embeddings.1), meta.1$nCount_ATAC)
    names(corr_ATAC_1) <- paste("Topic", seq(1,topics))
    corr_ATAC_2 <- cor(t(cell.embeddings.2), meta.2$nCount_ATAC)
    names(corr_ATAC_2) <- paste("Topic", seq(1,topics2))
    
    # Plot topic correlation with nCount_ATAC
    
    png(filename = file.path(dirname(CISTOPIC1), paste(CISTOPIC1.NAME, topics,"correlation_nCount_ATAC_barplot.png", sep="_")),
        width=20, 
        height=15,
        units="cm",
        res=300)
    barplot(corr_ATAC_1[order(abs(corr_ATAC_1))], las=2, col="cornflowerblue", main=paste(CISTOPIC1.NAME, "Topics correlation with nCount_ATAC"))
    dev.off()
    
    png(filename = file.path(dirname(CISTOPIC2), paste(CISTOPIC2.NAME, topics2,"correlation_nCount_ATAC_barplot.png", sep="_")),
        width=18, 
        height=15,
        units="cm",
        res=300)
    barplot(corr_ATAC_2[order(abs(corr_ATAC_2))], 
            las=2, 
            col="cornflowerblue", 
            main=paste(CISTOPIC2.NAME, "Topics correlation with nCount_ATAC")
            )
    dev.off()
    
    cat("done!\n##\n")
    ################## Marker distance evaluation ##################
    # we first read the signatures file
    cat("## MARKER DISTANCE EVALUATION:\n##\n")
    signatures <- read.csv2(MARKERS.FILE)
    
    plots_dir <- paste(paste(dirname(CISTOPIC1),"overlapping_regions","",sep='/'))
    if (!dir.exists(plots_dir)){
      dir.create(plots_dir, recursive = T)
    }
    
    # Extract regions for each topic pair
    for (idx in seq(1,length(df.cross$overlap))) {
      topic1 <- df.cross$cisTopic1[[idx]]
      topic2 <- df.cross$cisTopic2[[idx]]
      cat(paste("## \tMatching topic", topic1, "&", topic2,"(Overlap size:", df.cross$overlap[idx],")","\n"))
      # obtain a dataframe with regions with respective ranks and scores and a mean.rank
      if (df.cross$overlap[idx] == 0) {
        next
      }
      matching.regions <- get.topics.mixed.ranks(cisTopicObject1, cisTopicObject2, topic1, topic2)
      
      # Report closest marker features and associated clusters
      marker.features.500k <- get.closest.marker.features(regions = matching.regions$regions1,
                                                          signatures = signatures,
                                                          gapsize = 5e+5)
      colnames(marker.features.500k) <- c("coords", "marker_500kb", "cluster_500kb")
      
      marker.features.50k <- get.closest.marker.features(regions = matching.regions$regions1,
                                                          signatures = signatures,
                                                          gapsize = 5e+4)
      colnames(marker.features.50k) <- c("coords", "marker_50kb", "cluster_50kb")
      
      marker.features.5k <- get.closest.marker.features(regions = matching.regions$regions1,
                                                          signatures = signatures,
                                                          gapsize = 5e+3)
      colnames(marker.features.5k) <- c("coords", "marker_5kb", "cluster_5kb")
      
      # merge features into the main topic-paired regions dataframe 
      merged.regions.marker.df <- merge(x=matching.regions,
                                       y=marker.features.500k,
                                       all.x=TRUE, 
                                       by.x="regions1",
                                       by.y="coords")
      merged.regions.marker.df <- merge(x=merged.regions.marker.df,
                                        y=marker.features.50k,
                                        all.x=TRUE, 
                                        by.x="regions1",
                                        by.y="coords")
      merged.regions.marker.df <- merge(x=merged.regions.marker.df,
                                        y=marker.features.5k,
                                        all.x=TRUE, 
                                        by.x="regions1",
                                        by.y="coords")
      write.csv2(merged.regions.marker.df, file = file.path(dirname(CISTOPIC1), 
                                                            "overlapping_regions",
                                                            paste(CISTOPIC1.NAME, topics,
                                                                  CISTOPIC2.NAME, topics2,"Topic",topic1, "Topic",topic2, "matching_regions.csv", 
                                                                  sep="_")),
                 quote = F, row.names = F, dec = ".")
      }
  }
}

