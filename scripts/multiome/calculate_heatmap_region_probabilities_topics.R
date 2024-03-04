#'
#' This script is structured as following
#' - overlaps regions found in the IDR with regions in the cisTopic input
#' - 
#'
#'
#'
#'
#'


library(cisTopic)
library(GenomicRanges)
library(tidyverse)
library(Signac)
library(ComplexHeatmap)

#setwd("~/Documents/multiome/")

cisTopic1C <- readRDS("/hpcnfs/FN/data/sprocaccia/cisTopic/cisTopic_all/1C_ARC_seurat_cisTopic.Rds")
cisTopic1E <- readRDS("/hpcnfs/FN/data/sprocaccia/cisTopic/cisTopic_all/1E_ARC_seurat_cisTopic.Rds")

cisTopic1C <- selectModel(cisTopic1C, type="maximum")
cisTopic1E <- selectModel(cisTopic1E, type="maximum")

cisTopic1C.region.emb <- modelMatSelection(cisTopic1C,
                                        target = "region",
                                        method = "Probability",
                                        all.regions = TRUE)
cisTopic1E.region.emb <- modelMatSelection(cisTopic1E,
                                           target = "region",
                                           method = "Probability",
                                           all.regions = TRUE)

topic.2.topic.16 <- read.csv2("/hpcnfs/FN/data/sprocaccia/cisTopic/cisTopic_all/1C_ARC_seurat_cisTopic_1E_ARC_seurat_cisTopic_2_16_IDR_0.05_with_markers.csv")
topic.15.topic.19 <- read.csv2("/hpcnfs/FN/data/sprocaccia/cisTopic/cisTopic_all/1C_ARC_seurat_cisTopic_1E_ARC_seurat_cisTopic_15_19_IDR_0.05_with_markers.csv")
topic.27.topic.14 <- read.csv2("/hpcnfs/FN/data/sprocaccia/cisTopic/cisTopic_all/1C_ARC_seurat_cisTopic_1E_ARC_seurat_cisTopic_27_14_IDR_0.05_with_markers.csv")
topic.31.topic.15 <- read.csv2("/hpcnfs/FN/data/sprocaccia/cisTopic/cisTopic_all/1C_ARC_seurat_cisTopic_1E_ARC_seurat_cisTopic_31_15_IDR_0.05_with_markers.csv")

topic.2.topic.16 %>%
  subset(gIDR > -log10(0.05)) %>%
  top_n(n=40, wt=gIDR) %>% 
  arrange(desc(gIDR)) %>%
  makeGRangesFromDataFrame() %>%
  findOverlaps(subject=cisTopic1C@region.ranges) -> topic.2.topic.16.top.idx

topic.27.topic.14 %>%
  subset(gIDR > -log10(0.05)) %>%
  top_n(n=40, wt=gIDR) %>% 
  arrange(desc(gIDR))%>%
  makeGRangesFromDataFrame() %>%
  findOverlaps(subject=cisTopic1C@region.ranges) -> topic.27.topic.14.top.idx

topic.15.topic.19 %>%
  subset(gIDR > -log10(0.05)) %>%
  top_n(n=40, wt=gIDR) %>% 
  arrange(desc(gIDR))%>%
  makeGRangesFromDataFrame() %>%
  findOverlaps(subject=cisTopic1C@region.ranges) -> topic.15.topic.19.top.idx

topic.31.topic.15 %>%
  subset(gIDR > -log10(0.05)) %>%
  top_n(n=40, wt=gIDR) %>% 
  arrange(desc(gIDR))%>%
  makeGRangesFromDataFrame() %>%
  findOverlaps(subject=cisTopic1C@region.ranges) -> topic.31.topic.15.top.idx

heatmap.features.idx <- c(topic.2.topic.16.top.idx@to, 
                          topic.15.topic.19.top.idx@to, 
                          topic.27.topic.14.top.idx@to,
                          topic.31.topic.15.top.idx@to)

topic.2.topic.16 %>%
  subset(gIDR > -log10(0.05)) %>%
  top_n(n=30, wt=gIDR) %>% 
  arrange(desc(gIDR))%>% select(marker_50kb) -> marker.2.16

topic.15.topic.19 %>%
  subset(gIDR > -log10(0.05)) %>%
  top_n(n=30, wt=gIDR)%>% 
  arrange(desc(gIDR)) %>% select(marker_50kb) -> marker.15.19

topic.27.topic.14 %>%
  subset(gIDR > -log10(0.05)) %>%
  top_n(n=30, wt=gIDR) %>% 
  arrange(desc(gIDR)) %>% select(marker_50kb) -> marker.27.14

topic.31.topic.15 %>%
  subset(gIDR > -log10(0.05)) %>%
  top_n(n=30, wt=gIDR) %>% 
  arrange(desc(gIDR)) %>% select(marker_50kb) -> marker.31.15

marker.names.anno.df <- c(marker.2.16[topic.2.topic.16.top.idx@from,1], 
                              marker.15.19[topic.15.topic.19.top.idx@from,1], 
                              marker.27.14[topic.27.topic.14.top.idx@from,1], 
                              marker.31.15[topic.31.topic.15.top.idx@from,1])

marker.names.anno.df[is.na(marker.names.anno.df)] <- ""
names(marker.names.anno.df) <- c("Marker")

heatmap.features <- GRangesToString(cisTopic1C@region.ranges[heatmap.features.idx])

heatmap.matrix <- cisTopic1C.region.emb[c(2,31,27,15),heatmap.features.idx]

rownames(heatmap.matrix) <- paste("Module",1:4)

signatures <- read.csv2("signatures.csv")

PAGE5.markers <- signatures[grep("PAGE5+", signatures$cluster),]
FEZ1.markers <- signatures[grep("FEZ1+", signatures$cluster),]
Cluster3.markers <- signatures[grep("Cluster3", signatures$cluster),]

total.regions.1C <- cisTopic1C@region.ranges

regions.1C.markers.overlap.PAGE5 <- findOverlaps(total.regions.1C, 
                                                 StringToGRanges(PAGE5.markers$coords, sep=c(":","-")),
                                                 maxgap = 50000)

regions.1C.markers.overlap.FEZ1 <- findOverlaps(total.regions.1C, 
                                                StringToGRanges(FEZ1.markers$coords, sep=c(":","-")),
                                                maxgap = 50000)

regions.1C.markers.overlap.Cluster3 <- findOverlaps(total.regions.1C, 
                                                    StringToGRanges(Cluster3.markers$coords, sep=c(":","-")),
                                                    maxgap = 50000)


total.regions.1C@elementMetadata$PAGE5 <- rep(FALSE, length(total.regions.1C))
total.regions.1C@elementMetadata$PAGE5[regions.1C.markers.overlap.PAGE5@from] <- TRUE

total.regions.1C@elementMetadata$FEZ1 <- rep(FALSE, length(total.regions.1C))
total.regions.1C@elementMetadata$FEZ1[regions.1C.markers.overlap.FEZ1@from] <- TRUE

total.regions.1C@elementMetadata$Cluster3 <- rep(FALSE, length(total.regions.1C))
total.regions.1C@elementMetadata$Cluster3[regions.1C.markers.overlap.Cluster3@from] <- TRUE

total.regions.1C@elementMetadata$markers <- rep("", length(total.regions.1C))
total.regions.1C@elementMetadata$markers[heatmap.features.idx] <- heatmap.features.idx


marker.anno <- rowAnnotation(df=total.regions.1C@elementMetadata[heatmap.features.idx, 
                                                              c("PAGE5","FEZ1","Cluster3")],
                             col=list(PAGE5=structure(c("grey90","#E41A1C"), names=c(FALSE,TRUE)),
                                      FEZ1=structure(c("grey90","#377EB8"), names=c(FALSE,TRUE)),
                                      Cluster3=structure(c("grey90","#4DAF4A"), names=c(FALSE,TRUE))))

marker.anno.names <- rowAnnotation(Marker=anno_mark(at=which(!is.na(marker.names.anno.df)), 
                               labels = marker.names.anno.df[!is.na(marker.names.anno.df)],
                              which="row",
                              labels_gp=gpar(fontsize=10),
                               side = "left" ))

png("heatmap_region_cistopic_probabilities.png",
    height=30,
    width = 15,
    units = "cm",
    res=400)
Heatmap(scale(t(as.matrix(heatmap.matrix))),
        cluster_rows = FALSE, cluster_columns = FALSE,
        name="Region\nprobability\nZ-score",
        row_names_gp = gpar(fontsize=6),
        right_annotation = marker.anno,
        left_annotation = marker.anno.names
        )
dev.off()

###################### Heatmap x poster ####################################

heatmap.features.idx <- c(topic.2.topic.16.top.idx@to,
                          topic.31.topic.15.top.idx@to,
                          topic.27.topic.14.top.idx@to,
                          topic.15.topic.19.top.idx@to
                          )


marker.names.anno.df <- c(marker.2.16[topic.2.topic.16.top.idx@from,1], 
                          marker.31.15[topic.31.topic.15.top.idx@from,1], 
                          marker.27.14[topic.27.topic.14.top.idx@from,1], 
                          marker.15.19[topic.15.topic.19.top.idx@from,1]
                          )

marker.names.anno.df[is.na(marker.names.anno.df)] <- ""

heatmap.features <- GRangesToString(cisTopic1C@region.ranges[heatmap.features.idx])

heatmap.matrix <- cisTopic1C.region.emb[c(2,31,27,15),heatmap.features.idx]

rownames(heatmap.matrix) <- paste("Module",1:4)

signatures <- read.csv2("signatures.csv")

PAGE5.markers <- signatures[grep("PAGE5+", signatures$cluster),]
FEZ1.markers <- signatures[grep("FEZ1+", signatures$cluster),]
Cluster3.markers <- signatures[grep("Cluster3", signatures$cluster),]

total.regions.1C <- cisTopic1C@region.ranges

regions.1C.markers.overlap.PAGE5 <- findOverlaps(total.regions.1C, 
                                                 StringToGRanges(PAGE5.markers$coords, sep=c(":","-")),
                                                 maxgap = 50000)

regions.1C.markers.overlap.FEZ1 <- findOverlaps(total.regions.1C, 
                                                StringToGRanges(FEZ1.markers$coords, sep=c(":","-")),
                                                maxgap = 50000)

regions.1C.markers.overlap.Cluster3 <- findOverlaps(total.regions.1C, 
                                                    StringToGRanges(Cluster3.markers$coords, sep=c(":","-")),
                                                    maxgap = 50000)


total.regions.1C@elementMetadata$S1 <- rep(FALSE, length(total.regions.1C))
total.regions.1C@elementMetadata$S1[regions.1C.markers.overlap.PAGE5@from] <- TRUE

total.regions.1C@elementMetadata$S3 <- rep(FALSE, length(total.regions.1C))
total.regions.1C@elementMetadata$S3[regions.1C.markers.overlap.FEZ1@from] <- TRUE

total.regions.1C@elementMetadata$S2 <- rep(FALSE, length(total.regions.1C))
total.regions.1C@elementMetadata$S2[regions.1C.markers.overlap.Cluster3@from] <- TRUE

total.regions.1C@elementMetadata$markers <- rep("", length(total.regions.1C))
total.regions.1C@elementMetadata$markers[heatmap.features.idx] <- heatmap.features.idx


marker.anno <- rowAnnotation(df=total.regions.1C@elementMetadata[heatmap.features.idx, 
                                                                 c("S1","S3","S2")],
                             annotation_name_gp= gpar(fontsize = 18),
                             col=list(S1=structure(c("grey90","#FE7EDF"), names=c(FALSE,TRUE)),
                                      S3=structure(c("grey90","#00C3EF"), names=c(FALSE,TRUE)),
                                      S2=structure(c("grey90","#00C9A5"), names=c(FALSE,TRUE))),
                             annotation_legend_param = list(
                               S1=list( 
                                 title_gp = gpar(fontsize = 18, 
                                                 fontface = "bold"), 
                                 labels_gp = gpar(fontsize = 18)),
                               S3=list( 
                                 title_gp = gpar(fontsize = 18, 
                                                 fontface = "bold"), 
                                 labels_gp = gpar(fontsize = 18)),
                               S2=list( 
                                 title_gp = gpar(fontsize = 18, 
                                                 fontface = "bold"), 
                                 labels_gp = gpar(fontsize = 18))
                             ),
                             show_legend=FALSE)

marker.anno.names <- rowAnnotation(Marker=anno_mark(at=which(marker.names.anno.df != ""), 
                                                    labels = sub(",", ",\n",marker.names.anno.df[marker.names.anno.df != ""]),
                                                    which="row",
                                                    labels_gp=gpar(fontsize=12),
                                                    side = "left" ))


png("heatmap_region_cistopic_probabilities_sorted_poster.png",
    height=26,
    width = 9,
    units = "cm",
    res=100, type="quartz")
h <- Heatmap(scale(t(as.matrix(heatmap.matrix))),
        cluster_rows = FALSE, cluster_columns = FALSE,
        show_row_names = FALSE, row_title = "Top 40 regions per epigenetic module", row_title_side = "right",
        name="Region probability Z-score",
        row_title_gp = gpar(fontsize=18),
        column_names_gp = gpar(fontsize=18),
        right_annotation = marker.anno,
        left_annotation = marker.anno.names,
        heatmap_legend_param = list(
            title_gp = gpar(fontsize = 18,
                            fontface = "bold"), 
            labels_gp = gpar(fontsize = 18),
            legend_direction = "horizontal",
            legend_width = unit(8, "cm"),
            title_position = "topcenter"
        )
)
ComplexHeatmap::draw(h, heatmap_legend_side="bottom", annotation_legend_side="bottom")
dev.off()



