library(Seurat)
library(Signac)
library(tidyverse)
library(GenomicRanges)

object1C <- readRDS("~/Documents/multiome/hpc/cisTopic_all/1C_ARC_seurat.Rds")
object1E <- readRDS("~/Documents/multiome/hpc/cisTopic_all/1E_ARC_seurat.Rds")

cisTopic1C <- readRDS("~/Documents/multiome/hpc/cisTopic_all/1C_ARC_seurat_cisTopic.Rds")
cisTopic1E <- readRDS("~/Documents/multiome/hpc/cisTopic_all/1E_ARC_seurat_cisTopic.Rds")

# setwd("~/Documents/multiome")

topic.2.topic.16 <- read.csv2("hpc/cisTopic_all/1C_ARC_seurat_cisTopic_1E_ARC_seurat_cisTopic_2_16_IDR_0.05_with_markers.csv")
topic.15.topic.19 <- read.csv2("hpc/cisTopic_all/1C_ARC_seurat_cisTopic_1E_ARC_seurat_cisTopic_15_19_IDR_0.05_with_markers.csv")
topic.27.topic.14 <- read.csv2("hpc/cisTopic_all/1C_ARC_seurat_cisTopic_1E_ARC_seurat_cisTopic_27_14_IDR_0.05_with_markers.csv")
topic.31.topic.15 <- read.csv2("hpc/cisTopic_all/1C_ARC_seurat_cisTopic_1E_ARC_seurat_cisTopic_31_15_IDR_0.05_with_markers.csv")

topic.2.topic.16 %>%
  subset(gIDR > -log10(0.05)) %>%
  top_n(n=200, wt=gIDR) %>% 
  arrange(desc(gIDR)) %>%
  makeGRangesFromDataFrame() %>%
  findOverlaps(subject=cisTopic1C@region.ranges) -> topic.2.topic.16.top.idx

topic.27.topic.14 %>%
  subset(gIDR > -log10(0.05)) %>%
  top_n(n=200, wt=gIDR) %>% 
  arrange(desc(gIDR))%>%
  makeGRangesFromDataFrame() %>%
  findOverlaps(subject=cisTopic1C@region.ranges) -> topic.27.topic.14.top.idx

topic.15.topic.19 %>%
  subset(gIDR > -log10(0.05)) %>%
  top_n(n=30, wt=gIDR) %>% 
  arrange(desc(gIDR))%>%
  makeGRangesFromDataFrame() %>%
  findOverlaps(subject=cisTopic1C@region.ranges) -> topic.15.topic.19.top.idx

topic.31.topic.15 %>%
  subset(gIDR > -log10(0.05)) %>%
  top_n(n=200, wt=gIDR) %>% 
  arrange(desc(gIDR))%>%
  makeGRangesFromDataFrame() %>%
  findOverlaps(subject=cisTopic1C@region.ranges) -> topic.31.topic.15.top.idx

png("umap_top_regions_topics_tfidf_1C_topic_2_16.png",
    width = 50,
    height = 32, 
    units = "cm",
    res=600)
FeaturePlot(object1C, features=GRangesToString(cisTopic1C@region.ranges[topic.2.topic.16.top.idx@to]), 
            reduction = "umap_rna", 
            slot="data", 
            ncol = 5,
            max.cutoff = "q95", order=TRUE) * labs(color="TF-IDF")
dev.off()

png("umap_top_regions_topics_tfidf_1C_topic_27_14.png",
    width = 50,
    height = 32, 
    units = "cm",
    res=600)
FeaturePlot(object1C, features=GRangesToString(cisTopic1C@region.ranges[topic.27.topic.14.top.idx@to[1:20]]), 
            reduction = "umap_rna", 
            slot="data",  
            ncol = 5,
            max.cutoff = "q95", order=TRUE) * labs(color="TF-IDF")
dev.off()

png("umap_top_regions_topics_tfidf_1C_topic_15_19.png",
    width = 50,
    height = 32, 
    units = "cm",
    res=600)
FeaturePlot(object1C, features=GRangesToString(cisTopic1C@region.ranges[topic.15.topic.19.top.idx@to]), 
            reduction = "umap_rna", 
            slot="data", 
            ncol = 5, 
            max.cutoff = "q95", order=TRUE) * labs(color="TF-IDF")
dev.off()

png("umap_top_regions_topics_tfidf_1C_topic_31_15.png",
    width = 50,
    height = 32, 
    units = "cm",
    res=600)
FeaturePlot(object1C, features=GRangesToString(cisTopic1C@region.ranges[topic.31.topic.15.top.idx@to[1:20]]), 
            reduction = "umap_rna", 
            slot="data",  
            ncol = 5,
            max.cutoff = "q95", order=TRUE) * labs(color="TF-IDF")
dev.off()


################ UMAP selection for poster ####################################

png("poster/region_tfidf_2_16.png",
    height = 7,
    width = 9,
    units = "cm",
    res=100)
# 2/16
FeaturePlot(object1C, features=GRangesToString(cisTopic1C@region.ranges[topic.2.topic.16.top.idx@to[29]]), 
            reduction = "umap_rna", 
            slot="data", 
            max.cutoff = "q95", order=TRUE) + theme_void() +
  labs(color="TF-IDF", title = "") + 
  theme(legend.title = element_text(size=18), legend.text = element_text(size=18))
dev.off()

png("poster/region_tfidf_2_16_marker.png",
    height = 7,
    width = 7,
    units = "cm",
    res=100)
# 2/16
FeaturePlot(object1C, features=GRangesToString(cisTopic1C@region.ranges[topic.2.topic.16.top.idx@to[67]]), 
            reduction = "umap_rna", 
            slot="data", 
            max.cutoff = "q95", order=TRUE) + theme_void() +
  labs( title = "") + 
  theme(legend.position="none")
dev.off()

###### 27/14
png("poster/region_tfidf_27_14.png",
    height = 7,
    width = 9,
    units = "cm",
    res=100)
FeaturePlot(object1C, features="chr10-27959170-27960160", 
            reduction = "umap_rna", 
            slot="data", 
            max.cutoff = "q95", order=TRUE) + theme_void() +
  labs(color="TF-IDF", title = "") + 
  theme(legend.title = element_text(size=18), legend.text = element_text(size=18))
dev.off()

png("poster/region_tfidf_27_14_marker.png",
    height = 7,
    width = 7,
    units = "cm",
    res=100)
FeaturePlot(object1C, features=GRangesToString(cisTopic1C@region.ranges[topic.27.topic.14.top.idx@to[12]]), 
            reduction = "umap_rna", 
            slot="data", 
            max.cutoff = "q90", order=TRUE) + theme_void() +
  labs(title = "") + 
  theme(legend.position = "none")
dev.off()


########## 15/19
png("poster/region_tfidf_15_19.png",
    height = 7,
    width = 9,
    units = "cm",
    res=100)
FeaturePlot(object1C, features="chr11-119583953-119585111", 
            reduction = "umap_rna", 
            slot="data", 
            max.cutoff = "q95", order=TRUE) + theme_void() +
  labs(color="TF-IDF", title = "") + 
  theme(legend.title = element_text(size=18), legend.text = element_text(size=18))
dev.off()


png("poster/region_tfidf_15_19_marker.png",
    height = 7,
    width = 9,
    units = "cm",
    res=100)
FeaturePlot(object1C, features=GRangesToString(cisTopic1C@region.ranges[topic.15.topic.19.top.idx@to[28]]), 
            reduction = "umap_rna", 
            slot="data", 
            max.cutoff = "q95", order=TRUE) + theme_void() +
  labs(color="TF-IDF", title = "") + 
  theme(legend.title = element_text(size=18), legend.text = element_text(size=18))
dev.off()


#31/15
png("poster/region_tfidf_31_15.png",
    height = 7,
    width = 9,
    units = "cm",
    res=100)
FeaturePlot(object1C, features="chr15-22238356-22239653", 
            reduction = "umap_rna", 
            slot="data", 
            max.cutoff = "q95", order=TRUE) + theme_void() +
  labs(color="TF-IDF", title = "") + 
  theme(legend.title = element_text(size=18), legend.text = element_text(size=18))
dev.off()

png("poster/region_tfidf_31_15_marker.png",
    height = 7,
    width = 9,
    units = "cm",
    res=100)
FeaturePlot(object1C, features=GRangesToString(cisTopic1C@region.ranges[topic.31.topic.15.top.idx@to[25]]), 
            reduction = "umap_rna", 
            slot="data", 
            max.cutoff = "q95", order=TRUE) + theme_void() +
  labs(color="TF-IDF", title = "") + 
  theme(legend.title = element_text(size=18), legend.text = element_text(size=18))
dev.off()
