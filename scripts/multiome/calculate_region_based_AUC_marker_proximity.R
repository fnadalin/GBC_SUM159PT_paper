library(cisTopic)
library(Signac)
library(pROC)
library(ggplot2)
library(tidyverse)

# Find regions proximal to markers within 50kb
cisTopic1C <- readRDS("hpc/cisTopic_all/1C_ARC_seurat_cisTopic.Rds")
cisTopic1E <- readRDS("hpc/cisTopic_all/1E_ARC_seurat_cisTopic.Rds")

cisTopic1C <- selectModel(cisTopic1C,type = "maximum")
cisTopic1E <- selectModel(cisTopic1E,type = "maximum")

signatures <- read.csv2("signatures.csv")

# for each cluster, we retrieve genomic coordinates of all markers
PAGE5.markers <- signatures[grep("PAGE5+", signatures$cluster),]
FEZ1.markers <- signatures[grep("FEZ1+", signatures$cluster),]
Cluster3.markers <- signatures[grep("Cluster3", signatures$cluster),]

total.regions.1C <- cisTopic1C@region.ranges
total.regions.1E <- cisTopic1E@region.ranges

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

regions.1E.markers.overlap.PAGE5 <- findOverlaps(total.regions.1E, 
                                                 StringToGRanges(PAGE5.markers$coords, sep=c(":","-")),
                                                 maxgap = 50000)

regions.1E.markers.overlap.FEZ1 <- findOverlaps(total.regions.1E, 
                                                StringToGRanges(FEZ1.markers$coords, sep=c(":","-")),
                                                maxgap = 50000)

regions.1E.markers.overlap.Cluster3 <- findOverlaps(total.regions.1E, 
                                                    StringToGRanges(Cluster3.markers$coords, sep=c(":","-")),
                                                    maxgap = 50000)

total.regions.1E@elementMetadata$PAGE5 <- rep(FALSE, length(total.regions.1E))
total.regions.1E@elementMetadata$PAGE5[regions.1E.markers.overlap.PAGE5@from] <- TRUE

total.regions.1E@elementMetadata$FEZ1 <- rep(FALSE, length(total.regions.1E))
total.regions.1E@elementMetadata$FEZ1[regions.1E.markers.overlap.FEZ1@from] <- TRUE

total.regions.1E@elementMetadata$Cluster3 <- rep(FALSE, length(total.regions.1E))
total.regions.1E@elementMetadata$Cluster3[regions.1E.markers.overlap.Cluster3@from] <- TRUE

# calculate AUC between region-topic probability and marker proximality
cisTopic1C.region.emb <- modelMatSelection(cisTopic1C, target = "region", method="Probability", all.regions = TRUE)
cisTopic1E.region.emb <- modelMatSelection(cisTopic1E, target = "region", method="Probability", all.regions = TRUE)

# return results and hope for the best

topics.names <- paste("Topic", 1:40, sep = " ")

AUC.1C.regions.PAGE5 <- sapply(1:40, function(x) auc(total.regions.1C$PAGE5, cisTopic1C.region.emb[x,], direction="<"))
AUC.1C.regions.FEZ1 <- sapply(1:40, function(x) auc(total.regions.1C$FEZ1, cisTopic1C.region.emb[x,], direction="<"))
AUC.1C.regions.Cluster3 <- sapply(1:40, function(x) auc(total.regions.1C$Cluster3, cisTopic1C.region.emb[x,], direction="<"))

names(AUC.1C.regions.PAGE5) <- topics.names
names(AUC.1C.regions.FEZ1) <- topics.names
names(AUC.1C.regions.Cluster3) <- topics.names

AUC.1E.regions.PAGE5 <- sapply(1:40, function(x) auc(total.regions.1E$PAGE5, cisTopic1E.region.emb[x,], direction="<"))
AUC.1E.regions.FEZ1 <- sapply(1:40, function(x) auc(total.regions.1E$FEZ1, cisTopic1E.region.emb[x,], direction="<"))
AUC.1E.regions.Cluster3 <- sapply(1:40, function(x) auc(total.regions.1E$Cluster3, cisTopic1E.region.emb[x,], direction="<"))

names(AUC.1E.regions.PAGE5) <- topics.names
names(AUC.1E.regions.FEZ1) <- topics.names
names(AUC.1E.regions.Cluster3) <- topics.names

Topics.AUC.ranks <- as.data.frame(cbind(Topics=topics.names,
                          PAGE5_1C=rank(-AUC.1C.regions.PAGE5),
                          PAGE5_1E=rank(-AUC.1E.regions.PAGE5),
                          FEZ1_1C=rank(-AUC.1C.regions.FEZ1),
                          FEZ1_1E=rank(-AUC.1E.regions.FEZ1),
                          Cluster3_1C=rank(-AUC.1C.regions.Cluster3),
                          Cluster3_1E=rank(-AUC.1E.regions.Cluster3),
                          PAGE5_1C_AUC=AUC.1C.regions.PAGE5,
                          PAGE5_1E_AUC=AUC.1E.regions.PAGE5,
                          FEZ1_1C_AUC=AUC.1C.regions.FEZ1,
                          FEZ1_1E_AUC=AUC.1E.regions.FEZ1,
                          Cluster3_1C_AUC=AUC.1C.regions.Cluster3,
                          Cluster3_1E_AUC=AUC.1E.regions.Cluster3))

Topics.AUC.ranks <- transform(Topics.AUC.ranks,
                              PAGE5_1C=as.numeric(PAGE5_1C),
                              PAGE5_1E=as.numeric(PAGE5_1E),
                              FEZ1_1C=as.numeric(FEZ1_1C),
                              FEZ1_1E=as.numeric(FEZ1_1E),
                              Cluster3_1C=as.numeric(Cluster3_1C),
                              Cluster3_1E=as.numeric(Cluster3_1E),
                              PAGE5_1C_AUC=as.numeric(PAGE5_1C_AUC),
                              PAGE5_1E_AUC=as.numeric(PAGE5_1E_AUC),
                              FEZ1_1C_AUC=as.numeric(FEZ1_1C_AUC),
                              FEZ1_1E_AUC=as.numeric(FEZ1_1E_AUC),
                              Cluster3_1C_AUC=as.numeric(Cluster3_1C_AUC),
                              Cluster3_1E_AUC=as.numeric(Cluster3_1E_AUC))

png(
  "region_based_AUC_PAGE5_bubble_plot.png",
  height = 20,
  width = 10,
  units = "cm",
  res=300
)
ggplot(data=arrange(Topics.AUC.ranks, PAGE5_1C)) +
  lims(x=c(1, 2)) +
  geom_text(aes(x=1, PAGE5_1C, label=Topics), hjust="left", size=5, 
            position = position_nudge(x = 0.05)) +
  geom_text(aes(x=2, PAGE5_1E, label=Topics), hjust="right", size=5, 
            position = position_nudge(x = -0.05)) +
  geom_point(aes(x=1, PAGE5_1C, color=PAGE5_1C_AUC), size=5) +
  geom_point(aes(x=2, PAGE5_1E, color=PAGE5_1E_AUC), size=5) +
  scale_color_gradient2(low="blue", mid="grey", high="red", midpoint = 0.5, name="AUC PAGE5+") +
  scale_y_reverse() +
  annotate('text', x=1, y=-1, label='Multiome 1C', hjust='left', size=5, fontface=2) +
  annotate('text', x=2, y=-1, label='Multiome 1E', hjust='right', size=5, fontface=2) +
  theme_void()
dev.off()

png(
  "region_based_AUC_FEZ1_bubble_plot.png",
  height = 20,
  width = 10,
  units = "cm",
  res=300
)
ggplot(data=arrange(Topics.AUC.ranks, FEZ1_1C)) +
  lims(x=c(1, 2)) +
  geom_text(aes(x=1, FEZ1_1C, label=Topics), hjust="left", size=5, 
            position = position_nudge(x = 0.05)) +
  geom_text(aes(x=2, FEZ1_1E, label=Topics), hjust="right", size=5, 
            position = position_nudge(x = -0.05)) +
  geom_point(aes(x=1, FEZ1_1C, color=FEZ1_1C_AUC), size=5) +
  geom_point(aes(x=2, FEZ1_1E, color=FEZ1_1E_AUC), size=5) +
  scale_color_gradient2(low="blue", mid="grey", high="red", midpoint = 0.5, name="AUC FEZ1+") +
  scale_y_reverse() +
  annotate('text', x=1, y=-1, label='Multiome 1C', hjust='left', size=5, fontface=2) +
  annotate('text', x=2, y=-1, label='Multiome 1E', hjust='right', size=5, fontface=2) +
  theme_void()
dev.off()

png(
  "region_based_AUC_cluster3_bubble_plot.png",
  height = 20,
  width = 10,
  units = "cm",
  res=300
)
ggplot(data=arrange(Topics.AUC.ranks, Cluster3_1C)) +
  lims(x=c(1, 2)) +
  geom_text(aes(x=1, Cluster3_1C, label=Topics), hjust="left", size=5, 
            position = position_nudge(x = 0.05)) +
  geom_text(aes(x=2, Cluster3_1E, label=Topics), hjust="right", size=5, 
            position = position_nudge(x = -0.05)) +
  geom_point(aes(x=1, Cluster3_1C, color=Cluster3_1C_AUC), size=5) +
  geom_point(aes(x=2, Cluster3_1E, color=Cluster3_1E_AUC), size=5) +
  scale_color_gradient2(low="blue", mid="grey", high="red", midpoint = 0.5, name="AUC FEZ1+") +
  scale_y_reverse() +
  annotate('text', x=1, y=-1, label='Multiome 1C', hjust='left', size=5, fontface=2) +
  annotate('text', x=2, y=-1, label='Multiome 1E', hjust='right', size=5, fontface=2) +
  theme_void()
dev.off()








