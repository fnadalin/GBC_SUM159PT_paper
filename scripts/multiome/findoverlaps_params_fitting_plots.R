library(cisTopic)
library(tidyr)
library(GenomicRanges)

num.hits <- c()
num.duplicates <- c()
unique.occurrences <- c()

num.hits.color <- "#6DC0DE"
num.duplicates.color <- "#FFCC96"

#object1 <- readRDS("~/Documents/multiome/1E_ARC/1E_ARC_seurat.Rds")

seurat1.name <- sub("_seurat.Rds","",basename(SEURAT1))
cisTopicObject1 <- readRDS("Documents/multiome/hpc/cisTopic/1C_ARC_seurat_cisTopicObject.Rds")
cisTopicObject2 <- readRDS("~/Documents/multiome/hpc/cisTopic/1E_ARC_seurat_cisTopicObject.Rds")

cisTopicObject1 <- cisTopicObject1C
cisTopicObject2 <- cisTopicObject1E

hits.labeller <- as_labeller(c(
  `num.hits`="Total number of hits",
  `num.duplicates`="Number of duplicate hits")
)

for (distance in seq(0,1000, 10)) {
  hits <- findOverlaps(cisTopicObject1@region.ranges, cisTopicObject2@region.ranges, maxgap = distance)
  duplicates <- length(hits) - min(length(unique(queryHits(hits))),
                                   length(unique(subjectHits(hits))))
  num.hits <- c(num.hits, length(hits))
  num.duplicates <- c(num.duplicates, duplicates)}

dfa <- data.frame(num.hits=num.hits, num.duplicates=num.duplicates, num.unique=num.hits-num.duplicates, distance=seq(0,1000,10)) %>%
  reshape2::melt(select=c("num.hits", "num.duplicates", "num.unique"), id.var="distance")
num.hits <- c()
num.duplicates <- c()                 

p1 <- ggplot(data=dfa, aes(x=distance,y=value)) + 
  geom_line() +
  facet_grid(rows=vars(variable), scales = "free_y", labeller=hits.labeller) +
  theme_bw() +
  labs(x = "Maximum Gap allowed (in bp)", title = "Overlapping regions between 1C and 1E")

num.hits.1e <- c()
num.duplicates.1e <- c()

for (distance in seq(0,1000, 10)) {
  hits <- findOverlaps(cisTopicObject@region.ranges, cisTopicObject@region.ranges, maxgap = distance)
  duplicates <- length(hits) - min(length(unique(queryHits(hits))),
                                   length(unique(subjectHits(hits))))
  num.hits.1e <- c(num.hits.1e, length(hits))
  num.duplicates.1e <- c(num.duplicates.1e, sum(duplicates))}

df.1e <- data.frame(num.hits=num.hits.1e, 
                    num.duplicates=num.duplicates.1e, 
                    distance=seq(0,1000,10)) %>%
  reshape2::melt(select=c("num.hits", "num.duplicates"), id.var="distance")

p2 <- ggplot(data=df.1e, aes(x=distance, y=value)) +
  geom_line() +
  facet_grid(vars(variable), scales = "free_y", 
             labeller = hits.labeller) +
  theme_bw() +
  labs(x = "Maximum Gap allowed (in bp)", title = "Overlapping regions between 1E and 1E")

p1 + p2

df2 <- data.frame(num.hits=num.hits,
                  unique.hits=num.hits-num.duplicates, 
                  distance=seq(0,1000,10)) %>% 
  reshape2::melt(select=c("num.hits", "unique.hits"), 
                 id.var="distance")

p3 <-  ggplot(data=df2, aes(x=distance, y=value)) +
  geom_line(aes(color=variable), size=1.5) +
  theme_bw() + 
  labs(title="Overlapping regions between 1C and 1E", 
       x="Maximum Gap allowed (in bp)", 
       y="Number of hits") +
  theme(legend.title = element_blank(),
        legend.position = c(0.15, 0.85)) +
  scale_color_hue(labels = c("Total number of hits", "Number of unique hits", "Number of non-duplicate hits"))

(p1 + p2)/p3

############################## INCREASING OVERLAP #############################

for (overlap in seq(1000,0,-10)) {
  hits <- findOverlaps(cisTopicObject1@region.ranges, cisTopicObject2@region.ranges, minoverlap = overlap)
  duplicates <- length(hits) - min(length(unique(queryHits(hits))),
                                   length(unique(subjectHits(hits))))
  num.hits <- c(num.hits, length(hits))
  num.duplicates <- c(num.duplicates, duplicates)
  }

df1 <- data.frame(num.hits=num.hits, num.duplicates=num.duplicates, num.unique=num.hits-num.duplicates, distance=-seq(1000,0,-10)) %>%
  reshape2::melt(select=c("num.hits", "num.duplicates", "num.unique"), id.var="distance")

num.hits <- c()
num.duplicates <- c()

p4 <- ggplot(data=df1, aes(x=overlap,y=value)) + 
  geom_line() +
  facet_grid(rows=vars(variable), scales = "free_y", labeller=hits.labeller) +
  theme_bw() +
  scale_x_reverse() +
  labs(x = "Minimum overlap required (in bp)", title = "Overlapping regions between 1C and 1E")
p4 + p1

# dataframe containing all distances between peaks between replicates at varying params between -1000 and 1000
df1[, distance == 0] <- NULL
df.distances <- rbind(dfa, df1[-c(101,202,303),])
df.distance.cast <- reshape2::dcast(data=df.distances, 
                                    distance ~ variable, 
                                    value.var = "value")
df.distance.cast$dup.unique.ratio <- df.distance.cast$num.duplicates / df.distance.cast$num.unique
df.distance.cast$dup.tot.ratio <- df.distance.cast$num.duplicates / df.distance.cast$num.hits
df.distance.cast$unique.tot.ratio <- df.distance.cast$num.unique / df.distance.cast$num.hits
df.distance.cast$unique.dup.ratio <- df.distance.cast$num.unique / df.distance.cast$num.duplicates

df.distances <- reshape2::melt(data=df.distance.cast, 
                               select=c("num.hits", 
                                       "num.duplicates", 
                                       "num.unique", 
                                       "dup.unique.ratio", 
                                       "dup.tot.ratio",
                                       "unique.tot.ratio",
                                       "unique.dup.ratio"),
                               id.var="distance"
                               )

write.csv2(df.distances, 
           file="1C_1E_overlap_parameter_testing.csv",
           quote=F,row.names = F)

p5 <- ggplot(data=subset(df.distances, 
                         !(variable=="dup.tot.ratio" | variable=="dup.unique.ratio")), aes(x=distance, y=value)) +
  geom_line(aes(color=variable), size=1.2) +
  theme_bw() +
  labs(x = "Minoverlap/Maxgap (in bp)", 
       y = "Number of hits",
       title = "Overlapping regions between 1C and 1E") +
  theme(legend.title = element_blank(),
        legend.position = c(0.15, 0.85)) +
  scale_color_hue(labels = c("Total number of hits", 
                             "Number of duplicate hits", 
                             "Number of unique hits"))

p6 <- ggplot(data=subset(df.distances, 
                         variable=="dup.tot.ratio" | variable=="dup.unique.ratio"), aes(x=distance, y=value)) +
  geom_line(aes(color=variable), size=1.2) +
  theme_bw() +
  labs(x = "Minoverlap/Maxgap (in bp)", 
       y = "Number of hits",
       title = "Overlapping regions between 1C and 1E") +
  theme(legend.title = element_blank(),
        legend.position = c(0.15, 0.85)) +
  scale_color_hue(labels = c("Duplicates / Total hits",
                             "Duplicates / Unique hits"))

p7 <- ggplot(data=subset(df.distances, 
                         variable=="unique.tot.ratio" | variable=="unique.dup.ratio"), aes(x=distance, y=value)) +
  geom_line(aes(color=variable), size=1.2) +
  theme_bw() +
  labs(x = "Minoverlap/Maxgap (in bp)", 
       y = "Number of hits",
       title = "Overlapping regions between 1C and 1E") +
  theme(legend.title = element_blank(),
        legend.position = c(0.15, 0.85)) +
  scale_color_hue(labels = c("Unique / Total hits",
                             "UNique / Duplicates hits"))
########################### PEAKS WIDTH ######################################

peaks.width <- data.frame(Experiment=c(rep("1E", length(cisTopicObject1@region.ranges)),
                                       rep("1C", length(cisTopicObject2@region.ranges))),
                          Width=c(cisTopicObject1@region.ranges@ranges@width,
                                  cisTopicObject2@region.ranges@ranges@width))
ggplot(data=peaks.width, aes(x=Experiment)) +
  geom_boxplot(aes(y=Width), outlier.shape = NA) +
  coord_cartesian(ylim = quantile(peaks.width$Width, c(0, 0.92))) +
  theme_bw() +
  ggtitle("Peaks width")













