library(Seurat)
library(pROC)
library(tidyverse)
library(ggrepel)
library(ggplot2)
object1C <- readRDS("hpc/cisTopic_all/1C_ARC_seurat.Rds")
object1E <- readRDS("hpc/cisTopic_all/1E_ARC_seurat.Rds")


calc.roc.1C <- function(topic, label) {
  res <- roc(
    object1C$pca_rna_0.4_clusters == label,
    object1C@reductions$cisTopic@cell.embeddings[,topic],
    direction="<"
  )
  return(data.frame(topic=rep(topic, length(res[["specificities"]])), 
                    specificities=res[["specificities"]], 
                    sensitivities=res[["sensitivities"]]))
}

calc.roc.1E <- function(topic, label) {
  res <- roc(
    object1E$pca_rna_0.6_clusters == label,
    object1E@reductions$cisTopic@cell.embeddings[,topic],
    direction="<"
  )
  return(data.frame(topic=rep(topic, length(res[["specificities"]])), 
                    specificities=res[["specificities"]], 
                    sensitivities=res[["sensitivities"]]))
}


calc.auc.1C <- function(topic, label) {
  res <- roc(
    object1C$pca_rna_0.4_clusters == label,
    object1C@reductions$cisTopic@cell.embeddings[,topic],
    direction="<"
  )
  auc.vector <- ifelse(res$auc > 0.8, rep(paste("AUC =",round(res$auc,2)),length(res[["specificities"]])), 
                       rep("", length(res[["specificities"]])))
  return(data.frame(topic=rep(topic, length(res[["specificities"]])), 
                    AUC=auc.vector))
}

calc.auc.1E <- function(topic, label) {
  res <- roc(
    object1E$pca_rna_0.6_clusters == label,
    object1E@reductions$cisTopic@cell.embeddings[,topic],
    direction="<"
  )
  auc.vector <- ifelse(res$auc > 0.8, rep(paste("AUC =",round(res$auc,2)),length(res[["specificities"]])), 
                       rep("", length(res[["specificities"]])))
  return(data.frame(topic=rep(topic, length(res[["specificities"]])), 
                    AUC=auc.vector))
}

Topics.1C.FEZ1 <- mapply(calc.roc.1C, topic=c(2,27,15,31), label="4", SIMPLIFY = FALSE)
Topics.1C.FEZ1 <- Topics.1C.FEZ1 %>% purrr::reduce(full_join)
Topics.1C.FEZ1.auc <- mapply(calc.auc.1C, topic=c(2,27,15,31), label="4", SIMPLIFY = FALSE)
Topics.1C.FEZ1.auc <- Topics.1C.FEZ1.auc %>% purrr::reduce(full_join)
Topics.1C.FEZ1 <-  cbind(Topics.1C.FEZ1, AUC=Topics.1C.FEZ1.auc$AUC)
Topics.1C.FEZ1$replicate <- rep("1C", nrow(Topics.1C.FEZ1))
Topics.1C.FEZ1$cluster <- rep("FEZ1+", nrow(Topics.1C.FEZ1))
Topics.1C.FEZ1$topic <- factor(Topics.1C.FEZ1$topic)
levels(Topics.1C.FEZ1$topic) <- c("Topic 2/16", "Topic 15/19", "Topic 27/14", "Topic 31/15")

Topics.1C.PAGE5 <- mapply(calc.roc.1C, topic=c(2,27,15,31), label="6", SIMPLIFY = FALSE)
Topics.1C.PAGE5 <- Topics.1C.PAGE5 %>% purrr::reduce(full_join)
Topics.1C.PAGE5.auc <- mapply(calc.auc.1C, topic=c(2,27,15,31), label="6", SIMPLIFY = FALSE)
Topics.1C.PAGE5.auc <- Topics.1C.PAGE5.auc %>% purrr::reduce(full_join)
Topics.1C.PAGE5 <-  cbind(Topics.1C.PAGE5, AUC=Topics.1C.PAGE5.auc$AUC)
Topics.1C.PAGE5$replicate <- rep("1C", nrow(Topics.1C.PAGE5))
Topics.1C.PAGE5$cluster <- rep("PAGE5+", nrow(Topics.1C.PAGE5))
Topics.1C.PAGE5$topic <- factor(Topics.1C.PAGE5$topic)
levels(Topics.1C.PAGE5$topic) <- c("Topic 2/16", "Topic 15/19", "Topic 27/14", "Topic 31/15")

Topics.1C.C3 <- mapply(calc.roc.1C, topic=c(2,27,15,31), label="3", SIMPLIFY = FALSE)
Topics.1C.C3 <- Topics.1C.C3 %>% purrr::reduce(full_join)
Topics.1C.C3.auc <- mapply(calc.auc.1C, topic=c(2,27,15,31), label="3", SIMPLIFY = FALSE)
Topics.1C.C3.auc <- Topics.1C.C3.auc %>% purrr::reduce(full_join)
Topics.1C.C3 <- cbind(Topics.1C.C3, AUC=Topics.1C.C3.auc$AUC)
Topics.1C.C3$replicate <- rep("1C", nrow(Topics.1C.C3))
Topics.1C.C3$cluster <- rep("Cluster 3", nrow(Topics.1C.C3))
Topics.1C.C3$topic <- factor(Topics.1C.C3$topic)
levels(Topics.1C.C3$topic) <- c("Topic 2/16", "Topic 15/19", "Topic 27/14", "Topic 31/15")


Topics.1E.FEZ1 <- mapply(calc.roc.1E, topic=c(16,14,19,15), label="4", SIMPLIFY = FALSE)
Topics.1E.FEZ1 <- Topics.1E.FEZ1 %>% purrr::reduce(full_join)
Topics.1E.FEZ1.auc <- mapply(calc.auc.1E, topic=c(16,14,19,15), label="4", SIMPLIFY = FALSE)
Topics.1E.FEZ1.auc <- Topics.1E.FEZ1.auc %>% purrr::reduce(full_join)
Topics.1E.FEZ1 <- cbind(Topics.1E.FEZ1, AUC=Topics.1E.FEZ1.auc$AUC)
Topics.1E.FEZ1$replicate <- rep("1E", nrow(Topics.1E.FEZ1))
Topics.1E.FEZ1$cluster <- rep("FEZ1+", nrow(Topics.1E.FEZ1))
Topics.1E.FEZ1$topic <- factor(Topics.1E.FEZ1$topic)
levels(Topics.1E.FEZ1$topic) <- c("Topic 27/14","Topic 31/15","Topic 2/16","Topic 15/19")

Topics.1E.PAGE5 <- mapply(calc.roc.1E, topic=c(16,14,19,15), label="5", SIMPLIFY = FALSE)
Topics.1E.PAGE5 <- Topics.1E.PAGE5 %>% purrr::reduce(full_join)
Topics.1E.PAGE5.auc <- mapply(calc.auc.1E, topic=c(16,14,19,15), label="5", SIMPLIFY = FALSE)
Topics.1E.PAGE5.auc <- Topics.1E.PAGE5.auc %>% purrr::reduce(full_join)
Topics.1E.PAGE5 <- cbind(Topics.1E.PAGE5, AUC=Topics.1E.PAGE5.auc$AUC)
Topics.1E.PAGE5$replicate <- rep("1E", nrow(Topics.1E.PAGE5))
Topics.1E.PAGE5$cluster <- rep("PAGE5+", nrow(Topics.1E.PAGE5))
Topics.1E.PAGE5$topic <- factor(Topics.1E.PAGE5$topic)
levels(Topics.1E.PAGE5$topic) <- c("Topic 27/14","Topic 31/15","Topic 2/16","Topic 15/19")

Topics.1E.C3 <- mapply(calc.roc.1E, topic=c(16,14,19,15), label="3", SIMPLIFY = FALSE)
Topics.1E.C3 <- Topics.1E.C3 %>% purrr::reduce(full_join)
Topics.1E.C3.auc <- mapply(calc.auc.1E, topic=c(16,14,19,15), label="3", SIMPLIFY = FALSE)
Topics.1E.C3.auc <- Topics.1E.C3.auc %>% purrr::reduce(full_join)
Topics.1E.C3 <- cbind(Topics.1E.C3, AUC=Topics.1E.C3.auc$AUC)
Topics.1E.C3$replicate <- rep("1E", nrow(Topics.1E.C3))
Topics.1E.C3$cluster <- rep("Cluster 3", nrow(Topics.1E.C3))
Topics.1E.C3$topic <- factor(Topics.1E.C3$topic)
levels(Topics.1E.C3$topic) <- c("Topic 27/14","Topic 31/15","Topic 2/16","Topic 15/19")

Total.roc <- rbind(Topics.1C.FEZ1, Topics.1C.PAGE5, Topics.1C.C3, Topics.1E.FEZ1, Topics.1E.PAGE5, Topics.1E.C3)
Total.roc$cluster <- factor(Total.roc$cluster, levels = c("PAGE5+", "FEZ1+", "Cluster 3"), ordered = TRUE)
Total.roc$replicate <- factor(Total.roc$replicate)

write.csv2(Total.roc, "hpc/cisTopic_all/cistopic_clusters_roc_curve_plot_ggplot_df.csv", quote = F, row.names = F)

########################## To reproduce the plots only ########################

Total.roc <- read.csv2("hpc/cisTopic_all/cistopic_clusters_roc_curve_plot_ggplot_df.csv")
Total.roc$cluster <- factor(Total.roc$cluster, levels = c("PAGE5+", "FEZ1+", "Cluster 3"), ordered = TRUE)
Total.roc$topic <- factor(Total.roc$topic, levels = c("Topic 2/16","Topic 31/15","Topic 27/14","Topic 15/19"), ordered = TRUE)
Total.roc$replicate <- factor(Total.roc$replicate)
anno.roc <- unique(Total.roc[,c("topic","AUC","replicate","cluster")])

png("ROC_page5_fez1_s3_cistopics.png",
    height=20,
    width = 29,
    res=800,
    units="cm")
ggplot(data = Total.roc,
       aes(x=1-sensitivities, y=specificities, color=replicate, group=replicate)) +
  geom_line(key_glyph = "rect") +
  labs(x="1-Sensitivity", y="Specificity", color="Replicate") +
  facet_grid(cluster~topic) +
  geom_text(
    data=anno.roc,
    mapping = aes(label=AUC, x=0.8, y=0.1),
    position = position_stack()
  ) +
  theme_light() +
  theme(text = element_text(size=12)) +
  guides(label="none")
dev.off()

png("ROC_page5_fez1_s3_cistopics_by_clusters.png",
    height=20,
    width=40,
    res=800,
    units="cm")
ggplot(data = Total.roc,
       aes(x=1-sensitivities, y=specificities, color=cluster)) +
  geom_line(key_glyph = "rect") +
  scale_color_discrete() +
  labs(x="1-Sensitivity", y="Specificity", color="Subpopulation") +
  facet_grid(replicate~topic) +
  geom_text(
    data=anno.roc,
    mapping = aes(label=AUC, x=0.8, y=0.1),
    position = position_stack()
  ) +
  theme_light() +
  theme(text = element_text(size=12)) +
  guides(label="none")
dev.off()


png("poster/roc_curves.png",
    width = 8,
    height = 25,
    units = "cm",
    res=100)
ggplot(data = subset(Total.roc, replicate == "1C"),
       aes(x=1-sensitivities, y=specificities, color=cluster)) +
  geom_line(key_glyph = "rect") +
  scale_color_manual(values=c("#00c3ef","#fe7edf","#00c9a5"), labels=c("S1","S3","S2"))+
  facet_grid(topic~.) +
  geom_text(
    data=subset(anno.roc, replicate == "1C"),
    mapping = aes(label=AUC, x=0.7, y=0.1),
    position = position_stack(vjust = 0.5)
  ) +
  theme_cowplot() +
  labs(x="1-Sensitivity", y="Specificity", color="Subpopulation") +
  theme(text = element_text(size=18), 
        strip.text = element_text(vjust=1),
        strip.background = element_rect(fill="white"),
        axis.text.x = element_text(angle=90, hjust=0.5, vjust = 0.5),
        legend.position = "bottom",
        legend.direction = "vertical") +
  guides(label="none")
dev.off()

