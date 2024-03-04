
png("poster/page5_fez1_umap.png",
    height = 5,
    width = 5.5,
    units = "cm",
    res=100)
DimPlot(object1C, 
        reduction = "umap_rna", 
        cells.highlight = list(S1=names(object1C$pca_rna_0.4_clusters[object1C$pca_rna_0.4_clusters == 6]), S2=names(object1C$pca_rna_0.4_clusters[object1C$pca_rna_0.4_clusters == 4])),
        cols.highlight = c("#00c3ef","#fe7edf")) +theme_void() +
  theme(legend.position = "none")
dev.off()

png("poster/page5_umap.png",
    height = 5,
    width = 5.5,
    units = "cm",
    res=100)
DimPlot(object1C, 
        reduction = "umap_rna", 
        cells.highlight = list(S1=names(object1C$pca_rna_0.4_clusters[object1C$pca_rna_0.4_clusters == 6])),
        cols.highlight = c("#fe7edf")) +theme_void() +
  theme(legend.position = "none")
dev.off()

png("poster/fez1_umap.png",
    height = 5,
    width = 5.5,
    units = "cm",
    res=100)
DimPlot(object1C, 
        reduction = "umap_rna", 
        cells.highlight = list(S1=names(object1C$pca_rna_0.4_clusters[object1C$pca_rna_0.4_clusters == 4])),
        cols.highlight = c("#00c3ef")) +theme_void() +
  theme(legend.position = "none")
dev.off()


png("poster/cl3_umap.png",
    height = 5,
    width = 5.5,
    units = "cm",
    res=100)
DimPlot(object1C, 
        reduction = "umap_rna", 
        cells.highlight = list(S1=names(object1C$pca_rna_0.4_clusters[object1C$pca_rna_0.4_clusters == 3])),
        cols.highlight = c("#00c9a5")) +theme_void() +
  theme(legend.position = "none")
dev.off()

clusters <- releobject1C$pca_rna_0.4_clusters
clusters <- factor(x=object1C$pca_rna_0.4_clusters,
                           levels = c("Cluster 0","Cluster 1","Cluster 2","S2","S3","S1"),
                           ordered = TRUE)

DefaultAssay(object1C) <- "RNA"
Idents(object1C) <- object1C$pca_rna_0.4_clusters
object1C.markers <- FindAllMarkers(object1C, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
object1C.markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC) -> top10

png("poster/markers_heatmap.png",
    height = 16,
    width = 26,
    res=100,
    units="cm")
DoHeatmap(object1C, features = top10$gene) +
  scale_fill_gradient2(low="blue", mid="white", high = "red") +
  guides(label="none", labels="none", groups="none") +
  theme(axis.text.y = element_text(size=16),
        legend.title = element_text(size=18),
        legend.text = element_text(size=18))
dev.off()



