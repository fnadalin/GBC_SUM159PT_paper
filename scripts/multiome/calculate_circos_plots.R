library(circlize)
library(readxl)
###################### Circos plot Topics ###################################

topic.2.topic.16 <- read.csv2("hpc/cisTopic_all/1C_ARC_seurat_cisTopic_1E_ARC_seurat_cisTopic_2_16_IDR_0.05_with_markers.csv")
topic.15.topic.19 <- read.csv2("hpc/cisTopic_all/1C_ARC_seurat_cisTopic_1E_ARC_seurat_cisTopic_15_19_IDR_0.05_with_markers.csv")
topic.27.topic.14 <- read.csv2("hpc/cisTopic_all/1C_ARC_seurat_cisTopic_1E_ARC_seurat_cisTopic_27_14_IDR_0.05_with_markers.csv")
topic.31.topic.15 <- read.csv2("hpc/cisTopic_all/1C_ARC_seurat_cisTopic_1E_ARC_seurat_cisTopic_31_15_IDR_0.05_with_markers.csv")

cnv.1 <- read_excel("CNVkit_D15A_19168.xlsx")
cnv.2 <- read_excel("CNVkit_D15B_19169.xlsx")
cnv.3 <- read_excel("CNVkit_D15C_19170.xlsx")

######################### 15/19 #############################
cytoband = read.cytoband(species = "hg38")
df = cytoband$df
chromosome = cytoband$chromosome
chr.len = cytoband$chr.len

df_zoom = df[df[[1]] %in% c("chr7","chr11"), ]
df_zoom[[1]] = paste0(df_zoom[[1]], "_zoom")
df = rbind(df, df_zoom)

bed_cnv_1 = cnv.1[order(cnv.1$log2, decreasing = T)[1:20],1:3]
bed_cnv_1_zoom = bed_cnv_1[bed_cnv_1[[1]] %in% c("chr7","chr11"), ]
bed_cnv_1_zoom[[1]] = paste0(bed_cnv_1_zoom[[1]], "_zoom")
bed_cnv_1 = rbind(bed_cnv_1, bed_cnv_1_zoom)

bed_cnv_2 = cnv.2[order(cnv.2$log2, decreasing = T)[1:20],1:3]
bed_cnv_2_zoom = bed_cnv_2[bed_cnv_2[[1]] %in% c("chr7","chr11"), ]
bed_cnv_2_zoom[[1]] = paste0(bed_cnv_2_zoom[[1]], "_zoom")
bed_cnv_2 = rbind(bed_cnv_2, bed_cnv_2_zoom)

bed_cnv_3 = cnv.3[order(cnv.3$log2, decreasing = T)[1:20],1:3]
bed_cnv_3_zoom = bed_cnv_3[bed_cnv_3[[1]] %in% c("chr7","chr11"), ]
bed_cnv_3_zoom[[1]] = paste0(bed_cnv_3_zoom[[1]], "_zoom")
bed_cnv_3 = rbind(bed_cnv_3, bed_cnv_3_zoom)

bed_col = rep("white", 26)

bed_col[c(7, 25)] <- "#0000FF10"
bed_col[c(11, 26)] <- "#FF000010"

bed = topic.15.topic.19
bed_zoom = bed[bed[[1]] %in% c("chr7","chr11"), ]
bed_zoom[[1]] = paste0(bed_zoom[[1]], "_zoom")
bed = rbind(bed, bed_zoom)

png("circos_plot_topic_15_19_CNV.png",
    height=20,
    width = 19,
    units = "cm",
    res=600)
circos.par(start.degree = 90)
circos.initializeWithIdeogram(df, sort.chr = FALSE, sector.width = c(chr.len/sum(chr.len), 0.5, 0.5))
circos.genomicTrackPlotRegion(bed, panel.fun = function(region, value, ...) {
  circos.genomicPoints(region, value, pch = 20, cex = 0.7)
}, track.height=0.09, bg.col=add_transparency(bed_col, 0.8))

circos.yaxis(
  sector.index = "chr1",
  track.index = 3,
  labels.font = par("font"),
  labels.cex = 0.3)
circos.text(x=1, y=800, labels = "IDR",sector.index = "chr1",
            track.index = 3, facing = "reverse.clockwise", niceFacing = TRUE, adj = c(0.4, -0.3))


circos.genomicTrackPlotRegion(bed_cnv_1, stack = TRUE, panel.fun = function(region, value, ...) {
  circos.genomicRect(region, value, col = "red", border = NA)
}, bg.border = NA, track.height=0.05, track.margin = mm_h(h = c(0.001,0.001)))

circos.genomicTrackPlotRegion(bed_cnv_2, stack = TRUE, panel.fun = function(region, value, ...) {
  circos.genomicRect(region, value, col = "green", border = NA)
}, bg.border = NA, track.height=0.05, track.margin = mm_h(h = c(0.001,0.001)))

circos.genomicTrackPlotRegion(bed_cnv_3, stack = TRUE, panel.fun = function(region, value, ...) {
  circos.genomicRect(region, value, col = "blue", border = NA)
}, bg.border = NA, track.height=0.05, track.margin = mm_h(h = c(0.001,0.001)))

circos.link("chr7", get.cell.meta.data("cell.xlim", sector.index = "chr7"),
            "chr7_zoom", get.cell.meta.data("cell.xlim", sector.index = "chr7_zoom"), 
            col = "#0000FF10", border = NA)
circos.link("chr11", get.cell.meta.data("cell.xlim", sector.index = "chr11"),
            "chr11_zoom", get.cell.meta.data("cell.xlim", sector.index = "chr11_zoom"), 
            col = "#FF000010", border = NA)

title("Topic 15/19")
circos.clear()
dev.off()
############################ 2/16 #################################

bed = topic.2.topic.16
bed_zoom = bed[bed[[1]] %in% c("chr7","chr11"), ]
bed_zoom[[1]] = paste0(bed_zoom[[1]], "_zoom")
bed = rbind(bed, bed_zoom)

png("circos_plot_topic_2_16_CNV.png",
    height=20,
    width = 19,
    units = "cm",
    res=600)
circos.par(start.degree = 90)
circos.initializeWithIdeogram(df, sort.chr = FALSE, sector.width = c(chr.len/sum(chr.len), 0.5, 0.5))
circos.genomicTrackPlotRegion(bed, panel.fun = function(region, value, ...) {
  circos.genomicPoints(region, value, pch = 20, cex = 0.7)
}, track.height=0.09, bg.col=add_transparency(bed_col, 0.8))

circos.yaxis(
  sector.index = "chr1",
  track.index = 3,
  labels.font = par("font"),
  labels.cex = 0.3)
circos.text(x=1, y=800, labels = "IDR",sector.index = "chr1",
            track.index = 3, facing = "reverse.clockwise", niceFacing = TRUE, adj = c(0.4, -0.3))


circos.genomicTrackPlotRegion(bed_cnv_1, stack = TRUE, panel.fun = function(region, value, ...) {
  circos.genomicRect(region, value, col = "red", border = NA)
}, bg.border = NA, track.height=0.05, track.margin = mm_h(h = c(0.001,0.001)))

circos.genomicTrackPlotRegion(bed_cnv_2, stack = TRUE, panel.fun = function(region, value, ...) {
  circos.genomicRect(region, value, col = "green", border = NA)
}, bg.border = NA, track.height=0.05, track.margin = mm_h(h = c(0.001,0.001)))

circos.genomicTrackPlotRegion(bed_cnv_3, stack = TRUE, panel.fun = function(region, value, ...) {
  circos.genomicRect(region, value, col = "blue", border = NA)
}, bg.border = NA, track.height=0.05, track.margin = mm_h(h = c(0.001,0.001)))

circos.link("chr7", get.cell.meta.data("cell.xlim", sector.index = "chr7"),
            "chr7_zoom", get.cell.meta.data("cell.xlim", sector.index = "chr7_zoom"), 
            col = "#0000FF10", border = NA)
circos.link("chr11", get.cell.meta.data("cell.xlim", sector.index = "chr11"),
            "chr11_zoom", get.cell.meta.data("cell.xlim", sector.index = "chr11_zoom"), 
            col = "#FF000010", border = NA)

title("Topic 2/16")
circos.clear()
dev.off()

############################ 27/14 ###############################
bed = topic.27.topic.14
bed_zoom = bed[bed[[1]] %in% c("chr7","chr11"), ]
bed_zoom[[1]] = paste0(bed_zoom[[1]], "_zoom")
bed = rbind(bed, bed_zoom)

png("circos_plot_topic_27_14_CNV.png",
    height=20,
    width = 19,
    units = "cm",
    res=600)
circos.par(start.degree = 90)
circos.initializeWithIdeogram(df, sort.chr = FALSE, sector.width = c(chr.len/sum(chr.len), 0.5, 0.5))
circos.genomicTrackPlotRegion(bed, panel.fun = function(region, value, ...) {
  circos.genomicPoints(region, value, pch = 20, cex = 0.7)
}, track.height=0.09, bg.col=add_transparency(bed_col, 0.8))

circos.yaxis(
  sector.index = "chr1",
  track.index = 3,
  labels.font = par("font"),
  labels.cex = 0.3)
circos.text(x=1, y=800, labels = "IDR",sector.index = "chr1",
            track.index = 3, facing = "reverse.clockwise", niceFacing = TRUE, adj = c(0.4, -0.3))


circos.genomicTrackPlotRegion(bed_cnv_1, stack = TRUE, panel.fun = function(region, value, ...) {
  circos.genomicRect(region, value, col = "red", border = NA)
}, bg.border = NA, track.height=0.05, track.margin = mm_h(h = c(0.001,0.001)))

circos.genomicTrackPlotRegion(bed_cnv_2, stack = TRUE, panel.fun = function(region, value, ...) {
  circos.genomicRect(region, value, col = "green", border = NA)
}, bg.border = NA, track.height=0.05, track.margin = mm_h(h = c(0.001,0.001)))

circos.genomicTrackPlotRegion(bed_cnv_3, stack = TRUE, panel.fun = function(region, value, ...) {
  circos.genomicRect(region, value, col = "blue", border = NA)
}, bg.border = NA, track.height=0.05, track.margin = mm_h(h = c(0.001,0.001)))

circos.link("chr7", get.cell.meta.data("cell.xlim", sector.index = "chr7"),
            "chr7_zoom", get.cell.meta.data("cell.xlim", sector.index = "chr7_zoom"), 
            col = "#0000FF10", border = NA)
circos.link("chr11", get.cell.meta.data("cell.xlim", sector.index = "chr11"),
            "chr11_zoom", get.cell.meta.data("cell.xlim", sector.index = "chr11_zoom"), 
            col = "#FF000010", border = NA)

title("Topic 27/14")
circos.clear()
dev.off()

############################ 31/15 ###############################
bed = topic.31.topic.15
bed_zoom = bed[bed[[1]] %in% c("chr7","chr11"), ]
bed_zoom[[1]] = paste0(bed_zoom[[1]], "_zoom")
bed = rbind(bed, bed_zoom)

png("circos_plot_topic_31_15_CNV.png",
    height=20,
    width = 19,
    units = "cm",
    res=600)
circos.par(start.degree = 90)
circos.initializeWithIdeogram(df, sort.chr = FALSE, sector.width = c(chr.len/sum(chr.len), 0.5, 0.5))
circos.genomicTrackPlotRegion(bed, panel.fun = function(region, value, ...) {
  circos.genomicPoints(region, value, pch = 20, cex = 0.7)
}, track.height=0.09, bg.col=add_transparency(bed_col, 0.8))

circos.yaxis(
  sector.index = "chr1",
  track.index = 3,
  labels.font = par("font"),
  labels.cex = 0.3)
circos.text(x=1, y=800, labels = "IDR",sector.index = "chr1",
            track.index = 3, facing = "reverse.clockwise", niceFacing = TRUE, adj = c(0.4, -0.3))

circos.genomicTrackPlotRegion(bed_cnv_1, stack = TRUE, panel.fun = function(region, value, ...) {
  circos.genomicRect(region, value, col = "red", border = NA)
}, bg.border = NA, track.height=0.05, track.margin = mm_h(h = c(0.001,0.001)))

circos.genomicTrackPlotRegion(bed_cnv_2, stack = TRUE, panel.fun = function(region, value, ...) {
  circos.genomicRect(region, value, col = "green", border = NA)
}, bg.border = NA, track.height=0.05, track.margin = mm_h(h = c(0.001,0.001)))

circos.genomicTrackPlotRegion(bed_cnv_3, stack = TRUE, panel.fun = function(region, value, ...) {
  circos.genomicRect(region, value, col = "blue", border = NA)
}, bg.border = NA, track.height=0.05, track.margin = mm_h(h = c(0.001,0.001)))

circos.link("chr7", get.cell.meta.data("cell.xlim", sector.index = "chr7"),
            "chr7_zoom", get.cell.meta.data("cell.xlim", sector.index = "chr7_zoom"), 
            col = "#0000FF10", border = NA)
circos.link("chr11", get.cell.meta.data("cell.xlim", sector.index = "chr11"),
            "chr11_zoom", get.cell.meta.data("cell.xlim", sector.index = "chr11_zoom"), 
            col = "#FF000010", border = NA)

title("Topic 31/15")
circos.clear()
dev.off()

