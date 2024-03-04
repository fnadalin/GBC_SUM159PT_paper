#' Bubble plot with all clusters for a single replicate
#' 
#' Produces two bubble plots with NO consistency information of the first 30 topics \
#' sorted by reproducibility score and a bubble plot with just the three conserved clusters.
#' 
#' @param matchingtopics the file containing information about topic pairs, \
#' reproducibility scores, and Cluster-topic AUC
#' @param out Output file name.
#' 
#' @return two png files containing the bubble plot for single replicates, one png \
#' file with both replicates and consistency measure.
#' 


############################## OPTIONS MENU ###################################

library("optparse")

# specify our desired options in a list
# by default OptionParser will add an help option equivalent to 
# make_option(c("-h", "--help"), action="store_true", default=FALSE, 
#               help="Show this help message and exit")

cat("\n")
option_list <- list(
  make_option("--matchingtopics", type="character", 
              help="[REQUIRED] File containing the list of matching topics with repr. score and AUC information."),
  make_option("--out", type="character", 
              help="[REQUIRED]  Output file names separated by a comma (3 png files, replicate1, replicate2, both)")
)

############################## PARSE OPTIONS ##################################


# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults
opt <- parse_args(OptionParser(option_list=option_list))

if (is.null(opt$matchingtopics)) {
  write("Option --matchingtopics is required\nTry --help for help", stderr())
  q()
} else {
  MATCHING.TOPICS <- opt$matchingtopics
}

if (is.null(opt$out)) {
  write("Option --out is required\nTry --help for help", stderr())
  q()
} else {
  OUT <- unlist(strsplit(opt$out, ","))
}

################################# EXECUTION ###################################

library(ggplot2)
library(ggnewscale)
library(forcats)
library(scales)
library(cowplot)

matching.topics <- read.csv2(MATCHING.TOPICS)

# cisTopic1.blacklist <- c(1,7,16,19,30,33,38,39,40)
# cisTopic2.blacklist <- c(17,20,21,25,30,36)

top.30.topics <- with(matching.topics, order(repr.score, decreasing = T)[1:30])

## 1C bubble plot top 30 topic AUC for all clusters

matching.topics$module[order(matching.topics$repr.score, decreasing = T)] <- paste("Module", 1:nrow(matching.topics))

df.plot.1C <- reshape2::melt(matching.topics[top.30.topics,], id=c("cisTopic1","cisTopic2","module", "repr.score"), 
                                measure.vars=c("Cluster.0.cisTopic1", "Cluster.1.cisTopic1", "Cluster.2.cisTopic1",
                                               "Cluster.3.cisTopic1", "Cluster.4.cisTopic1", "Cluster.5.cisTopic1", "Cluster.6.cisTopic1"),
                                variable.name="Clusters",
                                value.name = "AUC")

is.first <- rep("", 30)

is.first[!duplicated(matching.topics[top.30.topics,"cisTopic1"]) & !duplicated(matching.topics[top.30.topics,"cisTopic2"])] <- "*"


df.plot.1C$is.first[order(df.plot.1C$repr.score, decreasing = T)] <- as.vector(sapply(is.first, rep, 7))

png(OUT[1],
    height = 20,
    width = 12,
    res= 300,
    units = "cm")
ggplot(data=df.plot.1C) +
  geom_point(aes(x=factor(Clusters), y=fct_reorder(factor(module), repr.score), 
                 color=AUC, size=5)) +
  theme_bw() +
  ggtitle("Modules clusters association (1C)") +
  labs(x="Clusters", y="Modules sorted by reproducibility") +
  scale_color_gradient2(low="blue", mid="grey90", high="red", midpoint = 0.5, name="AUC") +
  scale_x_discrete(labels=c("Cluster 0","Cluster 1","Cluster 2", "S2", "S3", "Cluster 5", "S1")) +
  theme(axis.text.x = element_text(angle=30, hjust = 1)) +
  guides(size="none") +
  new_scale_fill() +
  geom_tile(aes(y=fct_reorder(factor(module), repr.score),fill=repr.score, x=0.01, width=0.5)) +
  scale_fill_gradient2(name="Reproducibility\nscore", low = "yellow", high = "red") +
  geom_text(aes(y=fct_reorder(factor(module), repr.score), x=0.01, label=is.first))
dev.off()

## 1E bubble plot top 30 topic AUC for all clusters

matching.topics$module[order(matching.topics$repr.score, decreasing = T)] <- paste("Module", 1:nrow(matching.topics))

df.plot.1E <- reshape2::melt(matching.topics[top.30.topics,], id=c("cisTopic1","cisTopic2","module", "repr.score"), 
                             measure.vars=c("Cluster.0.cisTopic2", "Cluster.1.cisTopic2", "Cluster.2.cisTopic2",
                                            "Cluster.3.cisTopic2", "Cluster.4.cisTopic2", "Cluster.5.cisTopic2"),
                             variable.name="Clusters",
                             value.name = "AUC")

is.first <- rep("", 30)

is.first[!duplicated(matching.topics[top.30.topics,"cisTopic1"]) & !duplicated(matching.topics[top.30.topics,"cisTopic2"])] <- "*"


df.plot.1E$is.first[order(df.plot.1E$repr.score, decreasing = T)] <- as.vector(sapply(is.first, rep, 6))

png(OUT[2],
    height = 20,
    width = 12,
    res= 300,
    units = "cm")
ggplot(data=df.plot.1E) +
  geom_point(aes(x=factor(Clusters), y=fct_reorder(factor(module), repr.score), 
                 color=AUC, size=5)) +
  theme_bw() +
  ggtitle("Modules clusters association (1E)") +
  labs(x="Clusters", y="Modules sorted by reproducibility") +
  scale_color_gradient2(low="blue", mid="grey90", high="red", midpoint = 0.5, name="AUC") +
  scale_x_discrete(labels=c("Cluster 0","Cluster 1","Cluster 2", "S2", "S3", "S1")) +
  theme(axis.text.x = element_text(angle=30, hjust = 1)) +
  guides(size="none") +
  new_scale_fill() +
  geom_tile(aes(y=fct_reorder(factor(module), repr.score),fill=repr.score, x=0.01, width=0.5)) +
  scale_fill_gradient2(name="Reproducibility\nscore", low = "yellow", high = "red") +
  geom_text(aes(y=fct_reorder(factor(module), repr.score), x=0.01, label=is.first))
dev.off()

####################### Bubble plot with mean AUC & consistency ################


matching.topics$Cluster.3 <- rowMeans(matching.topics[,c("Cluster.3.cisTopic1","Cluster.3.cisTopic2")])
matching.topics$FEZ1 <- rowMeans(matching.topics[,c("Cluster.4.cisTopic1","Cluster.4.cisTopic2")])
matching.topics$PAGE5 <- rowMeans(matching.topics[,c("Cluster.6.cisTopic1","Cluster.5.cisTopic2")])


matching.topics$Cluster.3.delta <- with(matching.topics, abs(Cluster.3.cisTopic1 - Cluster.3.cisTopic2))
matching.topics$FEZ1.delta <- with(matching.topics, abs(Cluster.4.cisTopic1 - Cluster.4.cisTopic2))
matching.topics$PAGE5.delta <- with(matching.topics, abs(Cluster.6.cisTopic1 - Cluster.5.cisTopic2))


# melting a dataframe with mean AUC values and deltas
df.plot <- as.data.frame(cbind(reshape2::melt(matching.topics[top.30.topics,], id=c("cisTopic1","cisTopic2"), 
                                   measure.vars=c("Cluster.3", "FEZ1", "PAGE5"),
                                   variable.name="Clusters",
                                   value.name = "AUC"),
                 delta=reshape2::melt(matching.topics[top.30.topics,], id=c("cisTopic1","cisTopic2"), 
                                measure.vars=c("Cluster.3.delta", "FEZ1.delta", "PAGE5.delta"),
                                variable.name="Clusters",
                                value.name = "delta")[,4]))
# adding topic IDR scores
df.plot <- merge(df.plot, matching.topics[top.30.topics,c("cisTopic1","cisTopic2","repr.score", "module")])


# annotation for top matching pairs
is.first <- rep("", 30)
is.first[!duplicated(matching.topics[top.30.topics,"cisTopic1"]) & !duplicated(matching.topics[top.30.topics,"cisTopic2"])] <- "*"

df.plot$is.first[order(df.plot$repr.score, decreasing = T)] <- as.vector(sapply(is.first, rep, 3))


align_legend <- function(p, hjust = 0.5)
{
  # extract legend
  g <- cowplot::plot_to_gtable(p)
  grobs <- g$grobs
  legend_index <- which(sapply(grobs, function(x) x$name) == "guide-box")
  legend <- grobs[[legend_index]]
  
  # extract guides table
  guides_index <- which(sapply(legend$grobs, function(x) x$name) == "layout")
  
  # there can be multiple guides within one legend box  
  for (gi in guides_index) {
    guides <- legend$grobs[[gi]]
    
    # add extra column for spacing
    # guides$width[5] is the extra spacing from the end of the legend text
    # to the end of the legend title. If we instead distribute it by `hjust:(1-hjust)` on
    # both sides, we get an aligned legend
    spacing <- guides$width[5]
    guides <- gtable::gtable_add_cols(guides, hjust*spacing, 1)
    guides$widths[6] <- (1-hjust)*spacing
    title_index <- guides$layout$name == "title"
    guides$layout$l[title_index] <- 2
    
    # reconstruct guides and write back
    legend$grobs[[gi]] <- guides
  }
  
  # reconstruct legend and write back
  g$grobs[[legend_index]] <- legend
  g
}

png(OUT[3],
    height=24,
    width = 14,
    res=100,
    units = "cm")
p <- ggplot(data=df.plot) +
  geom_point(aes(x=factor(Clusters), y=fct_reorder(factor(module), repr.score), 
                 color=AUC, size=1-delta)) +
  theme_bw() +
  labs(x="Clusters", y="Modules sorted by reproducibility") +
  scale_color_gradient2(low="blue", mid="lightgrey", high="red", midpoint = 0.5, name="AUC") +
  scale_x_discrete(labels=c("S2", "S3", "S1")) +
  scale_size_area(trans="logit") +
  labs(size=stringr::str_wrap("logit(Consistency)", width = 1)) +
  guides(color = guide_colorbar(title.position = "top", 
                             title.hjust = 1,
                             label.position = "left",
                             ),
         size = guide_legend(title.position = "top", 
                             title.hjust = -2,
                             label.position = "left")) +
  theme(axis.text.x = element_text(size=18, angle = 90, hjust=1, vjust = 0.5),
        legend.position = "left",
        legend.box.just = "right",
        legend.title = element_text(size=18),
        legend.text = element_text(size=18),
        legend.title.align = 1,
        axis.title.x = element_text(size=18),
        axis.text.y = element_text(size=16),
        axis.title.y = element_text(size=18)) +
  theme( # remove the vertical grid lines
    panel.grid.major.x = element_blank() ,
    # explicitly set the horizontal lines (or they will disappear too)
    panel.grid.major.y = element_line(size=.05, color="black")
  ) +
  new_scale_fill() +
  geom_tile(aes(y=fct_reorder(factor(module), repr.score),fill=repr.score, x=0.01, width=0.5)) +
  scale_fill_gradient2(name="Reproducibility\nscore", low = "yellow", high = "red") +
  geom_text(aes(y=fct_reorder(factor(module), repr.score), x=0.01, label=is.first))

ggdraw(align_legend(p, hjust=1))
dev.off()




