
set.seed(1234)
ALGORITHM <- 3 # SLM algorithm is set as default
RESOLUTION <- 0.5 # clustering resolution is set to 0.5
############################## OPTIONS MENU ###################################

library("optparse")

# specify our desired options in a list
# by default OptionParser will add an help option equivalent to 
# make_option(c("-h", "--help"), action="store_true", default=FALSE, 
#               help="Show this help message and exit")

cat("\n")
option_list <- list(
  make_option("--Seurat", type="character", 
              help="[REQUIRED] Seurat Object path"))


############################# PARSE OPTIONS ###################################

opt <- parse_args(OptionParser(option_list=option_list))

if (is.null(opt$Seurat)){
  write("Option --Seurat is required\nTry --help for help", stderr())
  q()
} else {
  SEURAT <- opt$Seurat
}

############################### EXECUTION #####################################
library(Seurat)
library(Signac)
library(ggplot2)
library(cisTopic)

# load Seurat

object <- readRDS(SEURAT)

plots_dir <- paste(paste(dirname(SEURAT),"plots","",sep='/'))

if (!dir.exists(plots_dir)){
  dir.create(plots_dir, recursive = T)
}

exp_name <- sub("_seurat.Rds","", basename(SEURAT))

# features plots with Topics scores
topics <- colnames(object@reductions$cisTopic@cell.embeddings)
n.topics <- length(object@reductions$cisTopic@misc$cisTopicObject@selected.model$topic_sums)

for (topic in topics) {
  plot.name <- paste(c(exp_name,strsplit(topic, "_")[[1]]), collapse = " ")
  plot_path <- paste0(plots_dir, paste(exp_name,n.topics,"Topics", "featureplot", topic, "umap_rna.png", sep='_'))
  p1 <- FeaturePlot(object, reduction = "umap_rna", features = topic) +
    ggtitle(plot.name) +
    labs(colour="Z-score")
  ggsave(plot = p1, filename = plot_path, units = "cm", width = 15, height = 12, 
         device = "png")
}


saveRDS(object, SEURAT)
q()
