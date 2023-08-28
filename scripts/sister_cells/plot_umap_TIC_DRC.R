
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
    cat("\nUsage: <indir> <outdir> <meta> <params>\n\n")
    q()
}

INDIR <- args[1]
OUTDIR <- args[2]
META <- args[3]
PARAMS <- args[4]

library("Seurat")

params <- as.matrix(read.table(PARAMS, sep="="))
disps <- as.numeric(unlist(strsplit(params[params[,1]=="disps", 2], split=",")))
nfeats <- as.numeric(unlist(strsplit(params[params[,1]=="nfeats", 2], split=",")))

OBJECT <- file.path(INDIR, "object.Rds")
object <- readRDS(OBJECT)

dr_vec <- NULL
for (feature_method in c("mean.var.plot", "vst")) {
    if (feature_method == "mean.var.plot") {
        for (ymin in disps) {
            case <- paste0("mean.var.plot_disp", ymin)
            dr <- paste("pca", case, sep="_")
            dr_vec <- c(dr_vec, dr)
        }
    } else { 
        for (nfeat in nfeats) {
            case <- paste0("vst_top", nfeat)
            dr <- paste("pca", case, sep="_")
            dr_vec <- c(dr_vec, dr)
        }
    }
}

for (dr in dr_vec) {
    outdir <- file.path(OUTDIR,dr)
    dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
    out <- file.path(outdir, paste0("umap",META,".pdf"))
    umap_dr <- umap_dr <- paste0("umap_", dr)
    pdf(out)
    print(DimPlot(object, reduction = umap_dr, group.by = META, cols=c("gray90","red"), pt.size=1, order="1"))
    dev.off()
}

q()

