

if [[ $# -lt "3" ]]
then
    echo ""
    echo "Usage: bash $0 <script_dir> <OUT_DIR> <param_dir>"
    echo ""
    exit
fi

SCRIPT_DIR="$1"
OUT_DIR="$2"
PARAM_DIR="$3"

R_SCRIPT_HVG_PCA="${SCRIPT_DIR}/Seurat_hvg_pca.R"
R_SCRIPT_CLUSTERING="${SCRIPT_DIR}/Seurat_clustering.R"
R_SCRIPT_SILHOUETTE="${SCRIPT_DIR}/Seurat_silhouette.R"
R_SCRIPT_CLUSTER_COMPOSITION="${SCRIPT_DIR}/Seurat_cluster_composition.R"
R_SCRIPT_PLOT="${SCRIPT_DIR}/Seurat_plot.R"

HPC_DIR="${OUT_DIR}/hvg_pca_clust"
SILH_DIR="${OUT_DIR}/silhouette"
COMP_DIR="${OUT_DIR}/cluster_composition"
PLOT_DIR="${OUT_DIR}/plots"

IN_OBJ="${OUT_DIR}/object.Rds"
OUT_OBJ="${HPC_DIR}/object.Rds"

PARAMS_PCA="${PARAM_DIR}/param_hvg_pca.txt"
PARAMS_CL="${PARAM_DIR}/param_clust.txt"

# run feature selection and PCA
# do not run Jackstraw test here!
echo "run feature selection and PCA"
Rscript "${R_SCRIPT_HVG_PCA}" \
        "${IN_OBJ}" \
        "${HPC_DIR}" \
        "${PARAMS_PCA}" \
        1> "${OUT_DIR}/hvg_pca.STDOUT" \
        2> "${OUT_DIR}/hvg_pca.STDERR"

echo "run clustering on automatically selected PCs"
Rscript "${R_SCRIPT_CLUSTERING}" \
        "${OUT_OBJ}" \
        "${HPC_DIR}" \
        "${PARAMS_CL}" \
        1> "${OUT_DIR}/clustering.STDOUT" \
        2> "${OUT_DIR}/clustering.STDERR"

echo "validate the clusters with the silhouette score"
Rscript "${R_SCRIPT_SILHOUETTE}" \
        "${OUT_OBJ}" \
        "${HPC_DIR}" \
        "${SILH_DIR}" \
        "${PARAMS_CL}" \
        1> "${OUT_DIR}/silhouette.STDOUT" \
        2> "${OUT_DIR}/silhouette.STDERR"

echo "compute the cluster composition"
Rscript "${R_SCRIPT_CLUSTER_COMPOSITION}" \
        "${OUT_OBJ}" \
        "${COMP_DIR}" \
        "${PARAMS_CL}" \
        1> "${OUT_DIR}/cluster_composition.STDOUT" \
        2> "${OUT_DIR}/cluster_composition.STDERR"

echo "generate PCA and TSNE plots"
Rscript "${R_SCRIPT_PLOT}" \
        "${OUT_OBJ}" \
        "${HPC_DIR}" \
        "${PLOT_DIR}" \
        "${PARAMS_CL}" \
        1> "${OUT_DIR}/plot.STDOUT" \
        2> "${OUT_DIR}/plot.STDERR"




exit
