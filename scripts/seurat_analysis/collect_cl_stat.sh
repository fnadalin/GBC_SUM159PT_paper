

SCRIPT_DIR=$( cd $(dirname $0) ; pwd )

if [[ $# -lt "2" ]]
then
    echo ""
    echo "Usage: bash $0 <param_dir> <out_dir>"
    echo ""
    exit
fi

PARAM_DIR="$1"
OUT_DIR="$2"

PARAM="${PARAM_DIR}/param_clust.txt"
PCA_DIR="${OUT_DIR}/hvg_pca_clust"
SILH_DIR="${OUT_DIR}/silhouette"
OBJ="${PCA_DIR}/object.Rds"
OUT="${OUT_DIR}/STATS.tsv"

Rscript "${SCRIPT_DIR}/collect_cl_stat.R" ${PARAM} "${PCA_DIR}" "${SILH_DIR}" "${OBJ}" "${OUT}"

exit


