
if [[ $# -lt "4" ]]
then
    echo "bash $0 <script_dir> <matrix_list> <out_dir> <out_dir_sc>"
    exit
fi

SCRIPT_DIR="$1"
MATRIX_DIR_LIST="$2"
OUT_DIR="$3"
OUT_DIR_SC="$4"

LOCAL_SCRIPT_DIR=$( cd $(dirname $0) ; pwd ) 

R_SCRIPT_LOAD_OBJECT="${SCRIPT_DIR}/Seurat_load_multi_matrix.R"
R_SCRIPT_FILTER_CELLS="${LOCAL_SCRIPT_DIR}/select_epcam_and_filter.R"
R_SCRIPT_PREPR="${LOCAL_SCRIPT_DIR}/scissor_preprocessing.R"

[[ -d "${OUT_DIR}" ]] || mkdir -p "${OUT_DIR}"
[[ -d "${OUT_DIR_SC}" ]] || mkdir -p "${OUT_DIR_SC}"

Rscript ${R_SCRIPT_LOAD_OBJECT} --matrix_dir_list ${MATRIX_DIR_LIST} --out_RDS ${OUT_DIR}/object_all.Rds
Rscript ${R_SCRIPT_FILTER_CELLS} ${OUT_DIR}/object_all.Rds ${OUT_DIR}
Rscript ${R_SCRIPT_PREPR} ${OUT_DIR}/object.Rds ${OUT_DIR_SC}/sc_dataset.Rds

exit

