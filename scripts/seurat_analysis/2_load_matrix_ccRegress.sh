
module purge
module load r-4.1.0-gcc-9.3.0-wvnko7v gmp-6.1.2-gcc-9.3.0-hicntdj
R_LIBS_USER="/hps/software/users/marioni/francesca/R_libs"
export R_LIBS_USER

MIN_GENES="0"
MIN_CELLS="3"

if [[ $# -lt "4" ]]
then
    echo "bash $0 <script_dir> <matrix_list> <metadata> <in_dir> <out_dir>"
    exit
fi

SCRIPT_DIR="$1"
MATRIX_LIST="$2"
METADATA="$3"
IN_DIR="$4"
OUT_DIR="$5"

R_SCRIPT_LOAD_MATRIX="${SCRIPT_DIR}/Seurat_load_multi_matrix.R"
R_SCRIPT_ADD_META="${SCRIPT_DIR}/Seurat_add_meta_data.R"

# generate cell info
echo "generate cell info"
MATRIX_LIST_CELLS="${OUT_DIR}/matrix_dir_cells.tsv"
IFS=$'\n'
for line in $(cat ${MATRIX_LIST})
do
    S=$(echo $line | cut -f 1)
    echo "${line}	${IN_DIR}/cell_list_${S}.txt"
done > ${MATRIX_LIST_CELLS}

# create the object
echo "create the Seurat object"
Rscript "${R_SCRIPT_LOAD_MATRIX}" \
        --matrix_dir_list "${MATRIX_LIST_CELLS}" \
        --out_RDS "${OUT_DIR}/object.Rds" \
        --min_genes ${MIN_GENES} \
        --min_cells ${MIN_CELLS} \
        --vars_to_regress "S.Score,G2M.Score" \
        1> "${OUT_DIR}/Seurat_load_multi_matrix.STDOUT" \
        2> "${OUT_DIR}/Seurat_load_multi_matrix.STDERR"

# add metadata (clone info for each cell)
echo "add metadata (clone info for each cell)"
Rscript "${R_SCRIPT_ADD_META}" \
        --object "${OUT_DIR}/object.Rds" \
        --table "${METADATA}" \
        1> "${OUT_DIR}/Seurat_add_meta_data.STDOUT" \
        2> "${OUT_DIR}/Seurat_add_meta_data.STDERR"


exit


