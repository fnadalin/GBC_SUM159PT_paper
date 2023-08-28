

SCRIPT_DIR=$( cd $(dirname $0) ; pwd )

if [[ $# -lt "4" ]]
then
    echo ""
    echo "Usage: bash $0 <object.Rds> <out_dir> <seeds> <samples>"
    echo ""
    exit
fi

OBJ="$1"
OUT_DIR="$2"
SEEDS="$3"
SAMPLES="$4"

Rscript "${SCRIPT_DIR}/gene_clone_matrix_shuffled.R" "${OBJ}" "${OUT_DIR}" "${SEEDS}" "${SAMPLES}"

exit


