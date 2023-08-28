

SCRIPT_DIR=$( cd $(dirname $0) ; pwd )

if [[ $# -lt "2" ]]
then
    echo ""
    echo "Usage: bash $0 <object.Rds> <out_dir>"
    echo ""
    exit
fi

OBJ="$1"
OUT_DIR="$2"

Rscript "${SCRIPT_DIR}/gene_clone_matrix.R" "${OBJ}" "${OUT_DIR}"

exit


