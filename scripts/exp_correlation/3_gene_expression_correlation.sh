

SCRIPT_DIR=$( cd $(dirname $0) ; pwd )

if [[ $# -lt "3" ]]
then
    echo ""
    echo "Usage: bash $0 <log2fc.tsv> <mtx_dir> <out_dir>"
    echo ""
    exit
fi

CL_DIR="$1"
MTX_DIR="$2"
OUT_DIR="$3"

Rscript "${SCRIPT_DIR}/gene_expression_correlation.R" "${CL_DIR}" "${MTX_DIR}" "${OUT_DIR}"

exit


