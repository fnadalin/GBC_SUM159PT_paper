

SCRIPT_DIR=$( cd $(dirname $0) ; pwd )

if [[ $# -lt "5" ]]
then
    echo ""
    echo "Usage: bash $0 <cl_dir> <mtx_dir> <corr_dir> <corr_perm_dir> <seeds>"
    echo ""
    exit
fi

OBJECT="$1"
SAMPLE="$2"
CL_DIR="$3"
CORR_DIR="$4"
CORR_PERM_DIR="$5"
SEEDS="$6"

Rscript "${SCRIPT_DIR}/gene_expression_correlation_shuffled.R" "${OBJECT}" "${SAMPLE}" "${CL_DIR}" "${CORR_PERM_DIR}" "${SEEDS}"
Rscript "${SCRIPT_DIR}/pvalue.R" "${CL_DIR}" "${CORR_DIR}" "${CORR_PERM_DIR}"

exit


