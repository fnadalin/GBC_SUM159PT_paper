
DIR=$( cd $(dirname $0) ; pwd )

if [[ $# -lt "4" ]]
then
    echo "bash $0 <wes> <sc_name> <sc_dir> <cl_info>"
    exit
fi

WES_EXP_FILE=$1
SC_EXP_NAME=$2
SC_EXP_DIR=$3
CL_INFO=$4

Rscript ${DIR}/granges_comparison_scATAC_double.R "${WES_EXP_FILE}" "${SC_EXP_NAME}" "${SC_EXP_DIR}" "${CL_INFO}"
Rscript ${DIR}/heatmap_granges_comparison_scATAC_double.R "${SC_EXP_NAME}"

exit

