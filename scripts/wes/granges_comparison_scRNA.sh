
DIR=$( cd $(dirname $0) ; pwd )

if [[ $# -lt "7" ]]
then
    echo "bash $0 <wes> <sc_name> <sc_obj> <sc_samples> <features> <lin1_marks> <lin2_marks>"
    exit
fi

WES_EXP_FILE=$1
SC_EXP_NAME=$2
SC_EXP_OBJ=$3
SC_SAMPLES=$4
FEATURES=$5
LIN1_MARKS=$6
LIN2_MARKS=$7

Rscript ${DIR}/granges_comparison_scRNA_intersect.R "${WES_EXP_FILE}" "${SC_EXP_NAME}" "${SC_EXP_OBJ}" "${SC_SAMPLES}" \
        "${FEATURES}" "${LIN1_MARKS}" "${LIN2_MARKS}"
Rscript ${DIR}/heatmap_granges_comparison_scRNA_intersect.R "${SC_EXP_NAME}" "${SC_SAMPLES}"

exit

