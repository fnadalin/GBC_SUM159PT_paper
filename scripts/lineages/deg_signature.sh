

SCRIPT_DIR=$( cd $(dirname $0) ; pwd )

LOG2FC="0"

if [[ $# -lt "3" ]]
then
    echo ""
    echo "Usage: bash $0 <dea1> <dea2> <out>"
    echo ""
    exit
fi

DEA1="$1"
DEA2="$2"
OUT="$3"
[[ $# -le "3" ]] || LOG2FC="$4"

Rscript "${SCRIPT_DIR}/deg_signature_up.R" "${DEA1}" "${DEA2}" "${OUT}" "${LOG2FC}"
Rscript "${SCRIPT_DIR}/deg_signature_down.R" "${DEA1}" "${DEA2}" "${OUT}" "-${LOG2FC}"

exit


