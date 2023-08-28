

if [[ $# -lt "4" ]]
then
    echo ""
    echo "Usage: bash $0 <metadata> <outdir> <ctrl> <case>"
    echo ""
    exit
fi

METADATA="$1"
OUT_DIR="$2"
CTRL="$3"
CASE="$4"

SCRIPT_DIR=$( cd $(dirname $0) ; pwd )

R_SCRIPT_COUNT="${SCRIPT_DIR}/clone_count.R"
R_SCRIPT_NORM="${SCRIPT_DIR}/normalize.R"
R_SCRIPT_EXP="${SCRIPT_DIR}/expansion_rate.R"

Rscript "${R_SCRIPT_COUNT}" "${METADATA}" "${OUT_DIR}/counts.tsv"
Rscript "${R_SCRIPT_NORM}" "${OUT_DIR}/counts.tsv" "${OUT_DIR}/cpt.tsv"
Rscript "${R_SCRIPT_EXP}" "${OUT_DIR}/cpt.tsv" "${OUT_DIR}/log2fc.tsv" ${CTRL} ${CASE}

exit


