
if [[ $# -lt "7" ]]
then
    echo "bash $0 <script_dir> <indir> <outdir> <exp> <sample> <mode> <k>"
    exit
fi

SCRIPT_DIR="$1"
INDIR="$2"
OUTDIR="$3"
EXP="$4"
SAMPLE="$5"
MODE="$6"
K="$7"

R_GRAPH="${SCRIPT_DIR}/1_compute_nn_graph.R"
R_PROP_MATRIX="${SCRIPT_DIR}/2_compute_propensity_matrix.R"
R_LINEAGE_ANCHOR="${SCRIPT_DIR}/3_detect_lineages.R"
R_LINEAGE_AUGMENT="${SCRIPT_DIR}/4_augment_lineages.R"

# input
CLONES="${OUTDIR}/${EXP}/top_clones/${SAMPLE}.txt"

# output
OUTDIR_G="${OUTDIR}/${EXP}/graphs/${SAMPLE}"
OUTDIR_P="${OUTDIR}/${EXP}/pair_propensity/${SAMPLE}"

[[ -d "${OUTDIR_G}" ]] || mkdir -p "${OUTDIR_G}"
[[ -d "${OUTDIR_P}" ]] || mkdir -p "${OUTDIR_P}"

Rscript "${R_SCRIPT_GRAPH}" "${INDIR}" "${OUTDIR_G}" ${SAMPLE} ${MODE} $K
Rscript "${R_PROP_MATRIX}" "${CLONES}" "${OUTDIR_G}" "${OUTDIR_P}" $K
Rscript "${R_LINEAGE_ANCHOR}" "${CLONES}" "${OUTDIR_P}"
Rscript "${R_LINEAGE_AUGMENT}" "${OUTDIR_G}" "${OUTDIR_P}" $K

exit

