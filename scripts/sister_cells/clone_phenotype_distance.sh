

if [[ $# -lt "5" ]]
then
    echo "bash $0 <obj> <slot> <hvg_dir> <dr> <outdir>"
    exit
fi

OBJ="$1"
SLOT=$2
HVG_DIR="$3"
HVG=$4
OUT="$5"

DIR=$( cd $(dirname $0) ; pwd ) 

DR="pca_${HVG}"
NPC=$( cat "${HVG_DIR}/${HVG}/num_PCs.txt" )

Rscript "${DIR}/clone_phenotype_distance.R" "${OBJ}" "${SLOT}" "${DR}" "${NPC}" "${OUT}"

exit
