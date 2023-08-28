

if [[ $# -lt "5" ]]
then
    echo "bash $0 <obj> <id> <hvg_dir> <dr> <outdir>"
    exit
fi

OBJ="$1"
ID=$2
HVG_DIR="$3"
HVG=$4
OUT="$5"

DIR=$( cd $(dirname $0) ; pwd ) 

SLOT="sample.name"
DR="pca_${HVG}"
NPC=$( cat "${HVG_DIR}/${HVG}/num_PCs.txt" )

Rscript "${DIR}/clone_distance_intra_sample.R" "${OBJ}" "${SLOT}" "${ID}" "${DR}" "${NPC}" "${OUT}"

exit
