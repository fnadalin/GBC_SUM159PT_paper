
if [[ $# -lt "4" ]]
then
    echo "bash $0 <indir> <outdir> <meta> <params>"
    exit
fi

INDIR=$1
OUTDIR=$2
META=$3
PARAMS=$4

DIR=$( cd $(dirname $0) ; pwd ) 

Rscript "${DIR}/plot_umap_TIC_DRC.R" "${INDIR}" "${OUTDIR}" ${META} "${PARAMS}"

exit

