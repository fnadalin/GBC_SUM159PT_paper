
if [[ $# -lt "5" ]]
then
    echo "bash $0 <cpm_bulk.tsv> <CB_metadata_sc.tsv> <sample_IDs_bulk> <sample_IDs_sc> <outdir>"
    exit
fi

CPM_BULK=$1
CB_META=$2
SAMPLES_BULK=$3
SAMPLES_SC=$4
OUTDIR=$5

DIR=$( cd $(dirname $0) ; pwd ) 

Rscript "${DIR}/compare_cpm_cellcount.R" "${CPM_BULK}" "${CB_META}" "${SAMPLES_BULK}" "${SAMPLES_SC}" "${OUTDIR}"

exit

