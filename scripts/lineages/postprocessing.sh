
if [[ $# -lt "2" ]]
then
    echo "bash $0 <file_list> <outdir>"
    exit
fi

LINEAGE_INFO="$1"
OUTDIR="$2"

DIR=$( cd $(dirname $0) ; pwd ) 

# input
FILE_LIST=$( cut -f 3 "${LINEAGE_INFO}" | tr '\n' ',' | sed "s/,$//g" )
EXP=( $( cut -f 1 "${LINEAGE_INFO}" | sort | uniq ) )

# output
LIN="${OUTDIR}/robust_lineages.tsv"

Rscript "${DIR}/robust_lineages.R" "${FILE_LIST}" "${LIN}"
for exp in ${EXP[@]}
do
    OBJECT="${exp}/treated_filt/hvg_pca_clust/object.Rds"
    Rscript "${DIR}/add_lineage_info.R" "${OBJECT}" "${LIN}"
done

exit

