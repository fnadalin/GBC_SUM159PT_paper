
TOP_CLONES="20"

if [[ $# -lt "3" ]]
then
    echo "bash $0 <file_list> <dtc_union> <outdir>"
    exit
fi

CPT_INFO="$1"
DTC_UNION="$2"
OUTDIR="$3"

DIR=$( cd $(dirname $0) ; pwd ) 

# input
FILE_LIST=$( cut -f 2 "${CPT_INFO}" | tr '\n' ',' | sed "s/,$//g" )
SAMPLES=$( cut -f 3 "${CPT_INFO}" | tr '\n' ',' | sed "s/,$//g" )

# output
DF="${OUTDIR}/cpt.tsv"
CLONES="${OUTDIR}/top_clones.txt"

Rscript "${DIR}/merge_df.R" "${FILE_LIST}" "${DTC_UNION}" "${DF}" 
Rscript "${DIR}/top_clones.R" "${DF}" "${SAMPLES}" "${TOP_CLONES}" "${CLONES}"

exit

