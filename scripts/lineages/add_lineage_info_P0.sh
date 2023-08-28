
if [[ $# -lt "2" ]]
then
    echo "bash $0 <file_list> <outdir>"
    exit
fi

LINEAGE_INFO="$1"
OUT_DIR="$2"

DIR=$( cd $(dirname $0) ; pwd ) 

# input
EXP=( $( cut -f 1 "${LINEAGE_INFO}" | sort | uniq ) )
LIN="${OUT_DIR}/lineages_merge.tsv"

for exp in ${EXP[@]}
do
    IN_DIR="${exp}/P0_filt"
    HPC_DIR="${IN_DIR}/hvg_pca_clust"
    IN_OBJ="${HPC_DIR}/object.Rds"
    
    Rscript "${DIR}/add_lineage_info.R" "${IN_OBJ}" "${LIN}"

done


exit

