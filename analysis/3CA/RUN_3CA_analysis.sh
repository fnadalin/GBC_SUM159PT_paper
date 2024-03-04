SCRIPT_DIR=$( cd "../../scripts/3CA" ; pwd )
DIR=$( cd "../../analysis/gene_expression" ; pwd )

TABLE_3CA="../../data/3CA/meta_programs_2023-07-13.tsv"

EXP=( "1B" "1D" )

OBJ=( "${DIR}/${EXP[0]}_GEX/P0_ccRegress_filt_2ndRound/hvg_pca_clust/object.Rds" \
      "${DIR}/${EXP[1]}_GEX/P0_ccRegress_filt/hvg_pca_clust/object.Rds" )

MODE=( "clusters_pca_vst_top1000_k40_res0.5" \
       "clusters_pca_vst_top1000_k40_res0.5" )

SIGN_DIR="${DIR}/signatures/DEG_P0_ccRegress"

SUBPOP_ID=( "6,2,5,0,1,3,4" "6,2,4,0,1,3,5" )
SUBPOP_NAME=( "1,2,3,4,5,6,7" "1,2,3,4,5,6,7" )

MEM="20000"

for (( i=0; i<${#EXP[@]}; i++ ))
do
    ID="MP_${EXP[$i]}"
    COMMAND="${SCRIPT_DIR}/meta_programs_vs_signatures.sh"
    ARGS=( ${TABLE_3CA} ${OBJ[$i]} ${MODE[$i]} ${SUBPOP_ID[$i]} ${SUBPOP_NAME[$i]} ${EXP[$i]} ${SIGN_DIR} )
    bsub -J ${ID} -M ${MEM} -o ${ID}.STDOUT -e ${ID}.STDERR bash ${COMMAND} ${ARGS[@]}
done

ID="MP_merged"
COMMAND="${SCRIPT_DIR}/heatmap_auc_merged.sh"
ARGS=( ${EXP[@]} "merged" )
bsub -J ${ID} -M ${MEM} -o ${ID}.STDOUT -e ${ID}.STDERR bash ${COMMAND} ${ARGS[@]}

exit

