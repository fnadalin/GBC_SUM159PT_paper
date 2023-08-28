
SCRIPT_DIR=$( cd "../../scripts/sister_cells" ; pwd )
DIR="../gene_expression"

EXP=( "1B_GEX" "1D_GEX" )
HVG_MODE_P0=( "vst_top1000" "vst_top1000" )
HVG_MODE_TREATED=( "vst_top1000" "vst_top5000" )

SLOT="DRC_bulk_fisher"

MEM="16000"

DIR_P0="P0_ccRegress_filt_2ndRound"
for (( i=0; i<${#EXP[@]}; i++ ))
do

    OBJ_DIR="${DIR}/${EXP[$i]}/${DIR_P0}/hvg_pca_clust"
    OBJ="${OBJ_DIR}/object.Rds"
    HVG_DIR="${OBJ_DIR}/${HVG_MODE[$i]}"

    # intra-sample distance 
    ID="DTC_dist_treated_${EXP[$i]}"
    COMMAND="${SCRIPT_DIR}/clone_phenotype_distance.sh" 
    ARGS=( "${OBJ}" "${SLOT}" "${OBJ_DIR}" "${HVG_MODE_P0[$i]}" "${EXP[$i]}/dist_${SLOT}_P0" )
    bsub -J ${ID} -M ${MEM} -o ${ID}.STDOUT -e ${ID}.STDERR bash ${COMMAND} ${ARGS[@]}

    OBJ_DIR="${DIR}/${EXP[$i]}/treated_filt/hvg_pca_clust"
    OBJ="${OBJ_DIR}/object.Rds"
    HVG_DIR="${OBJ_DIR}/${HVG_MODE[$i]}"

    # intra-sample distance 
    ID="DTC_dist_treated_${EXP[$i]}"
    COMMAND="${SCRIPT_DIR}/clone_phenotype_distance.sh" 
    ARGS=( "${OBJ}" "${SLOT}" "${OBJ_DIR}" "${HVG_MODE_TREATED[$i]}" "${EXP[$i]}/dist_${SLOT}_treated" )
    bsub -J ${ID} -M ${MEM} -o ${ID}.STDOUT -e ${ID}.STDERR bash ${COMMAND} ${ARGS[@]}
    
    DIR_P0="P0_ccRegress_filt"
    
done

exit
