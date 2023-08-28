
SCRIPT_DIR=$( cd "../../scripts/sister_cells" ; pwd )
DIR="../gene_expression"

EXP=( "1B_GEX" "1D_GEX" )
SAMPLES=( "T13,T17" "T11_2_Multiseq,T13_2_Multiseq,T15_2_Multiseq" )
HVG_MODE=( "vst_top1000" "vst_top5000" )

MEM="16000"

for (( i=0; i<${#EXP[@]}; i++ ))
do
    OBJ_DIR="${DIR}/${EXP[$i]}/treated_filt/hvg_pca_clust"
    OBJ="${OBJ_DIR}/object.Rds"
    HVG_DIR="${OBJ_DIR}/${HVG_MODE[$i]}"

    samples=( ${SAMPLES[$i]//,/ } )

    for s in ${samples[@]}
    do
        # intra-sample distance 
        ID="sister_dist_treated_${EXP[$i]}"
        COMMAND="${SCRIPT_DIR}/clone_distance_intra_sample.sh" 
        ARGS=( "${OBJ}" $s "${OBJ_DIR}" "${HVG_MODE[$i]}" "${EXP[$i]}/dist_$s" )
        bsub -J ${ID} -M ${MEM} -o ${ID}.STDOUT -e ${ID}.STDERR bash ${COMMAND} ${ARGS[@]}
    done
    
done

exit
