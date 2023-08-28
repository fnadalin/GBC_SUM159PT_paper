
SCRIPT_DIR=$( cd "../../scripts/sister_cells" ; pwd )
GEX_DIR="../gene_expression"
MULTI_DIR="../multiome"

DIR=( "${GEX_DIR}" "${GEX_DIR}" )

EXP=( "1B_GEX" "1D_GEX" )
HVG_MODE=( "vst_top1000" "vst_top1000" )
CL_MODE=( "clusters_pca_vst_top1000_k40_res0.5" "clusters_pca_vst_top1000_k40_res0.5" )

MEM="20000"

for (( i=0; i<${#EXP[@]}-1; i++ ))
do
    OBJ1="${DIR[$i]}/${EXP[$i]}/P0_ccRegress_filt_2ndRound/hvg_pca_clust/object.Rds"
    HVG1="${DIR[$i]}/${EXP[$i]}/P0_ccRegress_filt_2ndRound/hvg_pca_clust/${HVG_MODE[$i]}/features.txt"

    for (( j=$((i+1)); j<${#EXP[@]}; j++ ))
    do
        OBJ2="${DIR[$j]}/${EXP[$j]}/P0_ccRegress_filt/hvg_pca_clust/object.Rds"
        HVG2="${DIR[$j]}/${EXP[$j]}/P0_ccRegress_filt/hvg_pca_clust/${HVG_MODE[$j]}/features.txt"
        
        # sister cells in cluster pairs
        ID="sisters_${EXP[$i]}_${EXP[$j]}_ccRegress"
        COMMAND="${SCRIPT_DIR}/cluster_clone_pairs.sh" 
        ARGS=( "${OBJ1}" "${OBJ2}" ${CL_MODE[$i]} ${CL_MODE[$j]} "${EXP[$i]}-${EXP[$j]}_ccRegress" )
        bsub -J ${ID} -M ${MEM} -o ${ID}.STDOUT -e ${ID}.STDERR bash ${COMMAND} ${ARGS[@]}
        
        # intra-clone distance on cca
        ID="sister_dist_cca_${EXP[$i]}_${EXP[$j]}_ccRegress"
        COMMAND="${SCRIPT_DIR}/clone_distance_withIntegration_cca.sh" 
        ARGS=( "${OBJ1}" "${OBJ2}" ${CL_MODE[$i]} ${CL_MODE[$j]} "${HVG1}" "${HVG2}" "${EXP[$i]}-${EXP[$j]}_ccRegress/dist_cca" )
        bsub -J ${ID} -M ${MEM} -o ${ID}.STDOUT -e ${ID}.STDERR bash ${COMMAND} ${ARGS[@]}
        
    done
    
done

exit
