SCRIPT_DIR=$( cd "../../scripts/seurat_analysis" ; pwd )
SCRIPT_DIR=$( cd "../../scripts/MDA" ; pwd )
R_SCRIPT_DIR=$( cd "../../../git/r_scripts" ; pwd )
GIT_DIR=$( cd "${R_SCRIPT_DIR}/gene_expression_analysis" ; pwd )

EXP=( "19924" )

# select the clustering with the best silhouette score 
CL_MODE_P0=( "clusters_pca_vst_top1000_k40_res0.1" )
# ...and increase resolution 
CL_MODE_P0_HIGHRES=( "clusters_pca_vst_top1000_k40_res0.5" )

FILT=( "filt" )

MEM1="25000"
MEM2="15000"

DEA="4_DEA.sh"

batch_job () {

    SET=$1
    SLOT=$2
    SCRIPT=$3

    OUT_OBJ="${exp}/${SET}"
    
    # Run DEA 
    ID1="DEA_${exp}_${SET}_${SLOT}"
    COMMAND="${SCRIPT_DIR}/${SCRIPT}" 
    ARGS=( "${GIT_DIR}" "${OUT_OBJ}" ${SLOT} )
#    bsub -J ${ID1} -M ${MEM1} -o ${ID1}.STDOUT -e ${ID1}.STDERR bash ${COMMAND} ${ARGS[@]}
    
    # plot signatures
    ID2="sign_${exp}_${SET}_${SLOT}"
    COMMAND="${SCRIPT_DIR}/misc_stats_and_plots.sh" 
    ARGS=( "${OUT_OBJ}" ${SLOT} )
#    bsub -J ${ID2} -M ${MEM2} -w "done(${ID1})" -o ${ID2}.STDOUT -e ${ID2}.STDERR bash ${COMMAND} ${ARGS[@]}
    bsub -J ${ID2} -M ${MEM2} -o ${ID2}.STDOUT -e ${ID2}.STDERR bash ${COMMAND} ${ARGS[@]}  

}

# run differential expression analysis and plot signatures

i=0
for exp in ${EXP[@]}
do  
    # T0
    batch_job ${exp}_ccRegress_${FILT[$i]} ${CL_MODE_P0[$i]} ${DEA}
    batch_job ${exp}_ccRegress_${FILT[$i]} ${CL_MODE_P0_HIGHRES[$i]} ${DEA}

    i=$((i+1))
done




exit

