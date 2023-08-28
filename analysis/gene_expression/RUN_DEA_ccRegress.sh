SCRIPT_DIR=$( cd "../../scripts/seurat_analysis" ; pwd )
R_SCRIPT_DIR=$( cd "../../../git/r_scripts" ; pwd )
GIT_DIR=$( cd "${R_SCRIPT_DIR}/gene_expression_analysis" ; pwd )

EXP=( "1B_GEX" "1D_GEX" )

# select the clustering with the best silhouette score & the S100A4 cluster detected
CL_MODE_P0=( "clusters_pca_vst_top1000_k40_res0.1" "clusters_pca_vst_top1000_k40_res0.1" )
# ...and increase resolution (for use in sister cell testing)
CL_MODE_P0_HIGHRES=( "clusters_pca_vst_top1000_k40_res0.5" "clusters_pca_vst_top1000_k40_res0.5" )

FILT=( "filt_2ndRound" "filt" )

CL_ID_EXP1B=( "6" "2" "5" )
CL_ID_EXP1D=( "6" "2" "4" )

CLONE_CLASS=( "DRC_bulk" "DRC_bulk_wilcox" "DRC_bulk_fisher" "DRC_bulk_gam" \
              "TIC" "TIC_wilcox" "TIC_fisher" "TIC_gam" "DRC_sc" )

SAMPLES_MIDDLE=( "T5_2_Multiseq" "T7_2_Multiseq" "T9_2_Multiseq" )

MEM="25000"
MEM_FA="10000"

MIN_LOG2FC="0"
MAX_PVAL="0.05"
TOP_GENES="100"


DEA="4_DEA.sh"
DEA_COND="4_DEA_cond.sh"

batch_job () {

    SET=$1
    SLOT=$2
    SCRIPT=$3

    OUT_OBJ="${exp}/${SET}"
    
    # Run DEA 
    ID="DEA_${exp}_${SET}_${SLOT}"
    COMMAND="${SCRIPT_DIR}/${SCRIPT}" 
    ARGS=( "${GIT_DIR}" "${OUT_OBJ}" ${SLOT} )
    bsub -J ${ID} -M ${MEM} -o ${ID}.STDOUT -e ${ID}.STDERR bash ${COMMAND} ${ARGS[@]}

}

# run differential expression analysis

i=0
for exp in ${EXP[@]}
do  
    # T0
    batch_job P0_ccRegress_${FILT[$i]} ${CL_MODE_P0[$i]} ${DEA}
    batch_job P0_ccRegress_${FILT[$i]} ${CL_MODE_P0_HIGHRES[$i]} ${DEA}
    for c in ${CLONE_CLASS[@]}
    do
        batch_job P0_ccRegress_${FILT[$i]} $c ${DEA}
    done
    
    # treated (T > 0)
    for c in ${CLONE_CLASS[@]}
    do
        batch_job treated_ccRegress_filt $c ${DEA_COND}
    done
    
    i=$((i+1))

done

# extract signatures for paired clusters

ID1="DEA_${EXP[0]}_P0_filt_${CL_MODE_P0_HIGHRES[0]}"
ID2="DEA_${EXP[1]}_P0_filt_${CL_MODE_P0_HIGHRES[1]}"

for (( i=0; i<${#CL_ID_EXP1B[@]}; i++ ))
do

    OUTDIR="signatures/DEG_P0_ccRegress"
    [ -d ${OUTDIR} ] || mkdir -p ${OUTDIR}

    DEA1="${EXP[0]}/P0_ccRegress_${FILT[0]}/DEA/${CL_MODE_P0_HIGHRES[0]}/MAST/DEG_MAST_cl${CL_ID_EXP1B[i]}-all.tsv"
    DEA2="${EXP[1]}/P0_ccRegress_${FILT[1]}/DEA/${CL_MODE_P0_HIGHRES[1]}/MAST/DEG_MAST_cl${CL_ID_EXP1D[i]}-all.tsv"
    
    ID="signature$((i+1))"
    OUT="${OUTDIR}/${ID}.txt"

    COMMAND="${SCRIPT_DIR}/deg_signature_up.sh" 
    ARGS=( ${MIN_LOG2FC} ${MAX_PVAL} "${OUT}" "${DEA1}" "${DEA2}" )
    bsub -J ${ID} -w "done(${ID1}) && done(${ID2})" -o ${ID}.STDOUT -e ${ID}.STDERR bash ${COMMAND} ${ARGS[@]}

    OUT_FA="${OUTDIR}/GO/${ID}"
    [ -d "${OUTDIR}/GO" ] || mkdir "${OUTDIR}/GO"

    IDF="${ID}_fa"
    COMMAND="${SCRIPT_DIR}/functional_annotation.sh"
    ARGS=( ${R_SCRIPT_DIR} ${OUT} ${OUT_FA} )
    bsub -J ${IDF} -w "done(${ID})" -o ${IDF}.STDOUT -e ${IDF}.STDERR bash ${COMMAND} ${ARGS[@]}
    
done

# extract signatures for clone classes (P0)

for (( i=0; i<${#CLONE_CLASS[@]}; i++ ))
do
    ID1="DEA_${EXP[0]}_P0_ccRegress_${FILT[0]}_${CLONE_CLASS[$i]}"
    ID2="DEA_${EXP[1]}_P0_ccRegress_${FILT[1]}_${CLONE_CLASS[$i]}"

    OUTDIR="signatures/DEG_P0_ccRegress"
    [ -d ${OUTDIR} ] || mkdir -p ${OUTDIR}

    DEG_SUFFIX="DEA/${CLONE_CLASS[$i]}/MAST/DEG_MAST_cl1-all.tsv"
    DEA_LIST=( "${EXP[0]}/P0_ccRegress_${FILT[0]}/${DEG_SUFFIX}" "${EXP[1]}/P0_ccRegress_${FILT[1]}/${DEG_SUFFIX}" )
    
    ID="signature${CLONE_CLASS[$i]}"
    OUT="${OUTDIR}/${ID}.txt"

    COMMAND="${SCRIPT_DIR}/deg_signature_up.sh" 
    ARGS=( ${MIN_LOG2FC} ${MAX_PVAL} "${OUT}" ${DEA_LIST[@]} )
    bsub -J ${ID} -w "done(${ID1}) && done(${ID2})" -o ${ID}.STDOUT -e ${ID}.STDERR bash ${COMMAND} ${ARGS[@]}

done

# extract signatures for clone classes (treated)

for (( i=0; i<${#CLONE_CLASS[@]}; i++ ))
do
    ID1="DEA_${EXP[1]}_treated_ccRegress_filt_${CLONE_CLASS[$i]}"

    OUTDIR="signatures/DEG_treated_ccRegress"
    [ -d ${OUTDIR} ] || mkdir -p ${OUTDIR}

    DEG_DIR="${EXP[1]}/treated_ccRegress_filt/DEA/${CLONE_CLASS[$i]}"
    DEA_LIST=()
    for s in ${SAMPLES_MIDDLE[@]}
    do
        DEA_LIST=( ${DEA_LIST[@]} "${DEG_DIR}/$s/MAST/DEG_MAST_cl1-all.tsv" )
    done
    
    ID="signature${CLONE_CLASS[$i]}"
    OUT="${OUTDIR}/${ID}_d5-9_up.txt"

    COMMAND="${SCRIPT_DIR}/deg_signature_up.sh" 
    ARGS=( ${MIN_LOG2FC} ${MAX_PVAL} "${OUT}" ${DEA_LIST[@]} )
    bsub -J ${ID} -w "done(${ID1})" -o ${ID}.STDOUT -e ${ID}.STDERR bash ${COMMAND} ${ARGS[@]}

done

exit

