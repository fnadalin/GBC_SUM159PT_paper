SCRIPT_DIR=$( cd "../../scripts/lineages" ; pwd )
DEA_SCRIPT_DIR=$( cd "../../scripts/seurat_analysis" ; pwd )
GIT_DIR=$( cd "../../../git/r_scripts" ; pwd )

CPT_INFO="cpt_late.tsv"
DTC_UNION="../clone_calling_sc/DRC/list_union.txt"

LIN_INFO="lineage_info.tsv"

EXP=( "1B_GEX" "1D_GEX" )
SAMPLES=( "T13,T17" "T13_2_Multiseq,T15_2_Multiseq" )

# MODES=( "mean.var.plot_disp1" "mean.var.plot_disp0.5" )
# K=( "40" "50" )
MODES=( "vst_top1000" "vst_top1000" )
K="40"

TOP_CLONES="20" # to anchor lineages

LOG2FC="0"
MAX_PVAL="0.05"

INDIR="gene_exp_correlation"

OUTDIR="lineages"
[[ -d ${OUTDIR} ]] || mkdir ${OUTDIR}

MEM="20000"

# compute the df with merged cpt counts and select the top 20 clones

ID="prepr_lineages"
COMMAND="${SCRIPT_DIR}/preprocessing.sh" 
ARGS=( "${CPT_INFO}" ${DTC_UNION} "${OUTDIR}" )
# bsub -J ${ID} -o ${ID}.STDOUT -e ${ID}.STDERR bash ${COMMAND} ${ARGS[@]}

# compute the nn graphs and the pair propensity matrices

depend=""
i=0
for exp in ${EXP[@]}
do
    IFS=","
    read -r -a samples <<< ${SAMPLES[$i]}
    IFS=$'\n' 
    
    mode=${MODES[$i]}

    for s in ${samples[@]}
    do  
        ID1="preprocessing_${exp}_$s"
        COMMAND="${SCRIPT_DIR}/preprocessing_v2.sh" 
        ARGS=( "${INDIR}" "${OUTDIR}" ${exp} $s ${TOP_CLONES} )
        bsub -J ${ID1} -w "done(${ID})" -o ${ID1}.STDOUT -e ${ID1}.STDERR bash ${COMMAND} ${ARGS[@]}

        OBJ_DIR="${exp}/treated_filt/hvg_pca_clust"
    
        ID2="lineages_${exp}_$s"
        COMMAND="${SCRIPT_DIR}/detect_lineages.sh"
        ARGS=( "${GIT_DIR}/sc_lineages" "${OBJ_DIR}" "${OUTDIR}" ${exp} $s ${mode} $K )
        bsub -J ${ID2} -w "done(${ID1})" -M ${MEM} -w "done(${ID1})" -o ${ID2}.STDOUT -e ${ID2}.STDERR bash ${COMMAND} ${ARGS[@]}
             
        depend="${depend}&&done(${ID2})"
    done
    
    i=$((i+1))
done
depend=${depend/&&/}

ID_POST="robust_lineages"
COMMAND="${SCRIPT_DIR}/postprocessing_v2.sh"
ARGS=( "${GIT_DIR}/gene_expression_analysis" "${GIT_DIR}/sc_lineages" "${LIN_INFO}" "${OUTDIR}" $K )
bsub -J ${ID_POST} -w "${depend}" -M ${MEM} -o ${ID_POST}.STDOUT -e ${ID_POST}.STDERR bash ${COMMAND} ${ARGS[@]}

ID_POST_P0="add_lineage_info_P0"
COMMAND="${SCRIPT_DIR}/add_lineage_info_P0.sh"
ARGS=( "${LIN_INFO}" "${OUTDIR}" )
bsub -J ${ID_POST_P0} -w "${depend}" -M ${MEM} -o ${ID_POST_P0}.STDOUT -e ${ID_POST_P0}.STDERR bash ${COMMAND} ${ARGS[@]}

ID_POST_P0="add_lineage_info_P0_ccRegress"
COMMAND="${SCRIPT_DIR}/add_lineage_info_P0_ccRegress.sh"
ARGS=( "${LIN_INFO}" "${OUTDIR}" )
bsub -J ${ID_POST_P0} -w "${depend}" -M ${MEM} -o ${ID_POST_P0}.STDOUT -e ${ID_POST_P0}.STDERR bash ${COMMAND} ${ARGS[@]}

i=0
depend=""
COMMAND="${DEA_SCRIPT_DIR}/4_DEA.sh" 
for exp in ${EXP[@]}
do
    IFS=","
    read -r -a samples <<< ${SAMPLES[$i]}
    IFS=$'\n' 
    
    OUT_OBJ="${exp}/treated_filt"

    ID="DEA_lineages_${exp}"
    ARGS=( "${GIT_DIR}/gene_expression_analysis" "${OUT_OBJ}" "lineage" )
    bsub -J ${ID} -w "done(${ID_POST})" -M ${MEM} -o ${ID}.STDOUT -e ${ID}.STDERR bash ${COMMAND} ${ARGS[@]}
    
    depend="${depend}&&done(${ID})"
    i=$((i+1))
done
depend=${depend/&&/}

OUTDIR="signatures/DEG_treated"
[ -d ${OUTDIR} ] || mkdir -p ${OUTDIR}

DEA1="${EXP[0]}/treated_filt/DEA/lineage/MAST/DEG_MAST_cl1-2.tsv"
DEA2="${EXP[1]}/treated_filt/DEA/lineage/MAST/DEG_MAST_cl1-2.tsv"

ID="signature_lineage1"
OUT="${OUTDIR}/${ID}.txt"

COMMAND="${DEA_SCRIPT_DIR}/deg_signature_up.sh"
ARGS=( ${LOG2FC} ${MAX_PVAL} "${OUT}" "${DEA1}" "${DEA2}" )
bsub -J ${ID} -w ${depend} -o ${ID}.STDOUT -e ${ID}.STDERR bash ${COMMAND} ${ARGS[@]}

ID="signature_lineage2"
OUT="${OUTDIR}/${ID}.txt"

COMMAND="${DEA_SCRIPT_DIR}/deg_signature_down.sh" 
ARGS=( -${LOG2FC} ${MAX_PVAL} "${OUT}" "${DEA1}" "${DEA2}" )
bsub -J ${ID} -w ${depend} -o ${ID}.STDOUT -e ${ID}.STDERR bash ${COMMAND} ${ARGS[@]}


exit

