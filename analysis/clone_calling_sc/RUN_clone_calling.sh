SCRIPT_DIR=$( cd "../../scripts/clone_calling" ; pwd )
SCRIPT_DIR_BULK=$( cd "../../scripts/clone_estimate" ; pwd )
GIT_DIR=$( cd "../../../git/r_scripts/clones_analysis" ; pwd )

EXP=( "1B_GEX" "1C_GEX" "1D_GEX" "1E_GEX" )
MEM1="10000"
MEM2="2000"

depend=""
for exp in ${EXP[@]}
do
    FILE_LIST="matrix_dir_${exp}.tsv"
    OUT="${exp}"
    PARAMS="params_${exp}"

    # Compute statistics on GBC read counts and UMI counts
    [[ -d "${OUT}" ]] || mkdir "${OUT}"
    COMMAND="${SCRIPT_DIR}/clone_analysis_step1.sh" 
    ARGS=( "${GIT_DIR}" "${FILE_LIST}" "${OUT}" )
    bsub -J step1_${exp} -M ${MEM1} -o step1_${exp}.STDOUT -e step1_${exp}.STDERR bash ${COMMAND} ${ARGS[@]}

    # Based on the GBC read count, determine which GBC are expressed in each cell
    # Then, classify the cells as single-infection, co-infection, or no-infection events, and identify the doublets
    # Finally, determine the clones as GBC sets and possibly assign each cell to a clone. 
    COMMAND="${SCRIPT_DIR}/clone_analysis_step2.sh"
    ARGS=( "${GIT_DIR}" "${PARAMS}" "${OUT}" )
    bsub -J step2_${exp} -M ${MEM2} -w "done(step1_${exp})" -o step2_${exp}.STDOUT -e step2_${exp}.STDERR bash ${COMMAND} ${ARGS[@]}
    depend="${depend}&&done(step2_${exp})"

done
depend=${depend/&&/}

[[ -d "DRC" ]] || mkdir "DRC"

# call DRCs
SUFFIX="merged_samples/surv_clones_union.tsv"
COMMAND="${SCRIPT_DIR}/call_DRC.sh"
ARGS=( "1B_GEX/${SUFFIX},1D_GEX/${SUFFIX}" "DRC/list.txt" )
bsub -J call_DRC -M ${MEM2} -w "${depend}" -o call_DRC.STDOUT -e call_DRC.STDERR bash ${COMMAND} ${ARGS[@]}

COMMAND="${SCRIPT_DIR}/call_DRC_union.sh"
ARGS=( "1B_GEX/${SUFFIX},1D_GEX/${SUFFIX}" "DRC/list_union.txt" )
bsub -J call_DRC_union -M ${MEM2} -w "${depend}" -o call_DRC_union.STDOUT -e call_DRC_union.STDERR bash ${COMMAND} ${ARGS[@]}

# compute fisher test
SUFFIX="merged_samples/clone_sample_CB_count.tsv"
COMMAND="${SCRIPT_DIR}/compute_fisher_wilcox.sh"
ARGS=( "1B_GEX/${SUFFIX},1D_GEX/${SUFFIX},1C_GEX/${SUFFIX},1E_GEX/${SUFFIX}" "all_clones_exp1" "files_exp_1_GBC_sc.tsv" )
bsub -J compute_fisher -M ${MEM2} -w "${depend}" -o compute_fisher.STDOUT -e compute_fisher.STDERR bash ${COMMAND} ${ARGS[@]}

exit

