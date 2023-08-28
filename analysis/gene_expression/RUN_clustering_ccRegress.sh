SCRIPT_DIR=$( cd "../../scripts/seurat_analysis" ; pwd )
GIT_DIR=$( cd "../../../git/r_scripts/gene_expression_analysis" ; pwd )
IN_DIR=$( cd "../clone_calling_sc" ; pwd )
WDIR=$(pwd)

EXP=( "1B_GEX" "1D_GEX" )

CLONE_INFO="clone_info.tsv"

MEM1="2000"
MEM2="25000"
MEM3="20000"

batch_job () {

    SET=$1
    LABEL=$2

    FILE_LIST="matrix_dir_${LABEL}.tsv"

    OUT_OBJ="${exp}/${SET}_ccRegress"
    [[ -d "${OUT_OBJ}" ]] || mkdir -p "${OUT_OBJ}"

    # Generate the seurat object 
    # Use only the cells in the metadata file (i.e. labelled with GBCs and quality pass)
    COMMAND="${SCRIPT_DIR}/2_load_matrix_ccRegress.sh" 
    ARGS=( "${GIT_DIR}" "${FILE_LIST}" "${exp}/CB_metadata.tsv" "${exp}" "${OUT_OBJ}" )
    bsub -J step2_${LABEL} -M ${MEM2} -w "done(step1_${exp})" -o step2_${LABEL}.STDOUT -e step2_${LABEL}.STDERR bash ${COMMAND} ${ARGS[@]}

    # Run the clustering
    COMMAND="${SCRIPT_DIR}/3_clustering.sh" 
    ARGS=( "${GIT_DIR}" "${OUT_OBJ}" "${WDIR}" )
    bsub -J step3_${LABEL} -M ${MEM2} -w "done(step2_${LABEL})" -o step3_${LABEL}.STDOUT -e step3_${LABEL}.STDERR bash ${COMMAND} ${ARGS[@]}
    
    # collect statistics on generated clustering solutions
    COMMAND="${SCRIPT_DIR}/collect_cl_stat.sh" 
    ARGS=( "${WDIR}" "${OUT_OBJ}" )
    bsub -J step4_${LABEL} -M ${MEM3} -w "done(step3_${LABEL})" -o step4_${LABEL}.STDOUT -e step4_${LABEL}.STDERR bash ${COMMAND} ${ARGS[@]}
}

for exp in ${EXP[@]}
do
    FILE_LIST="matrix_dir_${exp}.tsv"
    IN="${IN_DIR}/${exp}"
    
    # Generate metadata for clone-labelled cells
    [[ -d "${exp}" ]] || mkdir "${exp}"
    COMMAND="${SCRIPT_DIR}/1_create_metadata.sh" 
    ARGS=( "${GIT_DIR}" "${FILE_LIST}" "${IN}" "${exp}" "${CLONE_INFO}" )
    bsub -J step1_${exp} -M ${MEM1} -o step1_${exp}.STDOUT -e step1_${exp}.STDERR bash ${COMMAND} ${ARGS[@]}

    # all samples
    batch_job all ${exp}
    
    # T0
    batch_job P0 ${exp}_P0
    
    # treated (T > 0)
    batch_job treated ${exp}_treated

done


exit

