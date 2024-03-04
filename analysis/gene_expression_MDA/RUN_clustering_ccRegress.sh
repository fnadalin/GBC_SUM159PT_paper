SCRIPT_DIR=$( cd "../../scripts/seurat_analysis" ; pwd )
GIT_DIR=$( cd "../../../git/r_scripts/gene_expression_analysis" ; pwd )
WDIR=$(pwd)

EXP=( "19924" )

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
    COMMAND="${SCRIPT_DIR}/2_load_matrix_ccRegress_noGBC.sh" 
    ARGS=( "${GIT_DIR}" "${FILE_LIST}" "${OUT_OBJ}" )
    bsub -J step2_${LABEL} -M ${MEM2} -o step2_${LABEL}.STDOUT -e step2_${LABEL}.STDERR bash ${COMMAND} ${ARGS[@]}

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
    # single sample
    batch_job ${exp} ${exp}
done


exit

