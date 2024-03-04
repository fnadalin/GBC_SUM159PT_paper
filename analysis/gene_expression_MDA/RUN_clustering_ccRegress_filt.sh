SCRIPT_DIR=$( cd "../../scripts/seurat_analysis" ; pwd )
MDA_SCRIPT_DIR=$( cd "../../scripts/MDA" ; pwd )
GIT_DIR=$( cd "../../../git/r_scripts/gene_expression_analysis" ; pwd )

WDIR=$(pwd)

EXP=( "19924" )

# extract this info from */STATS.tsv
# take the solution with top silhouette where the low UMI cluster is identified
MODE_P0="clusters_pca_vst_top1000_k30_res0.1"
CL_P0="3"

MEM="25000"

batch_job () {

    SET=$1
    LABEL=$2

    FILE_LIST="matrix_dir_${LABEL}.tsv"

    OUT_OBJ="${exp}/${SET}"
    [[ -d "${OUT_OBJ}" ]] || mkdir -p "${OUT_OBJ}"

    # Run the clustering
    COMMAND="${SCRIPT_DIR}/3_clustering.sh" 
    ARGS=( "${GIT_DIR}" "${OUT_OBJ}" "${WDIR}" )
#    bsub -J step3_${LABEL} -M ${MEM} -w "done(step1_${LABEL})" -o step3_${LABEL}.STDOUT -e step3_${LABEL}.STDERR bash ${COMMAND} ${ARGS[@]}
    bsub -J step3_${LABEL} -M ${MEM} -o step3_${LABEL}.STDOUT -e step3_${LABEL}.STDERR bash ${COMMAND} ${ARGS[@]}
    
    # collect statistics on generated clustering solutions
    COMMAND="${SCRIPT_DIR}/collect_cl_stat.sh" 
    ARGS=( "${WDIR}" "${OUT_OBJ}" )
    bsub -J step4_${LABEL} -M ${MEM} -w "done(step3_${LABEL})" -o step4_${LABEL}.STDOUT -e step4_${LABEL}.STDERR bash ${COMMAND} ${ARGS[@]}
}

for exp in ${EXP[@]}
do
    
    # subset the cells (exclude low UMI cluster)
    [[ -d "${exp}" ]] || mkdir "${exp}"
    COMMAND="${MDA_SCRIPT_DIR}/cell_filtering_ccRegress.sh" 
    ARGS=( "${SCRIPT_DIR}" "${exp}" ${MODE_P0}" ${CL_P0}" )
#    bsub -J step1_${exp} -M ${MEM} -o step1_${exp}.STDOUT -e step1_${exp}.STDERR bash ${COMMAND} ${ARGS[@]}

    batch_job ${exp}_ccRegress_filt ${exp}
done


exit

