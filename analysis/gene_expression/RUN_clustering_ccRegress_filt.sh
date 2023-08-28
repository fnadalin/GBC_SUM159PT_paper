SCRIPT_DIR=$( cd "../../scripts/seurat_analysis" ; pwd )
GIT_DIR=$( cd "../../../git/r_scripts/gene_expression_analysis" ; pwd )

WDIR=$(pwd)

EXP=( "1B_GEX" "1D_GEX" )

# extract this info from */STATS.tsv
# take the solution with top silhouette where the low UMI cluster is identified
MODE_P0=( "clusters_pca_mean.var.plot_disp1.5_k30_res0.1" "clusters_pca_mean.var.plot_disp1.5_k50_res0.1" )
CL_P0=( "2" "2" )
MODE_T=( "clusters_pca_mean.var.plot_disp1.5_k30_res0.3" "clusters_pca_mean.var.plot_disp0.5_k40_res0.5" )
CL_T=( "6" "5" )

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
    bsub -J step3_${LABEL} -M ${MEM} -w "done(step1_${exp})" -o step3_${LABEL}.STDOUT -e step3_${LABEL}.STDERR bash ${COMMAND} ${ARGS[@]}
    
    # collect statistics on generated clustering solutions
    COMMAND="${SCRIPT_DIR}/collect_cl_stat.sh" 
    ARGS=( "${WDIR}" "${OUT_OBJ}" )
    bsub -J step4_${LABEL} -M ${MEM} -w "done(step3_${LABEL})" -o step4_${LABEL}.STDOUT -e step4_${LABEL}.STDERR bash ${COMMAND} ${ARGS[@]}
}

i=0
for exp in ${EXP[@]}
do
    
    # subset the cells (exclude low UMI cluster)
    [[ -d "${exp}" ]] || mkdir "${exp}"
    COMMAND="${SCRIPT_DIR}/cell_filtering_ccRegress.sh" 
    ARGS=( "${exp}" ${MODE_P0[$i]}" "${MODE_T[$i]}" ${CL_P0[$i]}" ${CL_T[$i]} )
    bsub -J step1_${exp} -M ${MEM} -o step1_${exp}.STDOUT -e step1_${exp}.STDERR bash ${COMMAND} ${ARGS[@]}

    # all samples
    batch_job all_ccRegress_filt ${exp}
    
    # T0
    batch_job P0_ccRegress_filt ${exp}_P0
    
    # treated (T > 0)
    batch_job treated_ccRegress_filt ${exp}_treated

    i=$((i+1))
done


exit

