SCRIPT_DIR=$( cd "../../scripts/Visvader" ; pwd )
GIT_DIR=$( cd "../../../git/r_scripts/gene_expression_analysis" ; pwd )
DIR=$( cd "../../analysis/gene_expression" ; pwd )

EXP=( "1B" "1D" )

OBJ=( "${DIR}/${EXP[0]}_GEX/P0_ccRegress_filt_2ndRound/hvg_pca_clust" \
      "${DIR}/${EXP[1]}_GEX/P0_ccRegress_filt/hvg_pca_clust" )

MODE=( "clusters_pca_vst_top1000_k40_res0.5" \
       "clusters_pca_vst_top1000_k40_res0.5" )
       
CL1=( "6" "6" )
CL2=( "2" "2" )
CL3=( "5" "4" )

MEM="20000"

# prepare Seurat object and scissor object
ID="1_prepare_object"
COMMAND="${SCRIPT_DIR}/1_prepare_object.sh"
ARGS=( "${GIT_DIR}" "matrix_dir_list.tsv" "input" "scissor/in" )
bsub -J ${ID} -M ${MEM} -o ${ID}.STDOUT -e ${ID}.STDERR bash ${COMMAND} ${ARGS[@]}

# calculate pseudobulk profiles and run Scissor
COMMAND="${SCRIPT_DIR}/2_run_Scissor.sh"
ID1="2_run_Scissor_S1"
ARGS=( "${OBJ[0]}/object.Rds" "${OBJ[1]}/object.Rds" ${MODE[0]} ${MODE[1]} ${CL1[0]} ${CL1[1]} ${EXP[0]} ${EXP[1]} "scissor" "S1" )
bsub -J ${ID1} -M ${MEM} -w "done(${ID})" -o ${ID1}.STDOUT -e ${ID1}.STDERR bash ${COMMAND} ${ARGS[@]}
ID2="2_run_Scissor_S2"
ARGS=( "${OBJ[0]}/object.Rds" "${OBJ[1]}/object.Rds" ${MODE[0]} ${MODE[1]} ${CL2[0]} ${CL2[1]} ${EXP[0]} ${EXP[1]} "scissor" "S2" )
bsub -J ${ID2} -M ${MEM} -w "done(${ID})" -o ${ID2}.STDOUT -e ${ID2}.STDERR bash ${COMMAND} ${ARGS[@]}
ID3="2_run_Scissor_S3"
ARGS=( "${OBJ[0]}/object.Rds" "${OBJ[1]}/object.Rds" ${MODE[0]} ${MODE[1]} ${CL3[0]} ${CL3[1]} ${EXP[0]} ${EXP[1]} "scissor" "S3" )
bsub -J ${ID3} -M ${MEM} -w "done(${ID})" -o ${ID3}.STDOUT -e ${ID3}.STDERR bash ${COMMAND} ${ARGS[@]}

depend="done(${ID1})&&done(${ID2})&&done(${ID3})"

# run Scissor postprocessing
COMMAND="${SCRIPT_DIR}/3_run_Scissor_postprocessing.sh"
ID="3_scissor_postprocessing"
ARGS=( "input/object.Rds" "scissor" "S1,S2,S3" )
bsub -J ${ID} -M ${MEM} -w "${depend}" -o ${ID}.STDOUT -e ${ID}.STDERR bash ${COMMAND} ${ARGS[@]}

exit

