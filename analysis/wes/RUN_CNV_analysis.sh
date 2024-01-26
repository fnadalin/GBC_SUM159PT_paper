
SCRIPT_DIR="../../scripts/wes"
CELL_CALLING_DIR="../cell_calling"
GENE_EXPR_DIR="../gene_expression"
MULTIOME_DATA_DIR="../../data/Multiome"

WES_EXP_FILE="in/CNV_intersection_80perc.bed"
SC_EXP_NAME_ARC=( "1C_ARC" "1E_ARC" )
SC_EXP_NAME_GEX=( "1B_GEX" "1D_GEX" )
SC_EXP_DIR=( "${CELL_CALLING_DIR}/1C_ARC/P0_Multiome/outs/filtered_feature_bc_matrix" \
             "${CELL_CALLING_DIR}/1E_ARC/P0_Multiome_2/outs/filtered_feature_bc_matrix" )
SC_EXP_OBJ=( "${GENE_EXPR_DIR}/1B_GEX/treated_filt/hvg_pca_clust/object.Rds" \
             "${GENE_EXPR_DIR}/1D_GEX/treated_filt/hvg_pca_clust/object.Rds")
SC_SAMPLES=( "T13,T17" "T13_2_Multiseq,T15_2_Multiseq" )

LIN1_MARKERS="${GENE_EXPR_DIR}/signatures/DEG_treated/signature_lineage1.txt"
LIN2_MARKERS="${GENE_EXPR_DIR}/signatures/DEG_treated/signature_lineage2.txt"
             
CL_INFO=( "${MULTIOME_DATA_DIR}/1C_ARC_pca_rna_0.4_clusters.csv" \
          "${MULTIOME_DATA_DIR}/1E_ARC_pca_rna_0.6_clusters.csv" )
          
MEM="10000"

ID_PREP="prepr"
bsub -J ${ID_PREP} -o ${ID_PREP}.STDOUT -e ${ID_PREP}.STDERR bash COMMAND_CNV_intersect.sh

for (( i=0 ; i<${#SC_EXP_NAME_ARC[@]} ; i++ ))
do

    ID="cnv_atac_${SC_EXP_NAME_ARC[$i]}"
    ARGS=( "${WES_EXP_FILE}" "${SC_EXP_NAME_ARC[$i]}" "${SC_EXP_DIR[$i]}" "${CL_INFO[$i]}" )
    COMMAND="bash ${SCRIPT_DIR}/granges_comparison_scATAC.sh ${ARGS[@]}"
    bsub -J ${ID} -w "done(${ID_PREP})" -o ${ID}.STDOUT -e ${ID}.STDERR ${COMMAND}
    
    ID="cnv_rna_${SC_EXP_NAME_GEX[$i]}"
    ARGS=( "${WES_EXP_FILE}" "${SC_EXP_NAME_GEX[$i]}" "${SC_EXP_OBJ[$i]}" "${SC_SAMPLES[$i]}" "${SC_EXP_DIR[$i]}/features.tsv.gz" \
           "${LIN1_MARKERS}" "${LIN2_MARKERS}" )
    COMMAND="bash ${SCRIPT_DIR}/granges_comparison_scRNA.sh ${ARGS[@]}"
    bsub -M ${MEM} -J ${ID} -w "done(${ID_PREP})" -o ${ID}.STDOUT -e ${ID}.STDERR ${COMMAND}

done


exit

