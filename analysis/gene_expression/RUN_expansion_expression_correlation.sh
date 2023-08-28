SCRIPT_DIR=$( cd "../../scripts/exp_correlation" ; pwd )
GIT_DIR=$( cd "../../../git/r_scripts/gene_expression_analysis" ; pwd )

EXP=( "1B_GEX" "1D_GEX" )
CTRL=( "T0_1,T0_2" "T0_1_Multiseq,T0_2_Multiseq" )
CASE=( "T13,T17" "T13_2_Multiseq,T15_2_Multiseq" )
SAMPLE=( "T9" "T7_2_Multiseq,T9_2_Multiseq,T11_2_Multiseq" )

SEEDS="SEEDS.txt"

MEM1="2000"
MEM2="25000"
MEM3="16000"
CORES="4"

i=0
for exp in ${EXP[@]}
do  
    METADATA="${exp}/CB_metadata_filt.tsv"
    OUT="gene_exp_correlation/${exp}"
    OUT_nobackup="${NOBACK}/${OUT}"
    
    OUT_CL="${OUT}/clones"
    [[ -d "${OUT_CL}" ]] || mkdir -p "${OUT_CL}"

    # compute log2FC on cpt
    ID1="clones_${exp}"
    COMMAND="${SCRIPT_DIR}/1_expansion_rate.sh" 
    ARGS=( "${METADATA}" "${OUT_CL}" ${CTRL[$i]} ${CASE[$i]} )
    bsub -J ${ID1} -M ${MEM1} -o ${ID1}.STDOUT -e ${ID1}.STDERR bash ${COMMAND} ${ARGS[@]}

    OBJ="${exp}/all_filt/object.Rds"
    OUT_MTX="${OUT}/gene_clone_matrix"
    [[ -d "${OUT_MTX}" ]] || mkdir -p "${OUT_MTX}"

    # compute clone-gene expression matrix (averaged across cells)
    ID2="mtx_${exp}"
    COMMAND="${SCRIPT_DIR}/2_gene_clone_matrix.sh" 
    ARGS=( "${OBJ}" "${OUT_MTX}" )
    bsub -J ${ID2} -M ${MEM2} -o ${ID2}.STDOUT -e ${ID2}.STDERR bash ${COMMAND} ${ARGS[@]}
    
    IFS=","
    read -r -a samples <<< ${SAMPLE[$i]}
    IFS=$'\n' 

    # compute gene expression correlation
    for s in ${samples[@]}
    do
        ID3="exp_corr_${exp}_$s"
        COMMAND="${SCRIPT_DIR}/3_gene_expression_correlation.sh"
        ARGS=( "${OUT_CL}" "${OUT_MTX}/$s" "${OUT}/$s" )
        bsub -J ${ID3} -n ${CORES} -M ${MEM3} -w "done(${ID1})&&done(${ID2})" -o ${ID3}.STDOUT -e ${ID3}.STDERR bash ${COMMAND} ${ARGS[@]}
        
        ID4="permutest_${exp}_$s"
        COMMAND="${SCRIPT_DIR}/4_permutation_test.sh"
        ARGS=( "${OBJ}" "$s" "${OUT_CL}" "${OUT}/$s" "${OUT}/$s" "${SEEDS}" )
        bsub -J ${ID4} -n ${CORES} -M ${MEM3} -w "done(${ID3})" -o ${ID4}.STDOUT -e ${ID4}.STDERR bash ${COMMAND} ${ARGS[@]}

    done

    i=$((i+1))
done




exit

