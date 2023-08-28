SCRIPT_DIR=$( cd "../../scripts/sister_cells" ; pwd )
CLONE_DIR=$( cd "../clone_estimate_bulk" ; pwd )
SC_DIR=$( cd "../gene_expression" ; pwd )

EXP_SC=( "1B_GEX" "1D_GEX" )
EXP_BULK_P0=( "exp_1_GBC_in_vitro" "exp_1_GBC_in_vivo" )
EXP_BULK_DRC=( "exp_1_GBC_in_vitro" )
PHEN=( "DRC_bulk" "DRC_bulk_wilcox" "DRC_bulk_fisher" "DRC_bulk_gam" "TIC" "TIC_wilcox" "TIC_fisher" "TIC_gam" "DRC_sc" )

OUT_DIR_COMP="cmp_cellcount_comparison"
OUT_DIR_PLOTS="plot_TIC_DRC_P0"

PARAMS="param_hvg_pca.txt"

MEM1="2000"
MEM2="20000"

for exp_sc in ${EXP_SC[@]}
do
    FILES_SC="${SC_DIR}/matrix_dir_${exp_sc}_P0.tsv"
    S_SC=$( cut -f 1 "${FILES_SC}" | tr '\n' ',' | sed "s/,$//g" )

    for exp_bulk in ${EXP_BULK_P0[@]}
    do
        OUT="${OUT_DIR_COMP}/P0/${exp_sc}-${exp_bulk}"
        [[ -d "${OUT}" ]] || mkdir -p "${OUT}"
    
        FILES_BULK="${CLONE_DIR}/files_${exp_bulk}.tsv"
        S_BULK=$( awk -F "	" -v e="early" -v p="parental" '{ if ($3 == e || $3 == p) print $1 }' "${FILES_BULK}" | tr '\n' ',' | sed "s/,$//g" )
        
        # Comparison between bulk & sc called clones
        ID="cpm_cellcount_${exp_sc}_${exp_bulk}"
        COMMAND="${SCRIPT_DIR}/compare_cpm_cellcount.sh" 
        ARGS=( "${CLONE_DIR}/${exp_bulk}/cpm.tsv" "${SC_DIR}/${exp_sc}/CB_metadata_filt.tsv" ${S_BULK} ${S_SC} "${OUT}" )
        bsub -J ${ID} -M ${MEM1} -o ${ID}.STDOUT -e ${ID}.STDERR bash ${COMMAND} ${ARGS[@]}

    done
    
    FILES_SC="${SC_DIR}/matrix_dir_${exp_sc}_treated.tsv"
    S_SC=$( cut -f 1 "${FILES_SC}" | tr '\n' ',' | sed "s/,$//g" )
    
    for exp_bulk in ${EXP_BULK_DRC[@]}
    do
        OUT="${OUT_DIR_COMP}/treated/${exp_sc}-${exp_bulk}"
        [[ -d "${OUT}" ]] || mkdir -p "${OUT}"
    
        FILES_BULK="${CLONE_DIR}/files_${exp_bulk}.tsv"
        S_BULK=$( awk -F "	" -v m="middle" -v l="late" '{ if ($3 == m || $3 == l) print $1 }' "${FILES_BULK}" | tr '\n' ',' | sed "s/,$//g" )
        
        # Comparison between bulk & sc called clones
        ID="cpm_cellcount_${exp_sc}_${exp_bulk}"
        COMMAND="${SCRIPT_DIR}/compare_cpm_cellcount.sh" 
        ARGS=( "${CLONE_DIR}/${exp_bulk}/cpm.tsv" "${SC_DIR}/${exp_sc}/CB_metadata_filt.tsv" ${S_BULK} ${S_SC} "${OUT}" )
        bsub -J ${ID} -M ${MEM1} -o ${ID}.STDOUT -e ${ID}.STDERR bash ${COMMAND} ${ARGS[@]}

    done
    
    # plot UMAPs
    INDIR="${SC_DIR}/${exp_sc}/P0_filt/hvg_pca_clust"
    for phen in ${PHEN[@]}
    do
        ID="plot_${exp_sc}_${phen}"
        COMMAND="${SCRIPT_DIR}/plot_umap_TIC_DRC.sh"
        ARGS=( "${INDIR}" "${OUT_DIR_PLOTS}/${exp_sc}" ${phen} "${PARAMS}" )
        bsub -J ${ID} -M ${MEM2} -o ${ID}.STDOUT -e ${ID}.STDERR bash ${COMMAND} ${ARGS[@]}
        
    done
    
done





exit

