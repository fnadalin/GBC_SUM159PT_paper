CMD_DIR=$( cd "../../scripts/GBC_calling" ; pwd )
SCRIPT_DIR=$( cd "$HOME/git/barcode_groups" ; pwd )

echo "SCRIPT_DIR=\"${SCRIPT_DIR}\""


#############################################################################
############################### COLLECT GBC SEQ #############################
#############################################################################

echo ""
echo "======================= COLLECT GBC SEQ ======================="
echo ""

FILES_DNA=( "files_counts_exp_1B_DNA_amplicon_v2.tsv" \
            "files_counts_exp_1C_DNA_amplicon_v2.tsv" \
            "files_counts_exp_1_invivo_A_amplicon_v2.tsv" \
            "files_counts_exp_1_invivo_B_amplicon_v2.tsv" \
            "files_counts_exp_1_invivo_C_amplicon_v2.tsv" \
            "files_counts_exp_1_invivo_D_amplicon_v2.tsv" )
FILES_RNA=( "files_counts_exp_1B_RNA_amplicon_v2.tsv" \
            "files_counts_exp_1C_RNA_amplicon_v2.tsv" \
            "files_counts_exp_1D_RNA_amplicon_v2.tsv" \
            "files_counts_exp_1E_RNA_amplicon_v2.tsv" )
FILES=( ${FILES_DNA[@]} ${FILES_RNA[@]} )

DIR="exp_1_amplicon_v2"
OUT_DNA=( "1B_DNA" "1C_DNA" "1_invivo_A" "1_invivo_B" "1_invivo_C" "1_invivo_D" )
OUT_RNA=( "1B_RNA" "1C_RNA" "1D_RNA" "1E_RNA" )
OUT=( ${OUT_DNA[@]} ${OUT_RNA[@]} )

# Collect DNA barcodes
for (( i=0; i<${#FILES_DNA[@]}; i++ ))
do
    [[ -d "${DIR}"/"${OUT_DNA[$i]}" ]] || mkdir -p "${DIR}"/"${OUT_DNA[$i]}"
    bash "${CMD_DIR}"/collect_DNA_seq_amplicon.sh "${FILES_DNA[$i]}" "${DIR}/${OUT_DNA[$i]}"
done

# Collect RNA barcodes
N_GBC_RNA=${#FILES_RNA[@]}
for (( i=0; i<${#FILES_RNA[@]}; i++ ))
do
    [[ -d "${DIR}"/"${OUT_RNA[$i]}" ]] || mkdir -p "${DIR}"/"${OUT_RNA[$i]}"
    bash "${CMD_DIR}"/collect_RNA_seq_amplicon.sh "${FILES_RNA[$i]}" "${DIR}"/"${OUT_RNA[$i]}"
done

# Name barcodes with a unique GBC id
[[ -d "${DIR}" ]] || mkdir -p "${DIR}"
FILES_LIST=$(IFS=','; echo "${FILES[*]}")
OUT_LIST=$(IFS=','; echo "${OUT[*]}")
bash "${CMD_DIR}"/assign_GBC_id.sh "${FILES_LIST}" "${OUT_LIST}" "${DIR}" 

echo "DONE!"

#############################################################################
############################## CREATE GBC GROUPS ############################
#############################################################################

echo ""
echo "====================== CREATE GBC GROUPS ======================"
echo ""

DIST="1"
FRAC="1"

# assign read count to barcodes
for (( i=0; i < ${#FILES[@]}; i++ ))
do
    bash "${CMD_DIR}"/group_GBC.sh "${SCRIPT_DIR}" "${FILES[$i]}" "${DIR}"/"${OUT[$i]}" ${DIST} ${FRAC}
done

INFILE="${DIR}/GBC.txt"

# collect GBC IDs and sequences across samples
for (( i=0; i < ${#FILES_RNA[@]}; i++ ))
do
    OUTFILE="${DIR}/GBC_${OUT_RNA[$i]}_d${DIST}_f${FRAC}.tsv"
    bash "${CMD_DIR}"/collect_hubs_for_CR.sh "${FILES_RNA[$i]}" "${DIR}"/"${OUT_RNA[$i]}" ${DIST} ${FRAC} "${INFILE}" "${OUTFILE}"
done

echo "DONE!"






exit


