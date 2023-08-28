
CMD_DIR=$( cd "../../scripts/GBC_calling" ; pwd )
SCRIPT_DIR=$( cd "$HOME/git/barcode_groups" ; pwd )

echo "SCRIPT_DIR=\"${SCRIPT_DIR}\""


#############################################################################
############################### COLLECT GBC SEQ #############################
#############################################################################

echo ""
echo "======================= COLLECT GBC SEQ ======================="
echo ""

FILES_DNA=( "files_counts_exp_5A_DNA_amplicon_v2.tsv" )
FILES=( ${FILES_DNA[@]} )

DIR="exp_5_amplicon_v2"
OUT_DNA=( "5A_DNA" )
OUT=( ${OUT_DNA[@]} )

# Collect DNA barcodes
for (( i=0; i<${#FILES_DNA[@]}; i++ ))
do
    [[ -d "${DIR}"/"${OUT_DNA[$i]}" ]] || mkdir -p "${DIR}"/"${OUT_DNA[$i]}"
    bash "${CMD_DIR}"/collect_DNA_seq_amplicon.sh "${FILES_DNA[$i]}" "${DIR}/${OUT_DNA[$i]}"
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

echo "DONE!"






exit




