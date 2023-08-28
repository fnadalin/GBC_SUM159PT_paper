
SCRIPT_DIR=$( cd "../../scripts/sister_cells" ; pwd )

EXP=( "1B_GEX" "1D_GEX" )
CL1="2,5,6"
CL2="2,4,6"

ID="ident_${EXP[0]}-${EXP[1]}"
DIR="${EXP[0]}-${EXP[1]}_ccRegress"
COMMAND="${SCRIPT_DIR}/transcriptionally_identifiable_clones.sh"
ARGS=( "${DIR}" "${CL1}" "${CL2}" "${DIR}/transcriptionally_identifiable_clones.txt" )
bsub -J ${ID} -o ${ID}.STDOUT -e ${ID}.STDERR bash ${COMMAND} ${ARGS[@]}

exit
