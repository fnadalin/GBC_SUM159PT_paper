
if [[ $# -lt "3" ]]
then
    echo "bash $0 <obj> <dir_sc> <pops>"
    exit
fi

OBJ=$1
DIR=$2
POPS=$3

LOCAL_SCRIPT_DIR=$( cd $(dirname $0) ; pwd ) 
R_SCRIPT_POST="${LOCAL_SCRIPT_DIR}/scissor_postprocessing.R"
R_SCRIPT_FISHER="${LOCAL_SCRIPT_DIR}/fisher_test.R"
R_SCRIPT_DEA="${LOCAL_SCRIPT_DIR}/DEA.R"

[[ -d "${DIR}/out/DEA" ]] || mkdir -p "${DIR}/out/DEA"

Rscript ${R_SCRIPT_POST} ${OBJ} "${DIR}/in/sc_dataset.Rds" "${DIR}/out" ${POPS}
Rscript ${R_SCRIPT_FISHER} "${DIR}/in/sc_dataset.Rds" ${POPS} "${DIR}/out/fisher"
Rscript ${R_SCRIPT_DEA} "${DIR}/in/sc_dataset.Rds" ${POPS} "${DIR}/out/DEA"


exit

