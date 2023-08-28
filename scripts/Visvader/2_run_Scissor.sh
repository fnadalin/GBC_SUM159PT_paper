
if [[ $# -lt "10" ]]
then
    echo "bash $0 <obj1> <obj2> <mode1> <mode2> <cl1> <cl2> <ID1> <ID2> <out_dir> <id>"
    exit
fi

OBJ1=$1
OBJ2=$2
MODE1=$3
MODE2=$4
CL1=$5
CL2=$6
ID1=$7
ID2=$8
OUT_DIR=$9
ID=${10}

LOCAL_SCRIPT_DIR=$( cd $(dirname $0) ; pwd ) 

R_SCRIPT_PSEUDOBULK="${LOCAL_SCRIPT_DIR}/extract_pseudobulk_phenotypes.R"
R_SCRIPT_MERGE="${LOCAL_SCRIPT_DIR}/merge_pseudobulk_phenotypes.R"
R_SCRIPT_SCISSOR="${LOCAL_SCRIPT_DIR}/scissor.R"

DIR_IN="${OUT_DIR}/in/${ID}"
DIR_OUT="${OUT_DIR}/out/${ID}"
[[ -d "${DIR_IN}" ]] || mkdir -p "${DIR_IN}"
[[ -d "${DIR_OUT}" ]] || mkdir -p "${DIR_OUT}"

Rscript ${R_SCRIPT_PSEUDOBULK} ${OBJ1} ${MODE1} ${CL1} ${DIR_IN}/${ID1}
Rscript ${R_SCRIPT_PSEUDOBULK} ${OBJ2} ${MODE2} ${CL2} ${DIR_IN}/${ID2}
Rscript ${R_SCRIPT_MERGE} ${DIR_IN}/${ID1},${DIR_IN}/${ID2} ${DIR_IN}/merge

Rscript ${R_SCRIPT_SCISSOR} "${OUT_DIR}/in/sc_dataset.Robj" ${DIR_IN}/merge ${ID} ${DIR_OUT}
        

exit

