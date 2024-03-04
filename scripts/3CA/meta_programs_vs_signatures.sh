module purge
module load r-4.1.0-gcc-9.3.0-wvnko7v gmp-6.1.2-gcc-9.3.0-hicntdj
R_LIBS_USER="/hps/software/users/marioni/francesca/R_libs"
export R_LIBS_USER

if [[ $# -lt "7" ]]
then
    echo "bash $0 <table3CA> <object> <cl_mode> <subpop_id> <subpop_name> <out_dir> <sign_dir>"
    exit
fi

TABLE_3CA=$1
OBJECT=$2
CL_MODE=$3
SUBPOP_ID=$4
SUBPOP_NAME=$5
OUT_DIR=$6
SIGN_DIR=$7

LOCAL_SCRIPT_DIR=$( cd $(dirname $0) ; pwd ) 
R_SCRIPT_PREPR="${LOCAL_SCRIPT_DIR}/preprocessing.R"
R_SCRIPT="${LOCAL_SCRIPT_DIR}/meta_programs_vs_signatures_general.R"

[[ -d "${OUT_DIR}" ]] || mkdir -p "${OUT_DIR}"

Rscript ${R_SCRIPT_PREPR} ${TABLE_3CA} ${OBJECT} ${OUT_DIR}
Rscript ${R_SCRIPT} "${OUT_DIR}/MP_table.tsv" ${OBJECT} ${CL_MODE} ${SUBPOP_ID} ${SUBPOP_NAME} ${OUT_DIR} ${SIGN_DIR}

exit

