module purge
module load r-4.1.0-gcc-9.3.0-wvnko7v gmp-6.1.2-gcc-9.3.0-hicntdj
R_LIBS_USER="/hps/software/users/marioni/francesca/R_libs"
export R_LIBS_USER

if [[ $# -lt "3" ]]
then
    echo "bash $0 <dir1> <dir2> <out_dir>"
    exit
fi

DIR1=$1
DIR2=$2
OUT_DIR=$3

LOCAL_SCRIPT_DIR=$( cd $(dirname $0) ; pwd ) 
R_SCRIPT="${LOCAL_SCRIPT_DIR}/heatmap_auc_merged.R"

[[ -d "${OUT_DIR}" ]] || mkdir -p "${OUT_DIR}"

Rscript ${R_SCRIPT} ${DIR1} ${DIR2} ${OUT_DIR} 

exit

