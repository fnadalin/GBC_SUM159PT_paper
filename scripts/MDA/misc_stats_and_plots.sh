module purge
module load r-4.1.0-gcc-9.3.0-wvnko7v gmp-6.1.2-gcc-9.3.0-hicntdj
R_LIBS_USER="/hps/software/users/marioni/francesca/R_libs"
export R_LIBS_USER

SCRIPT_DIR=$( cd $(dirname $0) ; pwd )

SIGN_DIR=$( cd $(pwd) ; cd "../../analysis/gene_expression/signatures/DEG_P0_ccRegress" ; pwd )
P0="19924_ccRegress"

if [[ $# -lt 2 ]]
then
    echo "Usage: bash $0 <dir> <slot>"
fi

DIR=$1
SLOT=$2

for i in 1 2 3
do
    SIGN_FILE="${SIGN_DIR}/signature${i}.txt"
    Rscript ${SCRIPT_DIR}/misc_stats_and_plots.R "${DIR}" "${SLOT}" "${SIGN_FILE}" "S${i}"
done

exit

