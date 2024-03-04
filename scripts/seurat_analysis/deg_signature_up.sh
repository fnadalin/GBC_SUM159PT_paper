
module purge
module load r-4.1.0-gcc-9.3.0-wvnko7v gmp-6.1.2-gcc-9.3.0-hicntdj
R_LIBS_USER="/hps/software/users/marioni/francesca/R_libs"
export R_LIBS_USER

if [[ $# -lt "3" ]]
then
    echo ""
    echo "Usage: bash $0 <dea1> <dea2> <out> <log2fc> <top>"
    echo ""
    exit
fi

SCRIPT_DIR=$( cd $(dirname $0) ; pwd )

Rscript "${SCRIPT_DIR}/deg_signature_up.R" $@

exit


