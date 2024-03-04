
module purge
module load r-4.1.0-gcc-9.3.0-wvnko7v gmp-6.1.2-gcc-9.3.0-hicntdj
R_LIBS_USER="/hps/software/users/marioni/francesca/R_libs"
export R_LIBS_USER

QVAL="1"

if [[ $# -lt "3" ]]
then
    echo ""
    echo "Usage: bash $0 <r_script_dir> <gene_list> <out_prefix>"
    echo ""
    exit
fi

R_SCRIPT_DIR=$1
GENE_LIST=$2
OUT_PREFIX=$3

Rscript "${R_SCRIPT_DIR}/functional_annotation.R" --genes ${GENE_LIST} --out_prefix ${OUT_PREFIX} --ann_qval ${QVAL}

exit


