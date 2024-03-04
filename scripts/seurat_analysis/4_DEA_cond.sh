
module purge
module load r-4.1.0-gcc-9.3.0-wvnko7v gmp-6.1.2-gcc-9.3.0-hicntdj
R_LIBS_USER="/hps/software/users/marioni/francesca/R_libs"
export R_LIBS_USER

TEST="MAST"

if [[ $# -lt "3" ]]
then
    echo ""
    echo "Usage: bash $0 <script_dir> <OUT_DIR> <cond_modes>"
    echo ""
    exit
fi

SCRIPT_DIR="$1"
OUT_DIR="$2"
COND_MODES="$3"

R_SCRIPT_DEA_SUB="${SCRIPT_DIR}/Seurat_DEA_subset.R"

HPC_DIR="${OUT_DIR}/hvg_pca_clust"
IN_OBJ="${HPC_DIR}/object.Rds"

IFS=$','
read -r -a conds <<< ${COND_MODES}     

# run differential expression analysis between conditions within samples
echo "run DEA between conditions within samples"
for cond in ${conds[@]}
do
    DEA_DIR="${OUT_DIR}/DEA/${cond}"
    [[ -d "${DEA_DIR}" ]] || mkdir -p "${DEA_DIR}"
    Rscript "${R_SCRIPT_DEA_SUB}" \
            --object "${IN_OBJ}" \
            --out_dir "${DEA_DIR}" \
            --test "${TEST}" \
            --cond_mode ${cond} \
            --subset "sample.name"
done

exit


