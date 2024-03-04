
module purge
module load r-4.1.0-gcc-9.3.0-wvnko7v gmp-6.1.2-gcc-9.3.0-hicntdj
R_LIBS_USER="/hps/software/users/marioni/francesca/R_libs"
export R_LIBS_USER

TEST="MAST"

if [[ $# -lt "3" ]]
then
    echo ""
    echo "Usage: bash $0 <script_dir> <OUT_DIR> <CL_MODE>"
    echo ""
    exit
fi

SCRIPT_DIR="$1"
OUT_DIR="$2"
CL_MODE="$3"

R_SCRIPT_DEA="${SCRIPT_DIR}/Seurat_DEA.R"
R_SCRIPT_DEA_SUB="${SCRIPT_DIR}/Seurat_DEA_subset.R"

HPC_DIR="${OUT_DIR}/hvg_pca_clust"
DEA_DIR="${OUT_DIR}/DEA/${CL_MODE}"
[[ -d "${DEA_DIR}" ]] || mkdir -p "${DEA_DIR}"

IN_OBJ="${HPC_DIR}/object.Rds"

# run differential expression analysis between clusters
echo "run DEA between clusters"
Rscript "${R_SCRIPT_DEA}" \
            --object "${IN_OBJ}" \
            --out_dir "${DEA_DIR}" \
            --test "${TEST}" \
            --cl_mode "${CL_MODE}"
            
exit


