
module purge
module load r-4.1.0-gcc-9.3.0-wvnko7v gmp-6.1.2-gcc-9.3.0-hicntdj
R_LIBS_USER="/hps/software/users/marioni/francesca/R_libs"
export R_LIBS_USER

SCRIPT_DIR=$( cd $(dirname $0) ; pwd )

P0="19924_ccRegress"

if [[ $# -lt 4 ]]
then
    echo "Usage: bash $0 <script_dir> <dir> <modeP0> <clP0>"
fi

SCRIPT_DIR=$1
dir=$2
modeP0=$3
clP0=$4

dir_P0_filt="${dir}/${P0}_filt"
[[ -d "${dir_P0_filt}" ]] || mkdir -p "${dir_P0_filt}"
cellsP0="${dir_P0_filt}/cells.txt"

# extract the cells
Rscript ${SCRIPT_DIR}/cell_subset.R ${dir}/${P0}/hvg_pca_clust/object.Rds ${modeP0} ${clP0} ${cellsP0}

# subset the object
Rscript ${SCRIPT_DIR}/obj_subset.R ${dir}/${P0}/object.Rds ${cellsP0} ${dir_P0_filt}/object.Rds

exit

