
module purge
module load r-4.1.0-gcc-9.3.0-wvnko7v gmp-6.1.2-gcc-9.3.0-hicntdj
R_LIBS_USER="/hps/software/users/marioni/francesca/R_libs"
export R_LIBS_USER

SCRIPT_DIR=$( cd $(dirname $0) ; pwd )

if [[ $# -lt 5 ]]
then
    echo "Usage: bash $0 <dir> <modeP0> <modeT> <clP0> <clT>"
fi

dir=$1
modeP0=$2
modeT=$3
clP0=$4
clT=$5

dir_P0_filt="${dir}/P0_filt"
[[ -d "${dir_P0_filt}" ]] || mkdir -p "${dir_P0_filt}"
cellsP0="${dir_P0_filt}/cells.txt"

dir_T_filt="${dir}/treated_filt"
[[ -d "${dir_T_filt}" ]] || mkdir -p "${dir_T_filt}"
cellsT="${dir_T_filt}/cells.txt"

dir_all_filt="${dir}/all_filt"
[[ -d "${dir_all_filt}" ]] || mkdir -p "${dir_all_filt}"
cellsAll="${dir_all_filt}/cells.txt"

# extract the cells
Rscript ${SCRIPT_DIR}/cell_subset.R ${dir}/P0/hvg_pca_clust/object.Rds ${modeP0} ${clP0} ${cellsP0}
Rscript ${SCRIPT_DIR}/cell_subset.R ${dir}/treated/hvg_pca_clust/object.Rds ${modeT} ${clT} ${cellsT}
cat ${cellsP0} ${cellsT} > ${cellsAll}

# subset the object
Rscript ${SCRIPT_DIR}/obj_subset.R ${dir}/P0/object.Rds ${cellsP0} ${dir_P0_filt}/object.Rds
Rscript ${SCRIPT_DIR}/obj_subset.R ${dir}/treated/object.Rds ${cellsT} ${dir_T_filt}/object.Rds
Rscript ${SCRIPT_DIR}/obj_subset.R ${dir}/all/object.Rds ${cellsAll} ${dir_all_filt}/object.Rds

# subset the metadata
Rscript ${SCRIPT_DIR}/meta_subset.R ${dir}/CB_metadata.tsv ${cellsAll} ${dir}/CB_metadata_filt.tsv

exit

