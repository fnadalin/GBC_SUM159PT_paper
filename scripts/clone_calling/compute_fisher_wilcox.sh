

if [[ $# -lt "3" ]]
then
    echo "bash $0 <files> <out_dir> <file_list>"
    exit
fi

FILES="$1"
OUT_DIR="$2"
FILE_LIST="$2"

SCRIPT_DIR=$( cd $(dirname $0) ; pwd )

[[ -d "${OUT_DIR}/early_passive_fisher" ]] || mkdir -p "${OUT_DIR}/early_passive_fisher"
[[ -d "${OUT_DIR}/early_passive_wilcox" ]] || mkdir -p "${OUT_DIR}/early_passive_wilcox"

Rscript ${SCRIPT_DIR}/collect_clones_sc.R ${FILES} ${OUT_DIR}
Rscript ${SCRIPT_DIR}/define_DRC_TIC_fisher_sc.R ${OUT_DIR}/cpm.tsv ${FILE_LIST} "early" "passive" "${OUT_DIR}/early_passive_fisher"
Rscript ${SCRIPT_DIR}/define_DRC_TIC_wilcox_sc.R ${OUT_DIR}/cpm.tsv ${FILE_LIST} "early" "passive" "${OUT_DIR}/early_passive_wilcox"

exit


