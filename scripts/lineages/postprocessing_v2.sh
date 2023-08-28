
TEST="MAST"
CL_MODE="lineage"

if [[ $# -lt "4" ]]
then
    echo "bash $0 <script_dir> <script_dir_1> <file_list> <outdir> <k>"
    exit
fi

SCRIPT_DIR="$1"
SCRIPT_DIR_1="$2"
LINEAGE_INFO="$3"
OUT_DIR="$4"
K="$5"

DIR=$( cd $(dirname $0) ; pwd ) 

R_SCRIPT_DEA="${SCRIPT_DIR}/Seurat_DEA.R"
R_PROP_MATRIX="${SCRIPT_DIR_1}/2_compute_propensity_matrix.R"

# input
FILE_LIST=$( cut -f 3 "${LINEAGE_INFO}" | tr '\n' ',' | sed "s/,$//g" )
EXP_LIST=$( cut -f 1 "${LINEAGE_INFO}" | tr '\n' ',' | sed "s/,$//g" )
EXP=( $( cut -f 1 "${LINEAGE_INFO}" | sort | uniq ) )

# output
LIN="${OUT_DIR}/lineages_merge.tsv"
LIN_DIR="${OUT_DIR}/lineage_merge/"

Rscript "${DIR}/collect_lineages.R" "${FILE_LIST}" "${EXP_LIST}" "${LIN}"

for exp in ${EXP[@]}
do
    SAMPLE=( $( ls "${OUT_DIR}/${exp}/graphs/" ) )
    
    for sample in ${SAMPLE[@]} 
    do
        OUTDIR_G="${OUT_DIR}/${exp}/graphs/${sample}"
        LIN_DIR_SAMPLE="${LIN_DIR}/${sample}"
        [[ -d ${LIN_DIR_SAMPLE} ]] || mkdir -p ${LIN_DIR_SAMPLE}
        
        Rscript "${R_PROP_MATRIX}" "${LIN}" "${OUTDIR_G}" "${LIN_DIR_SAMPLE}" $K
    done
    
done

for exp in ${EXP[@]}
do
    IN_DIR="${exp}/treated_filt"

    HPC_DIR="${IN_DIR}/hvg_pca_clust"
    DEA_DIR="${IN_DIR}/DEA/${CL_MODE}"
    [[ -d "${DEA_DIR}" ]] || mkdir -p "${DEA_DIR}"

    IN_OBJ="${HPC_DIR}/object.Rds"
    
    Rscript "${DIR}/add_lineage_info.R" "${IN_OBJ}" "${LIN}"

    # run differential expression analysis between lineages
    Rscript "${R_SCRIPT_DEA}" \
                --object "${IN_OBJ}" \
                --out_dir "${DEA_DIR}" \
                --test "${TEST}" \
                --cl_mode "${CL_MODE}" \
                --latent_vars "sample.name"
done


exit

