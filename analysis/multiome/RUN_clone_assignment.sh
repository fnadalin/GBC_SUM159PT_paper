SCRIPT_DIR=$( cd "../../scripts/seurat_analysis" ; pwd )
GIT_DIR=$( cd "../../../git/r_scripts/gene_expression_analysis" ; pwd )
IN_DIR=$( cd "../clone_calling_sc" ; pwd )
WDIR=$(pwd)

EXP=( "1C_GEX" "1E_GEX" )
EXP_ARC=( "1C_ARC" "1E_ARC" )

CLONE_INFO="clone_info.tsv"

MEM="2000"

i=0
for exp in ${EXP_ARC[@]}
do
    FILE_LIST="matrix_dir_${exp}.tsv"
    IN="${IN_DIR}/${EXP[$i]}"
    
    # Generate metadata for clone-labelled cells and cross with cellranger-arc results
    [[ -d "${exp}" ]] || mkdir "${exp}"
    COMMAND="${SCRIPT_DIR}/1_create_metadata_arc.sh" 
    ARGS=( "${GIT_DIR}" "${FILE_LIST}" "${IN}" "${exp}" "${CLONE_INFO}" )
    bsub -J step1_${exp} -M ${MEM} -o step1_${exp}.STDOUT -e step1_${exp}.STDERR bash ${COMMAND} ${ARGS[@]}

    i=$((i+1))
done


exit

