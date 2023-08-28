
if [[ $# -lt "3" ]]
then
    echo "bash $0 <script_dir> <matrix_list> <out_dir>"
    exit
fi

SCRIPT_DIR="$1"
MATRIX_LIST="$2"
OUT_DIR="$3"

R_SCRIPT_COMPUTE_STATS="${SCRIPT_DIR}/1_compute_matrix_stats.R"

# min read count to define a GBC as detected
MIN_READ_COUNT="2" 
# min UMI count to define a CB as detected
MIN_UMI_COUNT="5000"

echo "compute matrix stats"
IFS=$'\n'
for line in $(cat ${MATRIX_LIST})
do
    SAMPLE_ID=$(echo $line | cut -f 1)
    MATRIX_DIR=$(echo $line | cut -f 2)
    
    echo ${SAMPLE_ID}

    out_dir="${OUT_DIR}/${SAMPLE_ID}"
    [[ -d "${out_dir}" ]] || mkdir -p "${out_dir}"

    Rscript "${R_SCRIPT_COMPUTE_STATS}" \
            --matrix_dir "${MATRIX_DIR}" \
            --out_dir "${out_dir}"  \
            --min_read_count ${MIN_READ_COUNT} \
            --min_UMI_count ${MIN_UMI_COUNT} \
            --plot_title "${SAMPLE_ID}" \
            1> "${out_dir}/1_compute_matrix_stats.STDOUT" \
            2> "${out_dir}/1_compute_matrix_stats.STDERR"
done


exit

