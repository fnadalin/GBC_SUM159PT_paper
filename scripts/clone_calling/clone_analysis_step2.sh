
if [[ $# -lt "3" ]]
then
    echo "bash $0 <script_dir> <params> <out_dir>"
    exit
fi

SCRIPT_DIR="$1"
PARAMS="$2"
OUT_DIR="$3"

R_SCRIPT_CLASSIFY_CB="${SCRIPT_DIR}/2_classify_CB.R"
R_SCRIPT_EXTRACT_CLONES="${SCRIPT_DIR}/3_extract_clones.R"
R_SCRIPT_CROSS_STAT="${SCRIPT_DIR}/4_cross_sample_stats.R"
R_SCRIPT_SURV_STAT="${SCRIPT_DIR}/5_surv_clones_stats.R"

source ./${PARAMS}

# files containing the CB classification
RESULTS_DIR_LIST="${OUT_DIR}/results_dir_list.tsv"

# classify the cells based on GBC read counts
echo "classify the cells based on GBC read counts"
i=0
for SAMPLE in ${samples[@]}
do
    echo "${SAMPLE}"
    in_dir="${OUT_DIR}/${SAMPLE}"
    out_dir="${in_dir}/classification"
    [[ -d "${out_dir}" ]] || mkdir -p "${out_dir}"

    Rscript "${R_SCRIPT_CLASSIFY_CB}" \
            --in_dir "${in_dir}" \
            --out_dir "${out_dir}" \
            --doublet_rate ${doublet_rate[$i]} \
            --min_umi_count ${MIN_UMI_COUNT} \
            --max_umi_count ${max_umi_count[$i]} \
            --min_read_count_expr ${MIN_READ_COUNT_EXPR} \
            --min_read_frac_expr ${MIN_READ_FRAC_EXPR} \
            --min_rank1_read_frac_expr ${MIN_RANK1_READ_FRAC_EXPR} \
            --pval_cutoff ${pval_cutoff[$i]} \
            --doublet_prob_cutoff ${DOUBLET_PROB_CUTOFF} \
            --plot_title "${SAMPLE}" \
            1> "${out_dir}/2_classify_CB.STDOUT" 2> "${out_dir}/2_classify_CB.STDERR"

    Rscript "${R_SCRIPT_EXTRACT_CLONES}" \
            --in_dir "${out_dir}" \
            --out_dir "${out_dir}" \
            --plot_title "${SAMPLE}" \
            1> "${out_dir}/3_extract_clones.STDOUT" 2> "${out_dir}/3_extract_clones.STDERR"

    i=$((i+1))
done

# prepare the input tsv file
echo "prepare the input tsv file"
i=0
for SAMPLE in ${samples[@]}
do
    RESULTS_DIR="${OUT_DIR}/${SAMPLE}/classification"
    SURV=${surv[$i]}

    echo "${SAMPLE}	${RESULTS_DIR}	${SURV}"

    i=$((i+1))
done > "${RESULTS_DIR_LIST}"

# collect the results across samples
echo "collect the results across samples"
out_dir="${OUT_DIR}/merged_samples"
[[ -d ${out_dir} ]] || mkdir -p ${out_dir}
Rscript "${R_SCRIPT_CROSS_STAT}" \
        --results_dir_list "${RESULTS_DIR_LIST}" \
        --out_dir ${out_dir} \
        1> "${out_dir}/4_cross_sample_stats.STDOUT" \
        2> "${out_dir}/4_cross_sample_stats.STDERR"

# collect the results on surviving clones
echo "collect the results on surviving clones"
Rscript "${R_SCRIPT_SURV_STAT}" \
        --results_dir_list "${RESULTS_DIR_LIST}" \
        --out_dir ${out_dir} \
        1> "${out_dir}/5_surv_clones_stats.STDOUT" \
        2> "${out_dir}/5_surv_clones_stats.STDERR"


exit

