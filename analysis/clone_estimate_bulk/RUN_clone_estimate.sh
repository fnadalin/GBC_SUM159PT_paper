module purge
module load r-4.1.0-gcc-9.3.0-wvnko7v gmp-6.1.2-gcc-9.3.0-hicntdj
R_LIBS_USER="/hps/software/users/marioni/francesca/R_libs"
export R_LIBS_USER

SCRIPT_DIR=$( cd "../../scripts/clone_estimate/" ; pwd )

FILE_LIST_INVITRO=( "files_exp_1_GBC_in_vitro.tsv" \
                    "files_exp_4_GBC_in_vitro.tsv" \
                    "files_exp_5_GBC_in_vitro.tsv" )
FILE_LIST_INVIVO=( "files_exp_1_GBC_in_vivo.tsv" \
                     "files_exp_1_GBC_in_vivo_treated.tsv" )
FILE_LIST=( ${FILE_LIST_INVITRO[@]} ${FILE_LIST_INVIVO[@]} )

OUTDIR_INVITRO=( "exp_1_GBC_in_vitro" \
                 "exp_4_GBC_in_vitro" \
                 "exp_5_GBC_in_vitro" )
OUTDIR_INVIVO=( "exp_1_GBC_in_vivo" \
                "exp_1_GBC_in_vivo_treated" )
OUTDIR=( ${OUTDIR_INVITRO[@]} ${OUTDIR_INVIVO[@]} )

for (( i = 0 ; i < ${#OUTDIR[@]} ; i++ ))
do
    [[ -d ${OUTDIR[$i]} ]] || mkdir -p ${OUTDIR[$i]}

    Rscript ${SCRIPT_DIR}/collect_gbc.R ${FILE_LIST[$i]} ${OUTDIR[$i]}/counts.tsv
    Rscript ${SCRIPT_DIR}/normalize.R ${OUTDIR[$i]}/counts.tsv ${OUTDIR[$i]}/cpm.tsv

done

for (( i = 0 ; i < ${#OUTDIR_INVITRO[@]} ; i++ ))
do
    in="${OUTDIR_INVITRO[$i]}"
    out="${OUTDIR_INVITRO[$i]}/stat"
    out_passive="${OUTDIR_INVITRO[$i]}/stat_passive"
    files="${FILE_LIST_INVITRO[$i]}"
    
    [[ -d "${out}" ]] || mkdir -p "${out}"
    [[ -d "${out_passive}" ]] || mkdir -p "${out_passive}"
    
    # plot cumulative distribution
    bash ${SCRIPT_DIR}/plot_cumulative_cpm.sh ${files} ${in} ${out} "early" "late"
    bash ${SCRIPT_DIR}/plot_cumulative_cpm.sh ${files} ${in} ${out_passive} "early" "passive"

    # compute stats on GBCs
    Rscript ${SCRIPT_DIR}/expansion_rate.R ${files} ${in}/cpm.tsv ${out}/log2fc.tsv "early" "late"
    Rscript ${SCRIPT_DIR}/expansion_rate.R ${files} ${in}/cpm.tsv ${out_passive}/log2fc.tsv "early" "passive"
    Rscript ${SCRIPT_DIR}/ma_stats.R ${files} ${in}/cpm.tsv ${out}/ma.tsv "early" "late"
    Rscript ${SCRIPT_DIR}/ma_stats.R ${files} ${in}/cpm.tsv ${out_passive}/ma.tsv "early" "passive"
    
    # select GBCs
    Rscript ${SCRIPT_DIR}/gbc_select.R ${in}/cpm.tsv ${out}
    Rscript ${SCRIPT_DIR}/gbc_select.R ${in}/cpm.tsv ${out_passive}
    
    # compute significance
    size=$(grep -c . "${in}/cpm.tsv")
    Rscript ${SCRIPT_DIR}/pvalue.R ${files} ${out} "late" ${size} ${out}/pvalues.txt
    Rscript ${SCRIPT_DIR}/pvalue.R ${files} ${out_passive} "passive" ${size} ${out_passive}/pvalues.txt
done

for (( i = 0 ; i < ${#OUTDIR_INVIVO[@]} ; i++ ))
do
    in="${OUTDIR_INVIVO[$i]}"
    out="${OUTDIR_INVIVO[$i]}/stat"
    files="${FILE_LIST_INVIVO[$i]}"
    
    [[ -d "${out}" ]] || mkdir -p "${out}"
    
    # plot cumulative distribution
    bash ${SCRIPT_DIR}/plot_cumulative_cpm.sh ${files} ${in} ${out} "parental" "tumor"
    bash ${SCRIPT_DIR}/plot_cumulative_cpm.sh ${files} ${in} ${out_passive} "parental" "tumor"

    # compute stats on GBCs
    Rscript ${SCRIPT_DIR}/expansion_rate.R ${files} ${in}/cpm.tsv ${out}/log2fc.tsv "parental" "tumor"
    Rscript ${SCRIPT_DIR}/ma_stats.R ${files} ${in}/cpm.tsv ${out}/ma.tsv "parental" "tumor"
    
    # select GBCs
    Rscript ${SCRIPT_DIR}/gbc_select.R ${in}/cpm.tsv ${out}
    
    # compute significance
    size=$(grep -c . "${in}/cpm.tsv")
    Rscript ${SCRIPT_DIR}/pvalue.R ${files} ${out} "tumor" ${size} ${out}/pvalues.txt
done



