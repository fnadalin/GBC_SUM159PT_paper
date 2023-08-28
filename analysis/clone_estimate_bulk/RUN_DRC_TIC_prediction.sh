module purge
module load r-4.1.0-gcc-9.3.0-wvnko7v gmp-6.1.2-gcc-9.3.0-hicntdj
R_LIBS_USER="/hps/software/users/marioni/francesca/R_libs"
export R_LIBS_USER

SCRIPT_DIR=$( cd "../../scripts/clone_estimate/" ; pwd )

FILE_LIST_INVITRO=( "files_exp_1_GBC_in_vitro.tsv" \
                   "files_exp_4_GBC_in_vitro.tsv" )
FILE_LIST_INVIVO=( "files_exp_1_GBC_in_vivo.tsv" \
                   "files_exp_1_GBC_in_vivo_treated.tsv" )
FILE_LIST=( ${FILE_LIST_INVITRO[@]} ${FILE_LIST_INVIVO[@]} )

OUTDIR_INVITRO=( "exp_1_GBC_in_vitro" \
                 "exp_4_GBC_in_vitro" )
OUTDIR_INVIVO=( "exp_1_GBC_in_vivo" \
                "exp_1_GBC_in_vivo_treated" )
OUTDIR=( ${OUTDIR_INVITRO[@]} ${OUTDIR_INVIVO[@]} )

TYPE_INVITRO=( "DRC_bulk" "DRC_bulk" )
TYPE_INVIVO=( "TIC" "TIC" )
TYPE=( ${TYPE_INVITRO[@]} ${TYPE_INVIVO[@]} )

for (( i = 0 ; i < ${#OUTDIR[@]} ; i++ ))
do
    in="${OUTDIR[$i]}"
    stat="${in}/stat"

    out="${in}/${TYPE[$i]}"
    [[ -d "${out}" ]] || mkdir -p "${out}"

    # compute gain
    Rscript ${SCRIPT_DIR}/compute_cpm_gain.R "${stat}" "${out}"
    
    # select DRCs / TICs
    Rscript ${SCRIPT_DIR}/define_DRC_TIC.R "${stat}" "${out}"
    Rscript ${SCRIPT_DIR}/auc_table.R ${in}/cpm.tsv ${out}/list.txt ${out}/auc.txt

done

for (( i = 0 ; i < ${#OUTDIR_INVIVO[@]} ; i++ ))
do
    in="${OUTDIR_INVIVO[$i]}"   
    files="${FILE_LIST_INVIVO[$i]}"
    
    out="${in}/${TYPE_INVIVO[$i]}_wilcox"
    [[ -d "${out}" ]] || mkdir -p "${out}"
    
    # select DRCs / TICs with wilcoxon test
    Rscript ${SCRIPT_DIR}/define_DRC_TIC_wilcox.R ${in}/cpm.tsv ${files} "parental" "tumor" ${out}
    Rscript ${SCRIPT_DIR}/auc_table.R ${in}/cpm.tsv ${out}/list.txt ${out}/auc.txt

    out="${in}/${TYPE_INVIVO[$i]}_fisher"
    [[ -d "${out}" ]] || mkdir -p "${out}"
    
    # select DRCs / TICs with fisher test
    Rscript ${SCRIPT_DIR}/define_DRC_TIC_fisher.R ${in}/cpm.tsv ${files} "parental" "tumor" ${out} 
    Rscript ${SCRIPT_DIR}/auc_table.R ${in}/cpm.tsv ${out}/list.txt ${out}/auc.txt
    
    out="${in}/${TYPE_INVIVO[$i]}_fisher_paired"
    [[ -d "${out}" ]] || mkdir -p "${out}"
    
    # select DRCs / TICs with fisher test
    Rscript ${SCRIPT_DIR}/define_DRC_TIC_fisher_paired.R ${in}/cpm.tsv ${files} "parental" "tumor" ${out} 
    Rscript ${SCRIPT_DIR}/auc_table.R ${in}/cpm.tsv ${out}/list.txt ${out}/auc.txt

    out="${in}/${TYPE_INVIVO[$i]}_gam"
    [[ -d "${out}" ]] || mkdir -p "${out}"
    
    # select DRCs / TICs with gam
    Rscript ${SCRIPT_DIR}/define_DRC_TIC_gam.R ${in}/cpm.tsv ${files} "parental" "tumor" ${out}     
    Rscript ${SCRIPT_DIR}/auc_table.R ${in}/cpm.tsv ${out}/list.txt ${out}/auc.txt

done

for (( i = 0 ; i < ${#OUTDIR_INVITRO[@]} ; i++ ))
do
    in="${OUTDIR_INVITRO[$i]}"   
    files="${FILE_LIST_INVITRO[$i]}"
    
    out="${in}/${TYPE_INVITRO[$i]}_wilcox"
    [[ -d "${out}" ]] || mkdir -p "${out}"
    
    # select DRCs / TICs with wilcoxon test
    Rscript ${SCRIPT_DIR}/define_DRC_TIC_wilcox.R ${in}/cpm.tsv ${files} "early" "late" ${out}
    Rscript ${SCRIPT_DIR}/auc_table.R ${in}/cpm.tsv ${out}/list.txt ${out}/auc.txt

    out="${in}/${TYPE_INVITRO[$i]}_fisher"
    [[ -d "${out}" ]] || mkdir -p "${out}"
    
    # select DRCs / TICs with fisher test
    Rscript ${SCRIPT_DIR}/define_DRC_TIC_fisher.R ${in}/cpm.tsv ${files} "early" "late" ${out} 
    Rscript ${SCRIPT_DIR}/auc_table.R ${in}/cpm.tsv ${out}/list.txt ${out}/auc.txt

    out="${in}/${TYPE_INVITRO[$i]}_gam"
    [[ -d "${out}" ]] || mkdir -p "${out}"
    
    # select DRCs / TICs with gam
    Rscript ${SCRIPT_DIR}/define_DRC_TIC_gam.R ${in}/cpm.tsv ${files} "early" "late" ${out} 
    Rscript ${SCRIPT_DIR}/auc_table.R ${in}/cpm.tsv ${out}/list.txt ${out}/auc.txt

done





