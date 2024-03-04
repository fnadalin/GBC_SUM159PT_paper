SCRIPT_DIR=$( cd "../../scripts/3CA" ; pwd )
DIR=$( cd "../../analysis/gene_expression" ; pwd )
DIR_MDA=$( cd "../../analysis/gene_expression_MDA" ; pwd )

TABLE_3CA="../../data/3CA/meta_programs_2023-07-13.tsv"

OUT_DIR=( "19924_k40_res0.5" "19924_k40_res0.1" )

EXP="19924"
OBJ="${DIR_MDA}/${EXP}/${EXP}_ccRegress_filt/hvg_pca_clust/object.Rds"

MODE=( "clusters_pca_vst_top1000_k40_res0.5" "clusters_pca_vst_top1000_k40_res0.1" )

SIGN_DIR="${DIR}/signatures/DEG_P0_ccRegress"

SUBPOP_ID=( "0,1,2,3,4,5,6" "0,1,2" )
SUBPOP_NAME=( "0,1,2,3,4,5,6" "0,1,2" )

MEM="20000"

for (( i=0; i<${#OUT_DIR[@]}; i++ ))
do
    ID="MP_${OUT_DIR[$i]}"
    COMMAND="${SCRIPT_DIR}/meta_programs_vs_signatures.sh"
    ARGS=( ${TABLE_3CA} ${OBJ} ${MODE[$i]} ${SUBPOP_ID[$i]} ${SUBPOP_NAME[$i]} ${OUT_DIR[$i]} ${SIGN_DIR} )
    bsub -J ${ID} -M ${MEM} -o ${ID}.STDOUT -e ${ID}.STDERR bash ${COMMAND} ${ARGS[@]}
done

exit

