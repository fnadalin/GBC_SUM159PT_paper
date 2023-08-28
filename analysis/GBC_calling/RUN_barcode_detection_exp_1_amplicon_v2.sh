# extract the reads that match the GBC flanking region(s)

LOCALDIR=$( cd $( dirname $0 ) ; pwd )
# DIR="${LOCALDIR}/exp_1"
DIR="${LOCALDIR}/exp_1_amplicon_v2"

SCRIPTDIR=$( cd "../../scripts/GBC_calling/" ; pwd )

NCPU="4"
MEM="24G"

FILES_DNA=( "files_raw_exp_1B_DNA.tsv" \
            "files_raw_exp_1C_DNA.tsv" \
            "files_raw_exp_1_invivo_A.tsv" \
            "files_raw_exp_1_invivo_B.tsv" \
            "files_raw_exp_1_invivo_C.tsv" \
            "files_raw_exp_1_invivo_D.tsv" )
OUTDIR_DNA=( "${DIR}/1B_DNA" \
             "${DIR}/1C_DNA" \
             "${DIR}/1_invivo_A" \
             "${DIR}/1_invivo_B" \
             "${DIR}/1_invivo_C" \
             "${DIR}/1_invivo_D" )

FILES_RNA=( "files_raw_exp_1B_RNA.tsv" \
            "files_raw_exp_1C_RNA.tsv" \
            "files_raw_exp_1D_RNA.tsv" \
            "files_raw_exp_1E_RNA.tsv" )
OUTDIR_RNA=( "${DIR}/1B_RNA" \
             "${DIR}/1C_RNA" \
             "${DIR}/1D_RNA" \
             "${DIR}/1E_RNA" )


for (( i=0; i<${#FILES_DNA[@]}; i++ ))
do
#    COMMAND="bash ${SCRIPTDIR}/barcode_detection.sh ${FILES_DNA[$i]} ${OUTDIR_DNA[$i]} dna ${NCPU}"
    COMMAND="bash ${SCRIPTDIR}/barcode_detection_amplicon_v2.sh ${FILES_DNA[$i]} ${OUTDIR_DNA[$i]} dna ${NCPU}"
    ID=$( basename ${OUTDIR_DNA[$i]} )
    echo "cd ${LOCALDIR} ; $COMMAND" | qsub -N ${ID} -o ${ID}.STDOUT -e ${ID}.STDERR -l select=1:ncpus=${NCPU}:mem=${MEM} 
done

for (( i=0; i<${#FILES_RNA[@]}; i++ ))
do
#    COMMAND="bash ${SCRIPTDIR}/barcode_detection.sh ${FILES_RNA[$i]} ${OUTDIR_RNA[$i]} rna ${NCPU}"
    COMMAND="bash ${SCRIPTDIR}/barcode_detection_amplicon_v2.sh ${FILES_RNA[$i]} ${OUTDIR_RNA[$i]} rna ${NCPU}"
    ID=$( basename ${OUTDIR_RNA[$i]} )
    echo "cd ${LOCALDIR} ; $COMMAND" | qsub -N ${ID} -o ${ID}.STDOUT -e ${ID}.STDERR -l select=1:ncpus=${NCPU}:mem=${MEM}
done


exit



