# extract the reads that match the GBC flanking region(s)

LOCALDIR=$( cd $( dirname $0 ) ; pwd )
# DIR="${LOCALDIR}/exp_4"
DIR="${LOCALDIR}/exp_5_amplicon_v2"

SCRIPTDIR=$( cd "../../scripts/GBC_calling/" ; pwd )

NCPU="4"
MEM="24G"

FILES_DNA=( "files_raw_exp_5A_DNA.tsv" )
OUTDIR_DNA=( "${DIR}/5A_DNA" )

for (( i=0; i<${#FILES_DNA[@]}; i++ ))
do
#    COMMAND="bash ${SCRIPTDIR}/barcode_detection.sh ${FILES_DNA[$i]} ${OUTDIR_DNA[$i]} dna ${NCPU}"
    COMMAND="bash ${SCRIPTDIR}/barcode_detection_amplicon_v2.sh ${FILES_DNA[$i]} ${OUTDIR_DNA[$i]} dna ${NCPU}"
    ID=$( basename ${OUTDIR_DNA[$i]} )
    echo "cd ${LOCALDIR} ; $COMMAND" | qsub -N ${ID} -o ${ID}.STDOUT -e ${ID}.STDERR -l select=1:ncpus=${NCPU}:mem=${MEM} 
done


exit



