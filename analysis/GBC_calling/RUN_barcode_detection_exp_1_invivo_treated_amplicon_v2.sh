# extract the reads that match the GBC flanking region(s)

PETA_PATH="/hpcnfs/techunits/bioinformatics/software/petagene/petalink_1.3.13/bin/petalink.so"
PETA_REF="/hpcnfs/techunits/bioinformatics/software/petagene/petalink_1.3.13/species"

PREAMBLE="export LD_PRELOAD=${PETA_PATH} ; export PETASUITE_REFPATH=${PETA_REF}"

LOCALDIR=$( cd $( dirname $0 ) ; pwd )
# DIR="${LOCALDIR}/exp_1"
DIR="${LOCALDIR}/exp_1_amplicon_v2"

SCRIPTDIR=$( cd "../../scripts/GBC_calling/" ; pwd )

NCPU="4"
MEM="24G"

FILES_DNA=( "files_raw_exp_1_invivo_A_treated_NoRep.tsv" \
            "files_raw_exp_1_invivo_B_treated_NoRep.tsv" \
            "files_raw_exp_1_invivo_C_treated_NoRep.tsv" )
OUTDIR_DNA=( "${DIR}/1_invivo_A_treated" \
             "${DIR}/1_invivo_B_treated" \
             "${DIR}/1_invivo_C_treated" )

for (( i=0; i<${#FILES_DNA[@]}; i++ ))
do
#    COMMAND="bash ${SCRIPTDIR}/barcode_detection.sh ${FILES_DNA[$i]} ${OUTDIR_DNA[$i]} dna ${NCPU}"
    COMMAND="bash ${SCRIPTDIR}/barcode_detection_amplicon_v2.sh ${FILES_DNA[$i]} ${OUTDIR_DNA[$i]} dna ${NCPU}"
    ID=$( basename ${OUTDIR_DNA[$i]} )
    echo "cd ${LOCALDIR} ; $COMMAND" | qsub -N ${ID} -o ${ID}.STDOUT -e ${ID}.STDERR -l select=1:ncpus=${NCPU}:mem=${MEM} 
done


exit



