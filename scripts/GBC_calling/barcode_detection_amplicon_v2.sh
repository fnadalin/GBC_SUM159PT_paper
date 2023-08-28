
set -e -o pipefail
export LC_ALL=C

# collect the reads via pattern matching

PATTERN_FW="GGGTTTAAACGGGCCCTCTAGG"
PATTERN_RV="AATTAAGAATT"
S_DNA="23" # start
E_DNA="40" # end
S_RNA="12"
E_RNA="29"

if [[ $# -lt 3 ]]
then
    echo "Usage: <files.tsv> <outdir> <dna|rna> [<threads>]"
    exit
fi

FILES=$1
BASEDIR=$2
TYPE=$3
if [[ $# -lt 4 ]] ; then THREADS="4" ; else THREADS=$4 ; fi

DATA=( $(cut -f 2 "${FILES}") )
NAMES=( $(cut -f 1 "${FILES}") )

COUNTDIR="${BASEDIR}/counts"
[[ -d "${COUNTDIR}" ]] || mkdir -p "${COUNTDIR}"

SCRATCH="/scratch/ieo5369_$RANDOM"
TABLEDIR="${SCRATCH}/tables"
[[ -d "${TABLEDIR}" ]] || mkdir -p "${TABLEDIR}"

echo "########### barcode_calling ############"
echo "dir: ${TABLEDIR}"
echo "uname: $(uname -a)"

i=0
for dir in ${DATA[@]} 
do
    echo "########### PROCESSING DIR $dir ############"

    FILENAME=${NAMES[$i]}
    TABLE="${TABLEDIR}/${FILENAME}_table.txt"
    COUNTS="${COUNTDIR}/${FILENAME}_counts.txt"

    echo "sample ${NAMES[$i]}"
    
    ############################# FIND THE PATTERN AND CONSTRUCT THE TABLE

    touch "${TABLE}"
    for file in $(ls ${dir} | grep "R2_001.fastq.gz$")
    do
        echo "seqkit on ${dir}/${file}" # NEW: run amplicon
        if [ "${TYPE}" == "dna" ]
        then
            seqkit-v2.1.0 amplicon -j ${THREADS} -s -m 1 -F ${PATTERN_FW} -R ${PATTERN_RV} -r ${S_DNA}:${E_DNA} --bed "${dir}/${file}" >> "${TABLE}"
        else
            seqkit-v2.1.0 amplicon -j ${THREADS} -s -m 1 -F ${PATTERN_RV} -R ${PATTERN_FW} -P -r ${S_RNA}:${E_RNA} --bed "${dir}/${file}" >> "${TABLE}"
        fi
    done

    ############################# COUNT THE FREQUENCY OF EACH READ

    echo "extract counts"
    cut -f 7 "${TABLE}" | sort | uniq -c > "${COUNTS}"
    echo "done!"

    rm -f "${TABLE}"

    i=$((i+1))

done

rmdir "${TABLEDIR}"
rmdir "${SCRATCH}"

exit


