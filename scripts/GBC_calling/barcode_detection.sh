
set -e -o pipefail
export LC_ALL=C

# collect the reads via pattern matching

SEQ_DNA="GGGTTTAAACGGGCCCTCTAGG"
PATTERN_DNA="GGGTTTAAACGGGCCCTCTAGGNNNNNNNNNNNNNNNNNNAATTCTTAATT" # NEW: added an 'A' at the beginning of the right flanking region

SEQ_RNA="TAGCAAACTGGGGCACAAGCTTAATTAAGAATT"
PATTERN_RNA="TAGCAAACTGGGGCACAAGCTTAATTAAGAATTNNNNNNNNNNNNNNNNNN" # NEW: shorten the stretch of 'N's to adapt to different library sizes

if [[ $# -lt 3 ]]
then
    echo "Usage: <files.tsv> <outdir> <dna|rna> [<ncpu>]"
    exit
fi

FILES=$1
BASEDIR=$2
TYPE=$3
if [[ $# -gt 3 ]] ; then NCPU="4" ; else NCPU=$4 ; fi

if [ "${TYPE}" = "dna" ]
then
    SEQ="${SEQ_DNA}"
    PATTERN="${PATTERN_DNA}"
else 
    if [ "${TYPE}" = "rna" ]
    then
        SEQ="${SEQ_RNA}"
        PATTERN="${PATTERN_RNA}"
    else
        echo "Unknown type"
        exit
    fi
fi

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
echo "seq: ${SEQ}"
echo "pattern: ${PATTERN}"

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
        echo "seqkit on ${dir}/${file}"
        seqkit grep -j ${NCPU} -s -r -i -p "${SEQ}" "${dir}/${file}" | seqkit locate -j ${NCPU} -i -d -p "${PATTERN}" >> "${TABLE}"
    done

    ############################# COUNT THE FREQUENCY OF EACH READ

    echo "extract counts"
    cut -f 7 "${TABLE}" | grep -v "matched" | sort | uniq -c > "${COUNTS}"
    echo "done!"

    rm -f "${TABLE}"

    i=$((i+1))

done

rmdir "${TABLEDIR}"
rmdir "${SCRATCH}"

exit


