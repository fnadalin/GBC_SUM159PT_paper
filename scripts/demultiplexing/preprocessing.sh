
FASTQSPLIT="$HOME/fnadalin_scratch/exp_FINAL/scripts/cpp/fastqsplit/fastqsplit"
NSEQ="20000000"

if [ $# -lt 1 ]
then
    echo "Usage: bash preprocessing.sh <sample>"
    exit
fi

SAMPLE=$1

OUTDIR="/scratch/$(whoami)/$SAMPLE"
[[ -d "${OUTDIR}" ]] || mkdir -p "${OUTDIR}"

cd "${SAMPLE}"

LANES="lanes_orig.csv"
LANES_CHUNK="lanes.csv"
head -1 ${LANES} > ${LANES_CHUNK}

i=0
IFS=$','
while read -r lane r1 r2
do
    if [ "$i" -eq "0" ] ; then i=$((i+1)) ; continue ; fi
   
    echo "$lane $r1 $r2"

    R1="${OUTDIR}/R1.fastq"
    R2="${OUTDIR}/R2.fastq"
    gunzip -c "${r1}" > "${R1}"
    gunzip -c "${r2}" > "${R2}"

    out="${OUTDIR}/${i}"
    [[ -d "${out}" ]] || mkdir "${out}"
    $FASTQSPLIT ${R1} ${R2} $NSEQ "${out}"

    nchunks=$(ls "${out}" | grep -c "R1")
    for (( j=0; j<$nchunks; j++ ))  
    do
        echo "L$i,$out/R1_chunk$j.fastq,$out/R2_chunk$j.fastq"
    done >> ${LANES_CHUNK}

    i=$((i+1))
done < "${LANES}"

rm "${OUTDIR}/R1.fastq" "${OUTDIR}/R2.fastq"

exit


