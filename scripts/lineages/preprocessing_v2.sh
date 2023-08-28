
if [[ $# -lt "5" ]]
then
    echo "bash $0 <indir> <outdir> <exp> <sample> <top_clones>"
    exit
fi

INDIR="$1"
OUTDIR="$2"
EXP="$3"
SAMPLE="$4"
TOP_CLONES="$5"

DIR=$( cd $(dirname $0) ; pwd ) 

# input
DF="${INDIR}/${EXP}/clones/cpt.tsv"

# output
out="${OUTDIR}/${EXP}/top_clones"
[[ -d "${out}" ]] || mkdir "${out}"
CLONES="${out}/${SAMPLE}.txt"


Rscript "${DIR}/top_clones.R" "${DF}" "${SAMPLE}" "${TOP_CLONES}" "${CLONES}"

exit

