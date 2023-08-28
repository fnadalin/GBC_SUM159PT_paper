
# link read count information with the unique GBC identifier

SCRIPT_DIR=$1
FILEDIR=$2
OUTDIR=$3
DIST=$4
FRAC=$5

BARCODE_GROUPS="${SCRIPT_DIR}/barcode-groups"
BARCODE_GROUPS_PLOT="${SCRIPT_DIR}/barcode-groups-plots.R"

IFS=$'\n'
for line in $(cat "${FILEDIR}")
do
    IFS='	'
    read -a v <<< "${line}"
    SAMPLE=${v[0]}
    FILE=${v[1]} 

    echo ${SAMPLE}

    indir="${OUTDIR}/${SAMPLE}"
    outdir="${OUTDIR}/${SAMPLE}/out_d${DIST}_f${FRAC/./}"
    [[ -d "${outdir}" ]] || mkdir "${outdir}"

    "${BARCODE_GROUPS}" "${indir}/in_with_counts.GBC.txt" ${DIST} ${FRAC} "${outdir}" > "${outdir}"/STDOUT 2> "${outdir}"/STDERR

    MIN_GROUP_COUNT=$(grep "Min group count" "${outdir}"/STDOUT | sed "s/Min group count: //g")
    Rscript "${BARCODE_GROUPS_PLOT}" "${outdir}" "${outdir}"/plots "${SAMPLE}" ${MIN_GROUP_COUNT} > "${outdir}"/R.STDOUT 2> "${outdir}"/R.STDERR

    IFS=$'\n'

done


exit

