
# collect the hub IDs and sequences for cell ranger

FILEDIR=$1
OUTDIR=$2
DIST=$3
FRAC=$4
INFILE=$5
OUTFILE=$6

FRAC=${FRAC/./}

IFS=$'\n'
for line in $(cat "${FILEDIR}")
do
    IFS='	'
    read -a v <<< "${line}"
    SAMPLE="${v[0]}"

    outdir="${OUTDIR}/${SAMPLE}/out_d${DIST}_f${FRAC}"
    groups="${outdir}/groups.tsv"
    grep -v hub "${groups}" | cut -f 3
   
    IFS=$'\n'

done | sort | uniq > hubs

sort -k2 "${INFILE}" > all_gbc
join -1 2 -2 1 -o 1.1 1.2 -t '	' all_gbc hubs > "${OUTFILE}"

rm all_gbc hubs

exit

