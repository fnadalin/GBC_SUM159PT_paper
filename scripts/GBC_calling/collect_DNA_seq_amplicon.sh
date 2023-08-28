
# collect unique DNA barcode sequences 

FILEDIR=$1
OUTDIR=$2

IFS=$'\n'
for line in $(cat "${FILEDIR}")
do
    IFS='	'
    read -a v <<< "${line}"
    SAMPLE=${v[0]}
    FILE=${v[1]} 

    echo "${SAMPLE}"

    outdir="${OUTDIR}/${SAMPLE}"
    [[ -d "${outdir}" ]] || mkdir "${outdir}"

    sed "s/^ \+//g" "${FILE}" | sed "s/ \+/	/g" | cut -f 1 > "${outdir}/in_count_only.txt"
    sed "s/^ \+//g" "${FILE}" | sed "s/ \+/	/g" | cut -f 2 > "${outdir}/in_sequence_only.txt" 

    paste "${outdir}/in_sequence_only.txt" "${outdir}/in_count_only.txt" | sort > "${outdir}/in_with_counts.txt"

    IFS=$'\n'
    s=""
    c=0
    for line in $(cat "${outdir}/in_with_counts.txt")
    do
       IFS='	'
       read -a v <<< "${line}"
       s2=${v[0]}
       c2=${v[1]}
       if [ "$s" != "${s2}" ]
       then
           [ "$s" == "" ] || echo "$s	$c"
           s="${s2}"
           c=0
       fi
       let "c = c + c2"
       IFS=$'\n'
    done > "${outdir}/in_with_counts.unique.txt"
    echo "$s	$c" >> "${outdir}/in_with_counts.unique.txt"

done

exit

