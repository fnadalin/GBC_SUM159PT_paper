
# collect the sequences across samples and assign unique GBC identifiers

FILES_GBC=$1
OUTDIRS=$2
OUTDIR=$3

files=( ${FILES_GBC//,/ } )
outdirs=( ${OUTDIRS//,/ } )

if [[ ! ${#files[@]} -eq ${#outdirs[@]} ]] 
then
    echo "Different number of files (${#files[@]}) and dirs (${#outdirs[@]})"
    exit
fi

cut -f 1 "${OUTDIR}"/*/*/in_with_counts.*unique.txt | sort | uniq | awk '{ print $1 "	GBC" NR }' > "${OUTDIR}/GBC.txt"

for (( i=0; i < ${#files[@]}; i++ ))
do
    IFS=$'\n'
    for line in $(cat "${files[$i]}")
    do
        IFS='	'
        read -a v <<< "${line}"
        SAMPLE=${v[0]}
        FILE=${v[1]} 

        outdir="${OUTDIR}/${outdirs[$i]}/${SAMPLE}"
        join -j1 -t '	' -o 1.2,2.1,2.2 "${OUTDIR}/GBC.txt" "${outdir}/in_with_counts.unique.txt" > "${outdir}/in_with_counts.GBC.txt"

        IFS=$'\n'
    done
done

exit

