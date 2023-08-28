
# update existing GBC list:
# 1. extract the GBCs already used and match them with the sequences collected across samples
# 2. assign new unique GBC identifiers to the sequences that were not already present in the input GBC

IN_GBC_LIST=$1
FILES_GBC=$2
OUTDIRS=$3
OUTDIR=$4

files=( ${FILES_GBC//,/ } )
outdirs=( ${OUTDIRS//,/ } )

if [[ ! ${#files[@]} -eq ${#outdirs[@]} ]] 
then
    echo "Different number of files (${#files[@]}) and dirs (${#outdirs[@]})"
    exit
fi

# LAST_GBC=$(tail -1 ${IN_GBC_LIST} | cut -f 2 | sed "s/GBC//g")
# NEW: sort first!!
LAST_GBC=$(cut -f 2 ${IN_GBC_LIST} | sed "s/GBC//g" | sort -n | tail -1)
cut -f 1 "${OUTDIR}"/*/*/in_with_counts.*unique.txt | sort | uniq > sequences
cat ${IN_GBC_LIST} > seq_and_gbc
join -j1 -v1 sequences ${IN_GBC_LIST} | awk -v LAST_GBC=$LAST_GBC '{ print $1 "	GBC" (NR+LAST_GBC) }' >> seq_and_gbc
sort seq_and_gbc > pippo ; mv pippo seq_and_gbc
join -j1 -t '	' sequences seq_and_gbc > "${OUTDIR}/GBC.txt" 
rm seq_and_gbc sequences

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

