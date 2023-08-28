

if [[ $# -lt "4" ]]
then
    echo "bash $0 <dir> <cl_list_1> <cl_list_2> <out>"
    exit
fi

DIR="$1"
CL_LIST_1="$2" # comma-separated list of cluster IDs in sample 1 with high sharedness with clusters in sample 2
CL_LIST_2="$3" # comma-separated list of cluster IDs in sample 2 paired with CL_LIST_1
OUT="$4"

IFS=','
read -a CL1 <<< "${CL_LIST_1}"
read -a CL2 <<< "${CL_LIST_2}"
IFS=$'\n'

CL1_ALL=$( ls "${DIR}" | sed "s/GBC_list_//g" | sed "s/-.*//g" | sort | uniq )
CL2_ALL=$( ls "${DIR}" | sed "s/GBC_list_.*-//g" | sed "s/\.txt//g" | sort | uniq )

if [[ "${#CL1[@]}" -ne "${#CL2[@]}" ]]
then 
    echo "The number of clusters is not the same"
    exit
fi

in_cl="in_cl.txt"
out_cl="out_cl.txt"

PAIRS=()
for (( i=0; i<${#CL1}; i++ ))
do
    PAIRS=( ${PAIRS[@]} "${CL1[$i]}-${CL2[$i]}" )
done

# in pairs
for (( i=0; i<${#PAIRS[@]}; i++ ))
do
    cat "${DIR}/GBC_list_${PAIRS[$i]}.txt" 
done | sort | uniq > "${in_cl}" 

# NEW: scan all clusters, not only the transcriptionally identifiable ones!!!
# out of pairs
for (( i=0; i<${#CL1_ALL[@]}; i++ ))
do
    for (( j=0; j<${#CL2_ALL[@]}; j++ ))
    do 
        pair="${CL1_ALL[$i]}-${CL2_ALL[$j]}"
        [[ " ${PAIRS[@]} " =~ " ${pair} " ]] || cat "${DIR}/GBC_list_${pair}.txt"
    done 
done | sort | uniq > "${out_cl}" 

# keep only the GBCs that are found only in the pairs
comm -2 -3 "${in_cl}" "${out_cl}" > "${OUT}" 
 
rm "${in_cl}" "${out_cl}"

exit

