

if [[ $# -lt "5" ]]
then
    echo "bash $0 <file_list> <in_dir> <out_dir> <class1> <class2>"
    exit
fi

FILE_LIST="$1"
IN_DIR="$2"
OUT_DIR="$3"
CLASS1="$4"
CLASS2="$5"

DIR=$( cd $(dirname $0) ; pwd ) 

cpm="${IN_DIR}/cpm.tsv"
for s1 in $(awk -v c=${CLASS1} '{ if ($3 == c) print $1 }' ${FILE_LIST})
do 
    for s2 in $(awk -v c=${CLASS2} '{ if ($3 == c) print $1 }' ${FILE_LIST})
    do 
        out="${OUT_DIR}/cum_dist_${s1}-${s2}.pdf"
        Rscript ${DIR}/plot_cumulative_cpm.R "${cpm}" ${s2} ${s1} "${out}"
    done
done

exit


