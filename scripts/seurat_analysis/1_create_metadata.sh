

DIR=$( cd $(dirname $0) ; pwd )

if [[ $# -lt "5" ]]
then
    echo "bash $0 <script_dir> <matrix_list> <in_dir> <out_dir> <clone_info>"
    exit
fi

SCRIPT_DIR="$1"
MATRIX_LIST="$2"
IN_DIR="$3"
OUT_DIR="$4"
CLONE_INFO="$5" # table with clone_info_name (TIC, DRC_sc, DRC_bulk), the file containing the info, and the type (sc or bulk)

METADATA="${OUT_DIR}/CB_metadata.tsv"

# create metadata file
echo "create metadata file"
echo "sample.name	UMIcount	readcount	expr.GBC.num	expr.GBC.list" > ${METADATA}
IFS=$'\n'
for line in $(cat ${MATRIX_LIST})
do
    S=$(echo $line | cut -f 1)
    CB_CLASS_TSV="${IN_DIR}/$S/classification/CB_classification.tsv"

    grep -v "UMIcount" ${CB_CLASS_TSV} | \
    awk -F "	" -v s=$S -v t="TRUE" -v u="uninfected" -v d="doublet" \
        '{ if ($4 == t && $8 != u && $8 != d) print s "-" $1 "	" s "	" $2 "	" $3 "	" $5 "	" $6}'  
    
done >> ${METADATA}

IFS=$'\n'
for line in $(cat ${CLONE_INFO})
do
    NAME=$(echo $line | cut -f 1)
    FILE=$(echo $line | cut -f 2)
    TYPE=$(echo $line | cut -f 3)
    
    if [[ ${TYPE} == "sc" ]]
    then
        Rscript ${DIR}/add_clone_label_sc.R ${FILE} ${METADATA} ${NAME}
    elif [[ ${TYPE} == "bulk" ]]
    then
#        Rscript ${DIR}/add_clone_label_bulk.R ${FILE} ${METADATA} ${NAME}
        Rscript ${DIR}/add_clone_label_bulk_v2.R ${FILE} ${METADATA} ${NAME}
    fi

done

# generate cell lists
echo "generate cell lists"
IFS=$'\n'
for S in $(cut -f 1 ${MATRIX_LIST})
do
    CELL_LIST="${OUT_DIR}/cell_list_${S}.txt"
    
    grep -v "UMIcount" ${METADATA} | awk -v s=$S '{ if ($2 == s) print $1 }' | sed "s/^$S-//g" > ${CELL_LIST}
done 


exit


