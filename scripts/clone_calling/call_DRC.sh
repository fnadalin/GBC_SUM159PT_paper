
if [[ $# -lt "2" ]]
then
    echo "bash $0 <file_list> <out>"
    exit
fi

FILE_LIST="$1"
OUT="$2"

DIR=$( cd $(dirname $0) ; pwd ) 

Rscript "${DIR}/call_DRC.R" "${FILE_LIST}" "${OUT}" 

exit

