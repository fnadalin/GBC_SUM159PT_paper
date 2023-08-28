

if [[ $# -lt "4" ]]
then
    echo ""
    echo "Usage: bash $0 <log2fc> <top> <out> <deaN>"
    echo ""
    exit
fi

SCRIPT_DIR=$( cd $(dirname $0) ; pwd )

Rscript "${SCRIPT_DIR}/deg_signature_up_topSignificant.R" $@

exit


