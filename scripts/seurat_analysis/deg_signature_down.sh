

if [[ $# -lt "3" ]]
then
    echo ""
    echo "Usage: bash $0 <dea1> <dea2> <out> <log2fc> <top>"
    echo ""
    exit
fi

SCRIPT_DIR=$( cd $(dirname $0) ; pwd )

Rscript "${SCRIPT_DIR}/deg_signature_down.R" $@

exit


