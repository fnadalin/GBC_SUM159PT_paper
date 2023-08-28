
if [ $# -lt 1 ]
then
    echo "Usage: bash preprocessing.sh <sample>"
    exit
fi

SAMPLE=$1

OUTDIR="/scratch/$(whoami)/$SAMPLE"

cd $SAMPLE
LANES="lanes_orig.csv"

i=0
IFS=$','
while read -r lane r1 r2
do
    if [ "$i" -eq "0" ] ; then i=$((i+1)) ; continue ; fi

    out="${OUTDIR}/${i}"
    rm ${out}/R*
    rmdir ${out}
    i=$((i+1))

done < "${LANES}"

exit


