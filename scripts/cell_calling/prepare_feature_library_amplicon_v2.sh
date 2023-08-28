
# create the list of features in cellranger format

PATTERN="AATTAAGAATT(BC)CCTAGAGGGCCCGTTTAAACCC"
TYPE="CRISPR Guide Capture"
NAME="Non-Targeting"

if [[ $# -lt "2" ]]
then
    echo "Usage: <gbc_list> <out_features>"
fi

GBC=$1
FEATURES=$2

ran=$RANDOM

cut -f 2 ${GBC} > gbc${ran}.txt
awk '{ print ">" $2 "\n" $1 }' "${GBC}" | seqkit seq -p -r -v -t dna | grep -v '^>' > seq${ran}.rev.txt
paste seq${ran}.rev.txt gbc${ran}.txt > gbc_rev${ran}.tsv
rm gbc${ran}.txt seq${ran}.rev.txt

echo "id,name,read,pattern,sequence,feature_type,target_gene_id,target_gene_name" > "${FEATURES}"
awk -v p=${PATTERN} -v t="${TYPE}" -v n=${NAME} '{ print $2 "," $2 ",R2," p "," $1 "," t "," n "," n }' gbc_rev${ran}.tsv >> "${FEATURES}"
rm gbc_rev${ran}.tsv

exit

