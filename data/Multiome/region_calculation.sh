
DIR="../../analysis/multiome/cisTopic_all/IDR_files"
PREFIX="1C_ARC_seurat_cisTopicObject_1E_ARC_seurat_cisTopicObject"

SELECTED_TOPICS_1C=( 2 3 4 5 6 8 9 10 11 12 13 14 15 17 18 20 21 22 23 24 25 26 27 28 29 31 32 34 35 36 37 )
SELECTED_TOPICS_1E=( 1 2 3 4 5 6 8 7 9 10 11 12 13 14 15 16 18 19 22 23 24 26 27 28 29 31 32 34 35 37 38 39 40 )

TOPIC_1C=( 2 31 27 15 5 )
TOPIC_1E=( 16 15 14 19 1 )
MODULE=( 1 2 8 20 4 )

# N.B.: common peaks are the UNION of overlapping peaks between replicates

for i in ${SELECTED_TOPICS_1C[@]}
do
    for j in ${SELECTED_TOPICS_1E[@]}
    do
        f="${DIR}/${PREFIX}_${i}_${j}_IDR_0.05.bed"
        if [[ -e $f ]]
        then
            cut -f 1-3 $f
        fi
    done
done | sort | uniq > "cistopic_repr_regions_non_redundant.bed"

for i in 0 1 2 3 4
do
    f="${DIR}/${PREFIX}_${TOPIC_1C[$i]}_${TOPIC_1E[i]}_IDR_0.05.bed"
    cut -f 1-5,11,14 $f > MODULE_${MODULE[$i]}.bed
done

