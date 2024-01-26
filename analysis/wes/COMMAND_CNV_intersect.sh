CNVKIT_DIR="../../data/SUM159_PACLI_WES_ANALYSIS/CNVKIT"
BEDTOOLS="$HOME/francesca_sw/bedtools/bedtools"
OUT="in"

[[ -d ${OUT} ]] || mkdir ${OUT}

cut -f 1-3,5 ${CNVKIT_DIR}/D15A_markedDuplicates_19168.call.cns.txt | grep -v chromosome > ${OUT}/fileA.bed
cut -f 1-3,5 ${CNVKIT_DIR}/D15B_markedDuplicates_19169.call.cns.txt | grep -v chromosome > ${OUT}/fileB.bed 
cut -f 1-3,5 ${CNVKIT_DIR}/D15C_markedDuplicates_19170.call.cns.txt | grep -v chromosome > ${OUT}/fileC.bed 

$BEDTOOLS intersect -a ${OUT}/fileA.bed -b ${OUT}/fileB.bed > X
$BEDTOOLS intersect -a X -b ${OUT}/fileC.bed > Y
$BEDTOOLS intersect -F 0.8 -wa -a Y -b ${OUT}/fileA.bed > A
$BEDTOOLS intersect -F 0.8 -wa -a A -b ${OUT}/fileB.bed > AB
$BEDTOOLS intersect -F 0.8 -wa -a AB -b ${OUT}/fileC.bed > ABC
cut -f 1-3 ABC > ${OUT}/CNV_intersection_80perc.bed
rm X Y A AB ABC

$BEDTOOLS intersect -wb -a ${OUT}/CNV_intersection_80perc.bed -b ${OUT}/fileA.bed | cut -f 1-3,7 > X
$BEDTOOLS intersect -wb -a X -b ${OUT}/fileB.bed | cut -f 1-4,8 > Y
$BEDTOOLS intersect -wb -a Y -b ${OUT}/fileC.bed | cut -f 1-5,9 | grep -v chrY > ${OUT}/CNV_intersection_80perc.bed
rm X Y

exit

