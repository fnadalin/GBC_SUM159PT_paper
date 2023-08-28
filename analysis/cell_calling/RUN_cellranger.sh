# run cellranger 

SING_PATH="/hpcnfs/software/singularity/3.5.3/bin/"
BIOINFO_SING="/hpcnfs/techunits/bioinformatics/singularity/genomicsunit2_latest.sif" # cellranger 6.0

TSCRT="/hpcnfs/techunits/bioinformatics/refdata/cellranger-gex/GRCh38/2020-A/"
ARCREF="/hpcnfs/techunits/bioinformatics/refdata/cellranger-arc/GRCh38/2020-A-2.0.0/"

LOCALDIR=$( cd $( dirname $0 ) ; pwd )
GBCDIR=$( cd "../../analysis/GBC_calling/exp_1_amplicon_v2" ; pwd )
SCRIPTDIR=$( cd "../../scripts/cell_calling/" ; pwd )

NCPU="4"
MEM1="5G"
MEM2="64G"

DIST="1"
FRAC="1"

GBC_LIST=( "${GBCDIR}/GBC_1B_RNA_d${DIST}_f${FRAC}.tsv" \
           "${GBCDIR}/GBC_1D_RNA_d${DIST}_f${FRAC}.tsv" )
LIB_GEX=( "libraries_exp_1B.tsv" \
          "libraries_exp_1D.tsv" )
OUT_GEX=( "1B_GEX" \
          "1D_GEX" )

GBC_LIST_NUC=( "${GBCDIR}/GBC_1C_RNA_d${DIST}_f${FRAC}.tsv" \
               "${GBCDIR}/GBC_1E_RNA_d${DIST}_f${FRAC}.tsv" )
LIB_GEX_NUC=( "libraries_exp_1C.tsv" \
              "libraries_exp_1E.tsv" )
OUT_GEX_NUC=( "1C_GEX" \
              "1E_GEX" )

LIB_ARC=( "libraries_exp_1C_ARC.tsv" \
          "libraries_exp_1E_ARC.tsv" )
OUT_ARC=( "1C_ARC" \
          "1E_ARC" )


IFS=$'\n'
for (( i=0; i<${#LIB_GEX[@]}; i++ ))
do
    [[ -d "${OUT_GEX[$i]}" ]] || mkdir -p "${OUT_GEX[$i]}"

    #### PREPARE FEATURE LIBRARIES

    FEAT="${LOCALDIR}/${OUT_GEX[$i]}/features.csv"
    COMMAND="bash ${SCRIPTDIR}/prepare_feature_library_amplicon_v2.sh ${GBC_LIST[$i]} ${FEAT}"
    ID="prepare_feat_${OUT_GEX[$i]}"
    JOBID=$(echo "cd ${LOCALDIR} ; $COMMAND" | qsub -N ${ID} -o ${ID}.STDOUT -e ${ID}.STDERR -l select=1:ncpus=1:mem=${MEM1}) 

    # RUN CELLRANGER COUNT

    for line in $(cat ${LIB_GEX[$i]})
    do
        RAND=${RANDOM}
        TMPDIR="/scratch/ieo5369/${RAND}"

        SAMPLE_ID=$(echo $line | cut -f 1)
        LIBRARY=${LOCALDIR}/$(echo $line | cut -f 2)
        CELLS=$(echo $line | cut -f 3)        

        MEM=${MEM2/G/}
        COMMAND="cellranger count --nosecondary --localcores ${NCPU} --localmem ${MEM} --expect-cells ${CELLS} \
                 --id=${SAMPLE_ID} --transcriptome=${TSCRT} --libraries=${LIBRARY} --feature-ref=${FEAT} --no-bam"
        COMMAND_MV="mv ${TMPDIR}/* ${LOCALDIR}/${OUT_GEX[$i]}/ ; rmdir ${TMPDIR}"

        ID="cr_count_${OUT_GEX[$i]}_${SAMPLE_ID}"
        echo "export PATH=${SING_PATH}:\$PATH 
              [[ -d "${TMPDIR}" ]] || mkdir -p "${TMPDIR}"
              cd ${TMPDIR} 
              singularity exec --bind /hpcnfs/ ${BIOINFO_SING} ${COMMAND}
              ${COMMAND_MV}" | qsub -o ${ID}.STDOUT -e ${ID}.STDERR -l select=1:ncpus=${NCPU}:mem=${MEM2} -N ${ID} -W depend=afterok:${JOBID}

    done

done



IFS=$'\n'
for (( i=0; i<${#LIB_GEX_NUC[@]}; i++ ))
do
    [[ -d "${OUT_GEX_NUC[$i]}" ]] || mkdir -p "${OUT_GEX_NUC[$i]}"

    #### PREPARE FEATURE LIBRARIES

    FEAT="${LOCALDIR}/${OUT_GEX_NUC[$i]}/features.csv"
    COMMAND="bash ${SCRIPTDIR}/prepare_feature_library_amplicon_v2.sh ${GBC_LIST_NUC[$i]} ${FEAT}"
    ID="prepare_feat_${OUT_GEX_NUC[$i]}"
    JOBID=$(echo "cd ${LOCALDIR} ; $COMMAND" | qsub -N ${ID} -o ${ID}.STDOUT -e ${ID}.STDERR -l select=1:ncpus=1:mem=${MEM1})

    # RUN CELLRANGER COUNT ON NUCLEI

    for line in $(cat ${LIB_GEX_NUC[$i]})
    do
        RAND=${RANDOM}
        TMPDIR="/scratch/ieo5369/${RAND}"
        
        SAMPLE_ID=$(echo $line | cut -f 1)
        LIBRARY=${LOCALDIR}/$(echo $line | cut -f 2)
        CELLS=$(echo $line | cut -f 3)

        MEM=${MEM2/G/}
        COMMAND="cellranger count --nosecondary --localcores ${NCPU} --localmem ${MEM} --include-introns --chemistry=ARC-v1 --expect-cells ${CELLS} \
                 --id=${SAMPLE_ID} --transcriptome=${TSCRT} --libraries=${LIBRARY} --feature-ref=${FEAT} --no-bam"
        COMMAND_MV="mv ${TMPDIR}/* ${LOCALDIR}/${OUT_GEX_NUC[$i]}/ ; rmdir ${TMPDIR}"

        ID="cr_count_${OUT_GEX_NUC[$i]}_${SAMPLE_ID}"
        echo "export PATH=${SING_PATH}:\$PATH 
              [[ -d "${TMPDIR}" ]] || mkdir -p "${TMPDIR}"
              cd ${TMPDIR} 
              singularity exec --bind /hpcnfs/ ${BIOINFO_SING} ${COMMAND}
              ${COMMAND_MV}" | qsub -o ${ID}.STDOUT -e ${ID}.STDERR -l select=1:ncpus=${NCPU}:mem=${MEM2} -N ${ID} -W depend=afterok:${JOBID}

    done

done



IFS=$'\n'
for (( i=0; i<${#LIB_ARC[@]}; i++ ))
do
    [[ -d "${OUT_ARC[$i]}" ]] || mkdir -p "${OUT_ARC[$i]}"

    # RUN CELLRANGER ARC COUNT

    for line in $(cat ${LIB_ARC[$i]})
    do
        RAND=${RANDOM}
        TMPDIR="/scratch/ieo5369/${RAND}"

        SAMPLE_ID=$(echo $line | cut -f 1)
        LIBRARY=${LOCALDIR}/$(echo $line | cut -f 2)
        CELLS=$(echo $line | cut -f 3)

        MEM=${MEM2/G/}
        COMMAND="cellranger-arc count --localcores ${NCPU} --localmem ${MEM} --id=${SAMPLE_ID} --reference=${ARCREF} --libraries=${LIBRARY}"
        COMMAND_MV="mv ${TMPDIR}/* ${LOCALDIR}/${OUT_ARC[$i]}/ ; rm ${TMPDIR}"

        ID="cr_count_${OUT_ARC[$i]}_${SAMPLE_ID}"
        echo "export PATH=${SING_PATH}:\$PATH 
              [[ -d "${TMPDIR}" ]] || mkdir -p "${TMPDIR}"
              cd ${TMPDIR} 
              singularity exec --bind /hpcnfs/ ${BIOINFO_SING} ${COMMAND}
              ${COMMAND_MV}" | qsub -o ${ID}.STDOUT -e ${ID}.STDERR -l select=1:ncpus=${NCPU}:mem=${MEM2} -N ${ID}

    done
exit
done



exit



