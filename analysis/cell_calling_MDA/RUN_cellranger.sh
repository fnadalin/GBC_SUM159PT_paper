# run cellranger 

SING_PATH="/hpcnfs/software/singularity/3.5.3/bin/"
BIOINFO_SING="/hpcnfs/techunits/bioinformatics/singularity/genomicsunit2_latest.sif" # cellranger 6.0

TSCRT="/hpcnfs/techunits/bioinformatics/refdata/cellranger-gex/GRCh38/2020-A/"

LOCALDIR=$( cd $( dirname $0 ) ; pwd )
SCRIPTDIR=$( cd "../../scripts/cell_calling/" ; pwd )

NCPU="8"
MEM1="5G"
MEM2="64G"

LIBRARY="library_19924.csv"
OUT_GEX="MDA"
SAMPLE_ID="19924"
CELLS="3000"

[[ -d "${OUT_GEX}" ]] || mkdir -p "${OUT_GEX}"

# RUN CELLRANGER COUNT

RAND=${RANDOM}
TMPDIR="/scratch/ieo5369/${RAND}"

MEM=${MEM2/G/}
COMMAND="cellranger count --nosecondary --localcores ${NCPU} --localmem ${MEM} --expect-cells ${CELLS} \
         --id ${SAMPLE_ID} --transcriptome=${TSCRT} --libraries=${LIBRARY} --no-bam"
COMMAND_MV="mv ${TMPDIR}/* ${LOCALDIR}/${OUT_GEX}/ ; rmdir ${TMPDIR}"

ID="cr_count_${OUT_GEX}_${SAMPLE_ID}"
echo "export PATH=${SING_PATH}:\$PATH 
[[ -d \"${TMPDIR}\" ]] || mkdir -p \"${TMPDIR}\"
cp ${LOCALDIR}/${LIBRARY} ${TMPDIR}/
cd ${TMPDIR} 
singularity exec --bind /hpcnfs/ ${BIOINFO_SING} ${COMMAND}
${COMMAND_MV}" | qsub -o ${ID}.STDOUT -e ${ID}.STDERR -l select=1:ncpus=${NCPU}:mem=${MEM2} -N ${ID}

exit

