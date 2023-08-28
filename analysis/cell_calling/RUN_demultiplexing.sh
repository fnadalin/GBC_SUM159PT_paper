RPATH="$HOME/.local/R-3.6.3/bin"
RPATH2="/hpcnfs/data/FN/R403/bin"
PREAMBLE="source activate /hpcnfs/data/FN/R403 ; export R_LIBS_USER=/hpcnfs/data/FN/R403/lib/R/library"
SCRIPTDIR=$( cd "../../scripts/demultiplexing/" ; pwd )
DPX_SCRIPT_DIR="$HOME/git/r_scripts/clones_analysis"

DIR="1D_GEX/"
OUT="1D_GEX_demux/"

[[ -d "${OUT}" ]] || mkdir "${OUT}"

SAMPLES=( "P0" "A1" "A2" "B1" "B2" )
MBC=( "MBC_P0.txt" "MBC_TC.txt" "MBC_TC.txt" "MBC_TC.txt" "MBC_TC.txt" )

CORES="1"
MEM1="64G"
MEM2="20G"

i=0
depend="afterok"
for s in ${SAMPLES[@]}
do
    DPXDIR="${OUT}/$s/aln"
    [[ -d "${DPXDIR}" ]] || mkdir -p "${DPXDIR}"

    # copy input files in the MBC 
    cp "${DIR}/$s/outs/filtered_feature_bc_matrix/barcodes.tsv.gz" "${DPXDIR}"/
    cp "libraries_exp_1D/lanes_MBC_${s}.csv" "${DPXDIR}/lanes_orig.csv"
    cp "${MBC[$i]}" "${DPXDIR}/MBC.txt" 

    # split FASTQ files into chunks, run MBC calling (with deMULTIplex), remove files 
    ID="${s}_MBC_align"
    COMMAND="bash ${SCRIPTDIR}/preprocessing.sh \"${DPXDIR}\" 
             ${RPATH2}/Rscript ${SCRIPTDIR}/align_MBC.R \"${DPXDIR}\"
             bash ${SCRIPTDIR}/postprocessing.sh \"${DPXDIR}\""
    JOBID=$(echo "cd $(pwd) ; ${PREAMBLE} ; ${COMMAND}" -o | qsub -N ${ID} -o ${ID}.STDOUT -e ${ID}.STDERR -l select=1:ncpus=${CORES}:mem=${MEM1})

    INDIR="${DIR}/$s/outs/filtered_feature_bc_matrix"
    BARTABLE="${DPXDIR}/out/bar.table.tsv"
    OUTDIR="${OUT}/$s"

    # create gene-GBC-MBC matrix
    ID="${s}_MBC_matrix"
    COMMAND="${RPATH}/Rscript ${SCRIPTDIR}/build_mbc_matrix.R \"${INDIR}\" \"${BARTABLE}\" \"${OUTDIR}\""
    JOBID1=$(echo "cd $(pwd); ${COMMAND}" | qsub -N ${ID} -o ${ID}.STDOUT -e ${ID}.STDERR -l select=1:ncpus=1:mem=${MEM2} -W depend=afterok:${JOBID})

    OUTDIR_MBC="${OUT}/$s/MBC_call"
    [[ -d "${OUTDIR_MBC}" ]] || mkdir -p "${OUTDIR_MBC}"

    # demultiplex
    ID="${s}_demux"
    COMMAND="bash ${SCRIPTDIR}/demultiplexing_step1.sh \"${DPX_SCRIPT_DIR}\" \"${OUTDIR}\" \"${OUTDIR_MBC}\"
             bash ${SCRIPTDIR}/demultiplexing_step2.sh \"${DPX_SCRIPT_DIR}\" \"${OUTDIR_MBC}\"
             ${RPATH}/Rscript ${SCRIPTDIR}/demultiplexing_step3.R \"${OUTDIR_MBC}\" \"${OUT}/$s/MBC_calling_$s.tsv\""
    JOBID2=$(echo "cd $(pwd); ${COMMAND}" | qsub -N ${ID} -o ${ID}.STDOUT -e ${ID}.STDERR -l select=1:ncpus=1:mem=${MEM2} -W depend=afterok:${JOBID1})

    OUTDIR_DEMUX="${OUT}/$s/demux"
    [[ -d "${OUTDIR_DEMUX}" ]] || mkdir -p "${OUTDIR_DEMUX}"

    # create demultiplexed count matrices
    ID="${s}_demux_matrix"
    COMMAND="${RPATH}/Rscript ${SCRIPTDIR}/demultiplexing_step4.R \"${OUTDIR}\" \"${OUTDIR}/MBC_calling_${s}.tsv\" \"${OUTDIR_DEMUX}\""
    JOBID3=$(echo "cd $(pwd); ${COMMAND}" | qsub -N ${ID} -o ${ID}.STDOUT -e ${ID}.STDERR -l select=1:ncpus=1:mem=${MEM2} -W depend=afterok:${JOBID2})

    depend=":${JOBID3}"

    i=$((i+1))
done

# create count matrices per samples
COMMAND="${RPATH}/Rscript ${SCRIPTDIR}/demultiplexing_aggregate.R demux_info_exp_1D.tsv \"${OUT}/demux\""
echo "cd $(pwd); ${COMMAND}" | qsub -N demux_aggr -o demux_aggr.STDOUT -e demux_aggr.STDERR -l select=1:ncpus=1:mem=${MEM2} -W depend=${depend}

exit

