

if [[ $# -lt "5" ]]
then
    echo "bash $0 <obj1> <obj2> <cl_mode1> <cl_mode2> <outdir>"
    exit
fi

OBJ1=$1
OBJ2=$2
CL1=$3
CL2=$4
OUTDIR=$5

DIR=$( cd $(dirname $0) ; pwd )

Rscript ${DIR}/cluster_clone_pairs.R "${OBJ1}" "${OBJ2}" ${CL1} ${CL2} "${OUTDIR}"
Rscript ${DIR}/plot_cluster_clone_pairs.R "${OUTDIR}"

exit

