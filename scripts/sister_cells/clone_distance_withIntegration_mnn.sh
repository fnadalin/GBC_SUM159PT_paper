

if [[ $# -lt "7" ]]
then
    echo "bash $0 <obj1> <obj2> <cl_mode1> <cl_mode2> <hvg1> <hvg2> <outdir>"
    exit
fi

OBJ1=$1
OBJ2=$2
CL1=$3
CL2=$4
HVG1=$5
HVG2=$6
OUT=$7

DIR=$( cd $(dirname $0) ; pwd ) 

Rscript "${DIR}/clone_distance_withIntegration_mnn.R" "${OBJ1}" "${OBJ2}" ${CL1} ${CL2} "${HVG1}" "${HVG2}" "${OUT}"

exit

