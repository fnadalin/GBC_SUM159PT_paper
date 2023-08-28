# exp1
DIR="exp_1/1B_DNA"
OUT="${DIR}/exp_rates"
[ -d "${OUT}" ] || mkdir -p "${OUT}"
for i in D13 D17 D24
do
    Rscript COMMAND_compute_expansion_rates.R "${DIR}/D0/out_d1_f1/groups.tsv" "${DIR}/$i/out_d1_f1/groups.tsv" "${OUT}/D0-$i.tsv"
done

# exp2A
DIR="exp_2A/2A_DNA"
OUT="${DIR}/exp_rates"
[ -d "${OUT}" ] || mkdir -p "${OUT}"
for i in NT_D13 NT_D17 NT_D24 PACLI_D13-A PACLI_D13-B PACLI_D13-C PACLI_D17-A PACLI_D17-B PACLI_D24-A PACLI_D24-B PP-D20-A PP-D20-B  
do
    Rscript COMMAND_compute_expansion_rates.R "${DIR}/P0/out_d1_f1/groups.tsv" "${DIR}/$i/out_d1_f1/groups.tsv" "${OUT}/P0-$i.tsv"
done

# exp2B
DIR="exp_2B/2B_DNA"
OUT="${DIR}/exp_rates"
[ -d "${OUT}" ] || mkdir -p "${OUT}"
for i in Paclid17 Paclid20 
do
    Rscript COMMAND_compute_expansion_rates.R "${DIR}/P50R-P0/out_d1_f1/groups.tsv" "${DIR}/P50R-$i/out_d1_f1/groups.tsv" "${OUT}/P0-$i.tsv"
done



exit

