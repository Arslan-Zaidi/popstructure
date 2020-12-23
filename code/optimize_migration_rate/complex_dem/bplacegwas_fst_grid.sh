#!/bin/bash

#script to calculate Fst and run GWAS on birthplace given tau (time to panmixia) and m (migration rate among demes)

[ $# -ne 3 ] && { echo "Usage: $0 <tau {100,-9}> <sample_size> <migrate>"; exit 1; }

tau=${1} # time to panmixia
ss=${2} # sample size for EACH deme. NOT the total sample size
m=${3} #one-way migration rate among demes

mkdir -p output/tau${tau}/m${m}
cd output/tau${tau}/m${m}


echo "calculating Fst across demes"
plink2 --bfile genos_grid_d36_m${m}_s${ss}_t${tau} \
--freq \
--out genos_grid_d36_m${m}_s${ss}_t${tau}.frq

plink2 --bfile genos_grid_d36_m${m}_s${ss}_t${tau} \
--read-freq genos_grid_d36_m${m}_s${ss}_t${tau}.frq.afreq \
--maf 0.05 \
--indep-pairwise 100 10 0.1 \
--out genos_grid_d36_m${m}_s${ss}_t${tau}.pruning \

plink --bfile genos_grid_d36_m${m}_s${ss}_t${tau} \
--read-freq genos_grid_d36_m${m}_s${ss}_t${tau}.frq.afreq \
--extract genos_grid_d36_m${m}_s${ss}_t${tau}.pruning.prune.in \
--fst \
--within genos_grid_d36_m${m}_s${ss}_t${tau}.pop \
--out genos_grid_d36_m${m}_s${ss}_t${tau}.fst

echo "carrying out GWAS on deme identity"
plink2 --bfile genos_grid_d36_m${m}_s${ss}_t${tau} \
--read-freq genos_grid_d36_m${m}_s${ss}_t${tau}.frq.afreq \
--linear \
--pheno genos_grid_d36_m${m}_s${ss}_t${tau}.pop \
--pheno-name Latitude,Longitude \
--out genos_grid_d36_m${m}_s${ss}_t${tau}.bplace
