#!/bin/bash

#exit and print usage if arguments not provided
[ $# -ne 4 ] && { echo "Usage: $0 <tau {100,-9}> <chr_number> <sample_size> <migration rate>"; exit 1; }

tau=${1} # time to panmixia
chr=${2} # chromosome number (not the number of chromomes). should NOT exceed 22
ss=${3} # sample size for EACH deme. NOT the total sample size
m=${4}

# make directory to store results in
mkdir -p genotypes/tau${tau}

echo "generating vcf files"
#conda activate py3_clone --> source in terminal before running script

#generate genotypes given some tau for $nchr number of chromosomes, each of length 10Mb
#this will output the genotypes to csv format
python generate_genos_grid.py \
-s ${ss} \
-m ${m} \
-t ${tau} \
-d 36 \
-c ${chr} \
-o genotypes/tau${tau}/genos_grid_d36_m${m}_s${ss}_t${tau}

cd genotypes/tau${tau}

echo "converting to plink format"
plink2 --vcf genos_grid_d36_m${m}_s${ss}_t${tau}_chr${chr}.vcf \
--make-bed \
--double-id \
--allow-extra-chr \
--out genos_grid_d36_m${m}_s${ss}_t${tau}_chr${chr}

echo "cleaning up"

# mv genos_grid_d36_m1e-03_s${ss}_t${tau}_chr${chr}.bim tau${tau}
# mv genos_grid_d36_m1e-03_s${ss}_t${tau}_chr${chr}.bed tau${tau}
# mv genos_grid_d36_m1e-03_s${ss}_t${tau}_chr${chr}.fam tau${tau}
# rm genos_grid_d36_m1e-03_s${ss}_t${tau}_chr${chr}.vcf
# mv genos_grid_d36_m1e-03_s${ss}_t${tau}.pop tau${tau}
# rm genos_grid_d36_m1e-03_s${ss}_t${tau}_chr${chr}.log
#rm genos_grid_d36_m1e-03_s${ss}_t${tau}_chr${chr}.nosex
