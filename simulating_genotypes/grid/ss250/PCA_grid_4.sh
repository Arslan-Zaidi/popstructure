#!/bin/bash

#script will take in merged plink file and performs cleanup and pca

#exit and print usage if arguments not provided
[ $# -ne 4 ] && { echo "Usage: $0 <tau {100,-9}> <chr_number> <sample_size> <migration rate>"; exit 1; }

#enter parameters
#these are not directly used in anly analyses etc. but just to generate the filename
tau=${1} # time to panmixia
nchr=${2} # no. of chromosomes to simulate
ss=${3} # sample size for EACH deme. NOT the total sample size
m=${4} # migration rate. would be 0.05 for tau=100 and 0.07 for tau=-9

cd genotypes/tau${tau}

#calculate total length of the genome
chrlength=`expr $nchr \* 10`

#specify name of input file
input=genos_grid_d36_m${m}_s${ss}_t${tau}_gl${chrlength}


echo "filtering for MAF"
plink2 --bfile $input \
--freq \
--out ${input}.frq

#isolate genotypes for common variants
plink2 --bfile $input \
--read-freq ${input}.frq.afreq \
--maf 0.05 \
--make-bed \
--out ${input}.common | tee -a ${input}.process.log

#isolate genotypes for rare variants (doubletons - 4)
plink2 --bfile $input \
--read-freq ${input}.frq.afreq \
--min-ac 2 \
--max-mac 4 \
--make-bed \
--out ${input}.rare | tee -a ${input}.process.log

echo "LDpruning"

plink2 --bfile ${input}.common \
--read-freq ${input}.frq.afreq \
--indep-pairwise 100 10 0.1 \
--out tmp1 | tee -a ${input}.process.log

plink2 --bfile ${input}.rare \
--read-freq ${input}.frq.afreq \
--indep-pairwise 100 10 0.1 \
--out tmp2 | tee -a ${input}.process.log


plink2 --bfile ${input}.common \
--read-freq ${input}.frq.afreq \
--extract tmp1.prune.in \
--make-bed \
--out ${input}.common.ldpruned | tee -a ${input}.process.log

plink2 --bfile ${input}.rare \
--read-freq ${input}.frq.afreq \
--extract tmp2.prune.in \
--make-bed \
--out ${input}.rare.ldpruned | tee -a ${input}.process.log

echo "carrying out PCA"

plink2 --bfile ${input}.common.ldpruned \
--read-freq ${input}.frq.afreq \
--pca 200 \
--out ${input}.common.ldpruned | tee -a ${input}.process.log

plink2 --bfile ${input}.rare.ldpruned \
--read-freq ${input}.frq.afreq \
--pca 200 \
--out ${input}.rare.ldpruned | tee -a ${input}.process.log

echo "cleaning up"
rm tmp1*
rm tmp2*
rm *.nosex
