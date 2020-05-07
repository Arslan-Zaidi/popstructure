#!/bin/bash

#Run association between genotype and phenotype correcting for varying levels of correction

#exit and print usage if arguments not provided
[ $# -ne 4 ] && { echo "Usage: $0 <input genotype file path and prefix> <input phenotype file path> <output file prefix> <no. of Pcs to correct for>"; exit 1; }

input_geno=${1}
input_pheno=${2}
output=${3}
pcs=${4}

#calculate total length of the genome
#chrlength=`expr $nchr \* 10`
pc_cols=$(($pcs+2))

mkdir -p gwas_results

if [ $pcs == 0 ];
then
  plink2 --bfile ${input_geno} \
  --read-freq ${input_geno}.frq.afreq \
  --linear hide-covar \
  --pheno ${input_pheno} \
  --allow-no-sex \
  --out gwas_results/${output}_pcs${pcs}.all

else
  plink2 --bfile ${input_geno} \
  --read-freq ${input_geno}.frq.afreq \
  --linear hide-covar \
  --pheno ${input_pheno} \
  --covar ${input_geno}.common.ldpruned.eigenvec \
  --covar-col-nums 3-${pc_cols} \
  --allow-no-sex \
  --out gwas_results/${output}_common.pcs${pcs}.all

  plink2 --bfile ${input_geno} \
  --read-freq ${input_geno}.frq.afreq \
  --linear hide-covar \
  --pheno ${input_pheno} \
  --covar ${input_geno}.rare.ldpruned.eigenvec \
  --covar-col-nums 3-${pc_cols} \
  --allow-no-sex \
  --out gwas_results/${output}_rare.pcs${pcs}.all

fi
