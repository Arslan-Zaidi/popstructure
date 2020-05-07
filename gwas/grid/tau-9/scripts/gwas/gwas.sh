#!/bin/bash

#script generates genetic effects and phenotypes with different patterns of stratification
#
if [ "$#" -lt 3  ]; then
    echo "usage: bash gwas_ge.sh <path to pgen genotype file - prefix> <path to phenotype file> <path to output file> <eigenvector file - optional>"
    exit 1
fi


geno_file_prefix=${1}
phenotype_file=${2}
output_file_prefix=${3}
eigenvec_file=${4}


if [ "$#" == 3 ];
then
  plink2 --pfile ${geno_file_prefix} \
  --mac 1 \
  --glm hide-covar \
  --pheno ${phenotype_file} \
  --out ${output_file_prefix}
fi

if [ "$#" == 4 ];
then
  plink2 --pfile ${geno_file_prefix} \
  --mac 1 \
  --glm hide-covar \
  --pheno ${phenotype_file} \
  --covar ${eigenvec_file} \
  --covar-col-nums 3-102 \
  --out ${output_file_prefix}
fi
