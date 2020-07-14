#!/bin/bash

mkdir -p  train/gwas_results/mlma

freq=${1} # can either be "cm.200k" or "re.all"
pheno_ix=${2}
chr=${3}
out=${4}

#without PCs as fixed effects
gcta_b --bfile train/genotypes/beds/genos_gridt100_l1e7_ss750_m0.05_chr1_20.rmdup.train.snps.all \
--mlma \
--grm train/genotypes/grms/genos_gridt100_l1e7_ss750_m0.05_chr1_20.rmdup.train.${freq} \
--pheno phenotypes/noge/pheno_gridt100_noge_s9k.train.1.txt \
--threads 12 \
--chr ${chr} \
--out train/gwas_results/mlma/${out}.${chr}.gonly \
--mpheno ${pheno_ix} \

#with PCs as fixed effects
gcta_b --bfile train/genotypes/beds/genos_gridt100_l1e7_ss750_m0.05_chr1_20.rmdup.train.snps.all \
--mlma \
--grm train/genotypes/grms/genos_gridt100_l1e7_ss750_m0.05_chr1_20.rmdup.train.${freq} \
--pheno phenotypes/noge/pheno_gridt100_noge_s9k.train.1.txt \
--threads 12 \
--chr ${chr} \
--out train/gwas_results/mlma/${out}.${chr}.gwtpcs \
--mpheno ${pheno_ix} \
--qcovar train/genotypes/genos_gridt100_l1e7_ss750_m0.05_chr1_20.rmdup.train.${freq}.eigenvec
