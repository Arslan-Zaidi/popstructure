#!/bin/bash

mkdir -p  train/gwas_results/mlma/ge

freq=${1} # can either be "cm.200k" or "re.all"
phenotype=${2}
pheno_ix=${3}
chr=${4}
out=${5}

#without PCs as fixed effects
# gcta_b --bfile train/genotypes/beds/genos_gridt100_l1e7_ss750_m0.05_chr1_20.rmdup.train.snps.all \
# --mlma \
# --grm train/genotypes/grms/genos_gridt100_l1e7_ss750_m0.05_chr1_20.rmdup.train.${freq} \
# --pheno phenotypes/noge/pheno_gridt100_noge_s9k.train.1.txt \
# --threads 12 \
# --chr ${chr} \
# --out train/gwas_results/mlma/${out}.${chr}.gonly \
# --mpheno ${pheno_ix} \

#with PCs as fixed effects
gcta_b --bfile train/genotypes/beds/genos_gridt100_l1e7_ss750_m0.05_chr1_20.rmdup.train.snps.all \
--mlma \
--grm train/genotypes/grms/genos_gridt100_l1e7_ss750_m0.05_chr1_20.rmdup.train.${freq} \
--pheno ${phenotype} \
--threads 12 \
--chr ${chr} \
--out ${out} \
--mpheno ${pheno_ix} \
--qcovar train/genotypes/genos_gridt100_l1e7_ss750_m0.05_chr1_20.rmdup.train.${freq}.eigenvec \
--mlma-subtract-grm train/genotypes/grms/grm_gridt100_l1e7_ss750_m0.05_chr${chr}.rmdup.train.${freq}
