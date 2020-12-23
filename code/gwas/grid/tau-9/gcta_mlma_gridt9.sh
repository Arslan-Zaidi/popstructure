#!/bin/bash

mkdir -p  train/gwas_results/mlma

freq=${1} # can either be "cm.200k" or "re.1M"
pheno_ix=${2}
chr=${3}
out=${4}

#without PCs as fixed effects
# gcta_b --bfile train/genotypes/beds/genos_grid_d36_m0.07_s500_t9_chr1_20.rmdup.train.b${chr} \
# --mlma \
# --grm train/genotypes/grms/genos_grid_d36_m0.07_s500_t9_chr1_20.rmdup.train.snps.${freq}.pca.eigenvec \
# --pheno phenotypes/noge/pheno_gridt100_noge_s9k.train.1.txt \
# --threads 12 \
# --out train/gwas_results/mlma/${out}.${chr}.gonly \
# --mpheno ${pheno_ix} \

#with PCs as fixed effects
gcta_b --bfile train/genotypes/beds/genos_grid_d36_m0.07_s500_t9_chr1_20.rmdup.train.b${chr} \
--mlma \
--grm train/genotypes/grms/genos_grid_d36_m0.07_s500_t9_chr1_20.rmdup.train.snps.${freq} \
--pheno phenotypes/noge/pheno_gridt9_noge_s9k.train.1.txt \
--threads 12 \
--out train/gwas_results/mlma/${out}.${chr}.gwtpcs \
--mpheno ${pheno_ix} \
--qcovar train/genotypes/genos_grid_d36_m0.07_s500_t9_chr1_20.rmdup.train.snps.${freq}.pca.eigenvec
