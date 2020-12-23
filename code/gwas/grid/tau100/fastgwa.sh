#!/bin/bash


freq=${1} # can be 'cm' or 're'
phenotype=${2}
pheno_ix=${3}
out=${4}
grm=train/grms/genos_grid_d36_m0.05_s500_t100.rmdup.train.${freq}pruned.grm.sparse

#without PCs as fixed effects
gcta_b --bfile train/genotypes/genos_grid_d36_m0.05_s500_t100.rmdup.train.b \
--grm-sparse ${grm} \
--fastGWA-mlm \
--pheno ${phenotype} \
--threads 12 \
--out train/gwas_results/mlm/${out}.gonly \
--mpheno ${pheno_ix} \
--maf 5e-05

#with PCs as fixed effects
gcta_b --bfile train/genotypes/genos_grid_d36_m0.05_s500_t100.rmdup.train.b \
--grm-sparse ${grm} \
--fastGWA-mlm \
--pheno ${phenotype} \
--threads 12 \
--out train/gwas_results/mlm/${out}.gwtpcs \
--mpheno ${pheno_ix} \
--maf 5e-05 \
--qcovar train/plinkPCA/genos_grid_d36_m0.05_s500_t100.rmdup.train.${freq}pruned.pca.eigenvec
