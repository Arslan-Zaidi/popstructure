#!/bin/bash

# maf1 5e-4 1e-3
# maf2 1e-3 5e-3
# maf3 5e-3 1e-2
# maf4 1e-2 5e-2

label=${1}
maf=${2}
max_maf=${3}

plink2 \
--maf ${maf} \
--max-af ${max_maf} \
--pca 100 \
--pfile train/genotypes/genos_gridt100_l1e7_ss750_m0.05_chr1_20.rmdup.train \
--out train/genotypes/revisions/genos_gridt100_l1e7_ss750_m0.05_chr1_20.rmdup.train.${label}.pca

plink2 --pfile train/genotypes/genos_gridt100_l1e7_ss750_m0.05_chr1_20.rmdup.train \
--mac 1 \
--glm hide-covar \
--pheno phenotypes/noge/pheno_gridt100_noge_s9k.train.1.txt \
--pheno-name smooth,sharp \
--covar train/genotypes/revisions/genos_gridt100_l1e7_ss750_m0.05_chr1_20.rmdup.train.${label}.pca.eigenvec \
--covar-col-nums 3-102 \
--out train/gwas_results/fixed_effects/noge/revisions/gwas_gridt100_train.noge.1.${label}
