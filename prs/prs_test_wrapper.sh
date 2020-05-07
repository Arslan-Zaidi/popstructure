#!/bin/bash

phenotype=${1}


echo "picking significant associations"

Rscript clump_3.R \
train/genotypes/genos_complex_l1e7_ss500_m0.07_chr1_20.rmdup.train.snps.thinned_100kb.effects \
train/gwas_results/fixed_effects/ge/gwas_complex_train \
${phenotype} \
train/betas/est_effects.${phenotype}

echo "constructing PRS in test cases"

plink2 --pfile test/genotypes/genos_complex_l1e7_ss500_m0.07_chr1_20.rmdup.test \
--score train/betas/est_effects.${phenotype}.c.betas cols=dosagesum,scoresums \
--out test/prs/complexdem_prs_${phenotype}.c \
--score-col-nums 3,4,5


plink2 --pfile test/genotypes/genos_complex_l1e7_ss500_m0.07_chr1_20.rmdup.test \
--score train/betas/est_effects.${phenotype}.c.p.betas cols=dosagesum,scoresums \
--out test/prs/complexdem_prs_${phenotype}.c.p \
--score-col-nums 3,4,5

plink2 --pfile test/genotypes/genos_complex_l1e7_ss500_m0.07_chr1_20.rmdup.test \
--score train/betas/est_effects.${phenotype}.nc.betas cols=dosagesum,scoresums \
--out test/prs/complexdem_prs_${phenotype}.nc \
--score-col-nums 3,4,5
