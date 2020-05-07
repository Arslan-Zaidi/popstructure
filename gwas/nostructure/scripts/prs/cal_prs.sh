#!/bin/bash

phenotype=${1}
rep=${2}

echo "picking significant associations"

Rscript scripts/prs/clump.R \
train/genotypes/genos_nostr_s18k_l1e7_chr1_20.rmdup.train.${rep}.thinned_100kb.effects \
train/gwas_results/fixed_effects/ge/gwas_nostr_s9k_train.ge.${rep} \
${phenotype} \
train/betas/est_effects.${rep}.${phenotype}

echo "constructing PRS in test cases"

plink2 --pfile \
test/genotypes/genos_nostr_s18k_l1e7_chr1_20.rmdup.test \
--score train/betas/est_effects.${rep}.${phenotype}.c.betas cols=dosagesum,scoresums \
--out test/prs/${phenotype}/nostructure_prs_${phenotype}.${rep}.c \
--score-col-nums 3,4,5

plink2 --pfile \
test/genotypes/genos_nostr_s18k_l1e7_chr1_20.rmdup.test \
--score train/betas/est_effects.${rep}.${phenotype}.c.p.betas cols=dosagesum,scoresums \
--out test/prs/${phenotype}/nostructure_prs_${phenotype}.${rep}.c.p \
--score-col-nums 3,4,5

plink2 --pfile \
test/genotypes/genos_nostr_s18k_l1e7_chr1_20.rmdup.test \
--score train/betas/est_effects.${rep}.${phenotype}.nc.betas cols=dosagesum,scoresums \
--out test/prs/${phenotype}/nostructure_prs_${phenotype}.${rep}.nc \
--score-col-nums 3,4,5
