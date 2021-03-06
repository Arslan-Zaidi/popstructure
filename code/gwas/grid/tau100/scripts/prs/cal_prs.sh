
#!/bin/bash

phenotype=${1}
rep=${2}

mkdir -p train/betas
mkdir -p test/prs/${phenotype}

echo "picking significant associations"

Rscript scripts/prs/clump.R \
train/genotypes/genos_gridt100_l1e7_ss750_m0.05_chr1_20.rmdup.train.${rep}.thinned_100kb.effects \
train/gwas_results/fixed_effects/ge/gwas_gridt100_train.ge.${rep} \
${phenotype} \
train/betas/est_effects.${rep}.${phenotype}

echo "constructing PRS in test cases"

plink2 --pfile \
test/genotypes/genos_gridt100_l1e7_ss750_m0.05_chr1_20.rmdup.test \
--score train/betas/est_effects.${rep}.${phenotype}.c.betas cols=dosagesum,scoresums \
--out test/prs/${phenotype}/gridt100_prs_${phenotype}.${rep}.c \
--score-col-nums 3,4,5,6

plink2 --pfile \
test/genotypes/genos_gridt100_l1e7_ss750_m0.05_chr1_20.rmdup.test \
--score train/betas/est_effects.${rep}.${phenotype}.c.p.betas cols=dosagesum,scoresums \
--out test/prs/${phenotype}/gridt100_prs_${phenotype}.${rep}.c.p \
--score-col-nums 3,4,5,6

plink2 --pfile \
test/genotypes/genos_gridt100_l1e7_ss750_m0.05_chr1_20.rmdup.test \
--score train/betas/est_effects.${rep}.${phenotype}.nc.betas cols=dosagesum,scoresums \
--out test/prs/${phenotype}/gridt100_prs_${phenotype}.${rep}.nc \
--score-col-nums 3,4,5,6
