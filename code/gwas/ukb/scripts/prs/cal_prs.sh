
#!/bin/bash

phenotype=${1}
rep=${2}
w=${3}

mkdir -p test/prs/${phenotype}

echo "picking significant associations"

Rscript scripts/prs/clump.R \
train/genotypes/genos_ukb_l1e7_ss500_m0.08_${w}_chr1_20.rmdup.train.${rep}.thinned_100kb.effects \
train/gwas_results/fixed_effects/ge/gwas_ukb_train.ge.${rep}.${w} \
${phenotype} \
train/betas/est_effects.${w}.${rep}.${phenotype}

echo "constructing PRS in test cases"

plink2 --pfile \
test/genotypes/genos_ukb_l1e7_ss500_m0.08_${w}_chr1_20.rmdup.test \
--score train/betas/est_effects.${w}.${rep}.${phenotype}.c.betas cols=dosagesum,scoresums \
--out test/prs/${phenotype}/ukbdem_prs_${phenotype}.${rep}.${w}.c \
--score-col-nums 3,4,5

plink2 --pfile \
test/genotypes/genos_ukb_l1e7_ss500_m0.08_${w}_chr1_20.rmdup.test \
--score train/betas/est_effects.${w}.${rep}.${phenotype}.c.p.betas cols=dosagesum,scoresums \
--out test/prs/${phenotype}/ukbdem_prs_${phenotype}.${rep}.${w}.c.p \
--score-col-nums 3,4,5

plink2 --pfile \
test/genotypes/genos_ukb_l1e7_ss500_m0.08_${w}_chr1_20.rmdup.test \
--score train/betas/est_effects.${w}.${rep}.${phenotype}.nc.betas cols=dosagesum,scoresums \
--out test/prs/${phenotype}/ukbdem_prs_${phenotype}.${rep}.${w}.nc \
--score-col-nums 3,4,5
