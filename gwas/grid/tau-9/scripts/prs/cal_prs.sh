
#!/bin/bash

phenotype=${1}
rep=${2}

echo "picking significant associations"

Rscript scripts/prs/clump.R \
train/genotypes/genos_grid_d36_m0.07_s500_t9_chr1_20.rmdup.train.subsampled.${rep}.thinned_100kb.effects \
train/gwas_results/fixed_effects/ge/gwas_gridt9_train.ge.${rep} \
${phenotype} \
train/betas/est_effects.${rep}.${phenotype}

echo "constructing PRS in test cases"

plink2 --pfile \
test/genotypes/genos_grid_d36_m0.07_s500_t9_chr1_20.rmdup.test.subsampled \
--read-freq test/genotypes/genos_grid_d36_m0.07_s500_t9_chr1_20.rmdup.test.subsampled.frq.afreq \
--score train/betas/est_effects.${rep}.${phenotype}.c.betas cols=dosagesum,scoresums \
--out test/prs/${phenotype}/gridt9_prs_${phenotype}.${rep}.c \
--score-col-nums 3,4,5,6

plink2 --pfile \
test/genotypes/genos_grid_d36_m0.07_s500_t9_chr1_20.rmdup.test.subsampled \
--read-freq test/genotypes/genos_grid_d36_m0.07_s500_t9_chr1_20.rmdup.test.subsampled.frq.afreq \
--score train/betas/est_effects.${rep}.${phenotype}.c.p.betas cols=dosagesum,scoresums \
--out test/prs/${phenotype}/gridt9_prs_${phenotype}.${rep}.c.p \
--score-col-nums 3,4,5,6

plink2 --pfile \
test/genotypes/genos_grid_d36_m0.07_s500_t9_chr1_20.rmdup.test.subsampled \
--read-freq test/genotypes/genos_grid_d36_m0.07_s500_t9_chr1_20.rmdup.test.subsampled.frq.afreq \
--score train/betas/est_effects.${rep}.${phenotype}.nc.betas cols=dosagesum,scoresums \
--out test/prs/${phenotype}/gridt9_prs_${phenotype}.${rep}.nc \
--score-col-nums 3,4,5,6
