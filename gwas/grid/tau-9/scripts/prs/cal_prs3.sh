
#!/bin/bash

phenotype=${1}
rep=${2}


mkdir -p test/prs3/${phenotype}
mkdir -p train/betas3/${phenotype}

echo "picking significant associations"

Rscript scripts/prs/clump3.R \
train/genotypes/genos_grid_d36_m0.07_s500_t9_chr1_20.rmdup.train.${rep}.thinned_100kb.effects \
train/gwas_results/fixed_effects/ge/gwas_gridt9_train.ge.${rep} \
${phenotype} \
train/betas3/est_effects.${rep}.${phenotype}

echo "constructing PRS in test cases"

plink2 --pfile \
test/genotypes/genos_grid_d36_m0.07_s500_t9_chr1_20.rmdup.test \
--read-freq test/genotypes/genos_grid_d36_m0.07_s500_t9_chr1_20.rmdup.test.frq.afreq \
--score train/betas3/est_effects.${rep}.${phenotype}.c.p.betas cols=dosagesum,scoresums \
--out test/prs3/${phenotype}/gridt9_prs_${phenotype}.${rep}.c.p \
--score-col-nums 3

plink2 --pfile \
test/genotypes/genos_grid_d36_m0.07_s500_t9_chr1_20.rmdup.test \
--read-freq test/genotypes/genos_grid_d36_m0.07_s500_t9_chr1_20.rmdup.test.frq.afreq \
--score train/betas3/est_effects.${rep}.${phenotype}.nc.betas cols=dosagesum,scoresums \
--out test/prs3/${phenotype}/gridt9_prs_${phenotype}.${rep}.nc \
--score-col-nums 3
