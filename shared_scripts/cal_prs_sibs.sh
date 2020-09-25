
#!/bin/bash

phenotype=${1}
rep=${2}

mkdir -p train/prs/${phenotype}
mkdir -p test/sibs/betas

echo "picking significant associations"

Rscript scripts/prs/clump_sibs.R \
train/genotypes/genos_gridt100_l1e7_ss750_m0.05_chr1_20.rmdup.train.${rep}.thinned_100kb.effects \
sibs/gwas_results/e${rep}/gwas_gridt100_l1e7_ss750_m0.05_chr1_20_e${rep}_n9000_sibs_assort.sharp.glm.linear.gz \
${phenotype} \
sibs/betas/est_effects.p${phenotype}.e${rep}

echo "constructing PRS in test cases"

plink2 --pfile \
test/genotypes/genos_gridt100_l1e7_ss750_m0.05_chr1_20.rmdup.test \
--read-freq test/genotypes/genos_gridt100_l1e7_ss750_m0.05_chr1_20.rmdup.test.frq.afreq \
--score sibs/betas/est_effects.p${phenotype}.e${rep}.all.betas cols=dosagesum,scoresums \
--out test/sibs_prs/${phenotype}/gridt100_prs_p${phenotype}.e${rep}.all \
--score-col-nums 3,4,5
