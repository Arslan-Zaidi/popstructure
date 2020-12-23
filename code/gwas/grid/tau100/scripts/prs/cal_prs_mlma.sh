
#!/bin/bash

phenotype=${1}
freq=${2}
rep=${3}

mkdir -p test/prs_mlma/${phenotype}
mkdir -p train/betas_mlma

echo "picking significant associations"

Rscript scripts/prs/clump_mlma.R \
train/genotypes/genos_gridt100_l1e7_ss750_m0.05_chr1_20.rmdup.train.${rep}.thinned_100kb.effects \
train/gwas_results/mlma/ge/gwas_gridt100_train.ge.e${rep}.${freq}.${phenotype}.mlma.gz \
${phenotype} \
train/betas_mlma/est_effects.p${phenotype}.e${rep}.${freq}

echo "constructing PRS in test cases"

plink2 --pfile \
test/genotypes/genos_gridt100_l1e7_ss750_m0.05_chr1_20.rmdup.test \
--read-freq test/genotypes/genos_gridt100_l1e7_ss750_m0.05_chr1_20.rmdup.test.frq.afreq \
--score train/betas_mlma/est_effects.p${phenotype}.e${rep}.${freq}.all.betas cols=dosagesum,scoresums \
--out test/prs_mlma/${phenotype}/gridt100_prs_p${phenotype}.e${rep}.${freq} \
--score-col-nums 3,4,5
