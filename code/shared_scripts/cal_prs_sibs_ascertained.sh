
#!/bin/bash

phenotype=${1}
rep=${2}

mkdir -p test/sibs_prs/${phenotype}
mkdir -p sibs/betas

echo "get effect sizes of SNPs ascertained in GWAS"

Rscript scripts/prs/ascertain_effects.R \
sibs/gwas_results/e${rep}/gwas_gridt100_l1e7_ss750_m0.05_chr1_20_e${rep}_n9000_sibs_assort.${phenotype}.glm.linear.gz \
train/betas/est_effects.${rep}.${phenotype}.nc.betas \
sibs/betas/est_effects.p${phenotype}.e${rep}.asc.nc.betas

echo "constructing PRS in test cases"

plink2 --pfile \
test/genotypes/genos_gridt100_l1e7_ss750_m0.05_chr1_20.rmdup.test \
--read-freq test/genotypes/genos_gridt100_l1e7_ss750_m0.05_chr1_20.rmdup.test.frq.afreq \
--score sibs/betas/est_effects.p${phenotype}.e${rep}.asc.nc.betas cols=dosagesum,scoresums \
--out test/sibs_prs/${phenotype}/gridt100_prs_p${phenotype}.e${rep}.asc.nc \
--score-col-nums 3,4,5,6
