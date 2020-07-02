
#!/bin/bash

phenotype=${1}
mating=${2}
rep=${3}

mkdir -p train/prs/${phenotype}
mkdir -p test/sibs/betas

echo "get effect sizes of SNPs ascertained in GWAS"

Rscript scripts/prs/ascertain_effects.R \
test/sibs/gwas_results/gwas_complex_l1e7_ss500_m0.08_chr1_20_e${rep}_sibs_${mating}.${phenotype}.glm.linear.gz \
train/betas/est_effects.${rep}.${phenotype}.nc.betas \
test/sibs/betas/est_effects.p${phenotype}.m${mating}.e${rep}.asc.nc.betas

echo "constructing PRS in test cases"

plink2 --pfile \
train/genotypes/genos_complex_l1e7_ss500_m0.08_chr1_20.rmdup.train.snps \
--read-freq train/genotypes/genos_complex_l1e7_ss500_m0.08_chr1_20.rmdup.train.snps.frq.afreq \
--score test/sibs/betas/est_effects.p${phenotype}.m${mating}.e${rep}.asc.nc.betas cols=dosagesum,scoresums \
--out train/prs/${phenotype}/complexdem_prs_p${phenotype}.m${mating}.e${rep}.asc.nc \
--score-col-nums 3,4,5,6
