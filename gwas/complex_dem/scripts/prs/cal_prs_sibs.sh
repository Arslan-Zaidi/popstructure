
#!/bin/bash

phenotype=${1}
mating=${2}
rep=${3}

mkdir -p train/prs/${phenotype}
mkdir -p test/sibs/betas

echo "picking significant associations"

Rscript scripts/prs/clump_sibs.R \
train/genotypes/genos_complex_l1e7_ss500_m0.08_chr1_20.rmdup.train.${rep}.thinned_100kb.effects \
test/sibs/gwas_results/gwas_complex_l1e7_ss500_m0.08_chr1_20_e${rep}_sibs_${mating}.${phenotype}.glm.linear.gz \
${phenotype} \
test/sibs/betas/est_effects.p${phenotype}.m${mating}.e${rep}

echo "constructing PRS in test cases"

plink2 --pfile \
train/genotypes/genos_complex_l1e7_ss500_m0.08_chr1_20.rmdup.train.snps \
--read-freq train/genotypes/genos_complex_l1e7_ss500_m0.08_chr1_20.rmdup.train.snps.frq.afreq \
--score test/sibs/betas/est_effects.p${phenotype}.m${mating}.e${rep}.all.betas cols=dosagesum,scoresums \
--out train/prs/${phenotype}/complexdem_prs_p${phenotype}.m${mating}.e${rep}.all \
--score-col-nums 3,4,5

#
# plink2 --pfile \
# train/genotypes/genos_complex_l1e7_ss500_m0.08_chr1_20.rmdup.train.snps \
# --read-freq train/genotypes/genos_complex_l1e7_ss500_m0.08_chr1_20.rmdup.train.snps.frq.afreq \
# --score test/sibs/betas/est_effects.p${phenotype}.m${mating}.e${rep}.c.p.betas cols=dosagesum,scoresums \
# --out train/prs/${phenotype}/complexdem_prs_p${phenotype}.m${mating}.e${rep}.c.p \
# --score-col-nums 3
#
# plink2 --pfile \
# train/genotypes/genos_complex_l1e7_ss500_m0.08_chr1_20.rmdup.train.snps \
# --read-freq train/genotypes/genos_complex_l1e7_ss500_m0.08_chr1_20.rmdup.train.snps.frq.afreq \
# --score test/sibs/betas/est_effects.p${phenotype}.m${mating}.e${rep}.nc.betas cols=dosagesum,scoresums \
# --out train/prs/${phenotype}/complexdem_prs_p${phenotype}.m${mating}.e${rep}.nc \
# --score-col-nums 3
