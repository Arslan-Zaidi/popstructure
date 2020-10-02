#!/bin/bash

rep=${1}

source /home/aazaidi/anaconda3/etc/profile.d/conda.sh
conda activate rpkgs

mkdir -p test/phenotypes/ge
mkdir -p test/gwas_results/fixed_effects/ge/

echo "simulating phenotypes"
#use these genetic values to generate phenotypes
simphenotype_ge.R \
test/gvalue/genos_grid_d36_m0.05_s500_t100.rmdup.test.${rep}.gvalue.sscore \
iid_test.txt \
test/phenotypes/ge/pheno_gridt100_ge_s9k.test.${rep}.txt \
${rep}

echo "running GWAS on heritable phenotypes"
## genetic phenotypes
#1. no correction
gwas.sh \
test/genotypes/genos_gridt100_l1e7_ss750_m0.05_chr1_20.rmdup.test \
test/phenotypes/ge/pheno_gridt100_ge_s9k.test.${rep}.txt \
test/gwas_results/fixed_effects/ge/gwas_gridt100_test.ge.${rep}.pcs0

#2. common pca
gwas.sh \
test/genotypes/genos_gridt100_l1e7_ss750_m0.05_chr1_20.rmdup.test \
test/phenotypes/ge/pheno_gridt100_ge_s9k.test.${rep}.txt \
test/gwas_results/fixed_effects/ge/gwas_gridt100_test.ge.${rep}.cm \
test/genotypes/genos_gridt100_l1e7_ss750_m0.05_chr1_20.rmdup.test.cm.200k.pca.eigenvec

#3. rare pca
gwas.sh \
test/genotypes/genos_gridt100_l1e7_ss750_m0.05_chr1_20.rmdup.test \
test/phenotypes/ge/pheno_gridt100_ge_s9k.test.${rep}.txt \
test/gwas_results/fixed_effects/ge/gwas_gridt100_test.ge.${rep}.re \
test/genotypes/genos_gridt100_l1e7_ss750_m0.05_chr1_20.rmdup.test.re.all.pca.eigenvec

#4. both common and rare
# gwas.sh \
# test/genotypes/genos_gridt100_l1e7_ss750_m0.05_chr1_20.rmdup.test \
# phenotypes/ge/pheno_gridt100_ge_s9k.test.${rep}.txt \
# test/gwas_results/fixed_effects/ge/gwas_gridt100_test.ge.${rep}.cmre \
# test/genotypes/genos_gridt100_l1e7_ss750_m0.05_chr1_20.rmdup.test.cmre.pca.eigenvec

echo "Zipping glm results"
gzip -f test/gwas_results/fixed_effects/ge/*.${rep}.*.glm.linear
