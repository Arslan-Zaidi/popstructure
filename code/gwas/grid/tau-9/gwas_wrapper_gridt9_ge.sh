#!/bin/bash

rep=${1}

source /home/aazaidi/anaconda3/etc/profile.d/conda.sh
conda activate rpkgs

mkdir -p phenotypes/ge


echo "simulating heritable phenotypes"
#simulate genetic phenotype
bash scripts/simphenotype/simphenotype_ge_wrapper.sh \
train/genotypes/genos_grid_d36_m0.07_s500_t9_chr1_20.rmdup.train \
genos_grid_d36_m0.07_s500_t9.train.pop \
phenotypes/ge/pheno_gridt9_ge_s9k.train.${rep}.txt \
${rep}

echo "running GWAS on heritable phenotypes"
## genetic phenotypes
#1. no correction
bash scripts/gwas/gwas.sh \
train/genotypes/genos_grid_d36_m0.07_s500_t9_chr1_20.rmdup.train.snps \
phenotypes/ge/pheno_gridt9_ge_s9k.train.${rep}.txt \
train/gwas_results/fixed_effects/ge/gwas_gridt9_train.ge.${rep}.pcs0

#2. common pca
bash scripts/gwas/gwas.sh \
train/genotypes/genos_grid_d36_m0.07_s500_t9_chr1_20.rmdup.train.snps \
phenotypes/ge/pheno_gridt9_ge_s9k.train.${rep}.txt \
train/gwas_results/fixed_effects/ge/gwas_gridt9_train.ge.${rep}.cm \
train/genotypes/genos_grid_d36_m0.07_s500_t9_chr1_20.rmdup.train.snps.cm.200k.pca.eigenvec

#3. rare pca
bash scripts/gwas/gwas.sh \
train/genotypes/genos_grid_d36_m0.07_s500_t9_chr1_20.rmdup.train.snps \
phenotypes/ge/pheno_gridt9_ge_s9k.train.${rep}.txt \
train/gwas_results/fixed_effects/ge/gwas_gridt9_train.ge.${rep}.re \
train/genotypes/genos_grid_d36_m0.07_s500_t9_chr1_20.rmdup.train.snps.re.1M.pca.eigenvec

#4. both common and rare
bash scripts/gwas/gwas.sh \
train/genotypes/genos_grid_d36_m0.07_s500_t9_chr1_20.rmdup.train.snps \
phenotypes/ge/pheno_gridt9_ge_s9k.train.${rep}.txt \
train/gwas_results/fixed_effects/ge/gwas_gridt9_train.ge.${rep}.cmre \
train/genotypes/genos_grid_d36_m0.07_s500_t9_chr1_20.rmdup.train.snps.cmre.pca.eigenvec

echo "Zipping glm results"
gzip -f train/gwas_results/fixed_effects/ge/*.${rep}.*.glm.linear
