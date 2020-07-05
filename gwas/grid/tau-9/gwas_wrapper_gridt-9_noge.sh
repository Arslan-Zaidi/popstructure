#!/bin/bash

rep=${1}

source /home/aazaidi/anaconda3/etc/profile.d/conda.sh
conda activate rpkgs

mkdir -p phenotypes/noge

echo "simulating non-heritable phenotypes"
Rscript scripts/simphenotype/simphenotype_noge.R \
genos_grid_d36_m0.07_s500_t9.train.pop \
phenotypes/noge/pheno_gridt9_noge_s9k.train.${rep}.txt \
${rep}

echo "running GWAS on non-heritable phenotypes"
## genetic phenotypes
#1. no correction
bash scripts/gwas/gwas.sh \
train/genotypes/genos_grid_d36_m0.07_s500_t9_chr1_20.rmdup.train.snps \
phenotypes/noge/pheno_gridt9_noge_s9k.train.${rep}.txt \
train/gwas_results/fixed_effects/noge/gwas_gridt9_train.noge.${rep}.pcs0

#2. common pca
bash scripts/gwas/gwas.sh \
train/genotypes/genos_grid_d36_m0.07_s500_t9_chr1_20.rmdup.train.snps \
phenotypes/noge/pheno_gridt9_noge_s9k.train.${rep}.txt \
train/gwas_results/fixed_effects/noge/gwas_gridt9_train.noge.${rep}.cm \
train/genotypes/genos_grid_d36_m0.07_s500_t9_chr1_20.rmdup.train.snps.cm.200k.pca.eigenvec

#3. rare pca
bash scripts/gwas/gwas.sh \
train/genotypes/genos_grid_d36_m0.07_s500_t9_chr1_20.rmdup.train.snps \
phenotypes/noge/pheno_gridt9_noge_s9k.train.${rep}.txt \
train/gwas_results/fixed_effects/noge/gwas_gridt9_train.noge.${rep}.re \
train/genotypes/genos_grid_d36_m0.07_s500_t9_chr1_20.rmdup.train.snps.re.1M.pca.eigenvec


echo "Zipping glm results"
gzip -f train/gwas_results/fixed_effects/noge/*.${rep}.*.glm.linear
