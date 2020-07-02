#!/bin/bash

rep=${1}

source /home/aazaidi/anaconda3/etc/profile.d/conda.sh
conda activate rpkgs

mkdir -p phenotypes/ge


echo "simulating heritable phenotypes"
#simulate genetic phenotype
bash scripts/simphenotype/simphenotype_ge_wrapper.sh \
train/genotypes/genos_gridt100_l1e7_ss750_m0.05_chr1_20.rmdup.train \
iid_train.txt \
phenotypes/ge/pheno_gridt100_ge_s9k.train.${rep}.txt \
${rep}

echo "running GWAS on heritable phenotypes"
## genetic phenotypes
#1. no correction
bash scripts/gwas/gwas.sh \
train/genotypes/genos_gridt100_l1e7_ss750_m0.05_chr1_20.rmdup.train \
phenotypes/ge/pheno_gridt100_ge_s9k.train.${rep}.txt \
train/gwas_results/fixed_effects/ge/gwas_gridt100_train.ge.${rep}.pcs0

#2. common pca
bash scripts/gwas/gwas.sh \
train/genotypes/genos_gridt100_l1e7_ss750_m0.05_chr1_20.rmdup.train \
phenotypes/ge/pheno_gridt100_ge_s9k.train.${rep}.txt \
train/gwas_results/fixed_effects/ge/gwas_gridt100_train.ge.${rep}.cm \
train/genotypes/genos_gridt100_l1e7_ss750_m0.05_chr1_20.rmdup.train.cm.200k.eigenvec

#3. rare pca
bash scripts/gwas/gwas.sh \
train/genotypes/genos_gridt100_l1e7_ss750_m0.05_chr1_20.rmdup.train \
phenotypes/ge/pheno_gridt100_ge_s9k.train.${rep}.txt \
train/gwas_results/fixed_effects/ge/gwas_gridt100_train.ge.${rep}.re \
train/genotypes/genos_gridt100_l1e7_ss750_m0.05_chr1_20.rmdup.train.re.all.eigenvec

#4. both common and rare
bash scripts/gwas/gwas.sh \
train/genotypes/genos_gridt100_l1e7_ss750_m0.05_chr1_20.rmdup.train \
phenotypes/ge/pheno_gridt100_ge_s9k.train.${rep}.txt \
train/gwas_results/fixed_effects/ge/gwas_gridt100_train.ge.${rep}.cmre \
train/genotypes/genos_gridt100_l1e7_ss750_m0.05_chr1_20.rmdup.train.cmre.eigenvec

echo "Zipping glm results"
gzip -f train/gwas_results/fixed_effects/ge/*.${rep}.*.glm.linear
