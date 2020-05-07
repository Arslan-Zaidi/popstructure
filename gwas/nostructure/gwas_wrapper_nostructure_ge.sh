#!/bin/bash

rep=${1}

source /home/aazaidi/anaconda3/etc/profile.d/conda.sh
conda activate rpkgs

mkdir -p phenotypes/ge

echo "simulating heritable phenotypes"
#simulate genetic phenotype
bash scripts/simphenotype/simphenotype_ge_wrapper.sh \
train/genotypes/genos_nostr_s18k_l1e7_chr1_20.rmdup.train \
genos_complex_l1e7_ss500_m0.07.train.pop \
phenotypes/ge/pheno_nostructure_ge_s9k.train.${rep}.txt \
${rep}

echo "running GWAS on heritable phenotypes"
## genetic phenotypes
#1. no correction
bash scripts/gwas/gwas.sh \
train/genotypes/genos_nostr_s18k_l1e7_chr1_20.rmdup.train \
phenotypes/ge/pheno_nostructure_ge_s9k.train.${rep}.txt \
train/gwas_results/fixed_effects/ge/gwas_nostr_s9k_train.ge.${rep}.pcs0

#2. common pca
bash scripts/gwas/gwas.sh \
train/genotypes/genos_nostr_s18k_l1e7_chr1_20.rmdup.train \
phenotypes/ge/pheno_nostructure_ge_s9k.train.${rep}.txt \
train/gwas_results/fixed_effects/ge/gwas_nostr_s9k_train.ge.${rep}.cm \
train/genotypes/genos_nostr_s18k_l1e7_chr1_20.rmdup.train.cmpruned.pca.eigenvec

#3. rare pca
bash scripts/gwas/gwas.sh \
train/genotypes/genos_nostr_s18k_l1e7_chr1_20.rmdup.train \
phenotypes/ge/pheno_nostructure_ge_s9k.train.${rep}.txt \
train/gwas_results/fixed_effects/ge/gwas_nostr_s9k_train.ge.${rep}.re \
train/genotypes/genos_nostr_s18k_l1e7_chr1_20.rmdup.train.repruned.pca.eigenvec

#4. both common and rare
bash scripts/gwas/gwas.sh \
train/genotypes/genos_nostr_s18k_l1e7_chr1_20.rmdup.train \
phenotypes/ge/pheno_nostructure_ge_s9k.train.${rep}.txt \
train/gwas_results/fixed_effects/ge/gwas_nostr_s9k_train.ge.${rep}.cmre \
train/genotypes/genos_nostr_s18k_l1e7_chr1_20.rmdup.train.cmre.pca.eigenvec

echo "Zipping glm results"
gzip train/gwas_results/fixed_effects/ge/*.${rep}.*.glm.linear
