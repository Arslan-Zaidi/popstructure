#!/bin/bash

rep=${1}
w=${2}

source /home/aazaidi/anaconda3/etc/profile.d/conda.sh
conda activate rpkgs

mkdir -p phenotypes/ge
mkdir -p train/gwas_results/fixed_effects/ge


echo "simulating heritable phenotypes"
#simulate genetic phenotype
bash scripts/simphenotype/simphenotype_ge_wrapper.sh \
train/genotypes/genos_ukb_l1e7_ss500_m0.08_${w}_chr1_20.rmdup.train \
popfiles/ukb_ss500_d35.${w}.pop \
phenotypes/ge/pheno_ukb_ge_s9k_${w}.train.${rep}.txt \
${rep}

echo "running GWAS on heritable phenotypes"
## genetic phenotypes
#1. no correction
bash scripts/gwas/gwas.sh \
train/genotypes/genos_ukb_l1e7_ss500_m0.08_${w}_chr1_20.rmdup.train \
phenotypes/ge/pheno_ukb_ge_s9k_${w}.train.${rep}.txt \
train/gwas_results/fixed_effects/ge/gwas_ukb_train.ge.${rep}.${w}.pcs0

#2. common pca
bash scripts/gwas/gwas.sh \
train/genotypes/genos_ukb_l1e7_ss500_m0.08_${w}_chr1_20.rmdup.train \
phenotypes/ge/pheno_ukb_ge_s9k_${w}.train.${rep}.txt \
train/gwas_results/fixed_effects/ge/gwas_ukb_train.ge.${rep}.${w}.cm \
train/genotypes/genos_ukb_l1e7_ss500_m0.08_${w}_chr1_20.rmdup.train.cm.200k.pca.eigenvec

#3. rare pca
bash scripts/gwas/gwas.sh \
train/genotypes/genos_ukb_l1e7_ss500_m0.08_${w}_chr1_20.rmdup.train \
phenotypes/ge/pheno_ukb_ge_s9k_${w}.train.${rep}.txt \
train/gwas_results/fixed_effects/ge/gwas_ukb_train.ge.${rep}.${w}.re \
train/genotypes/genos_ukb_l1e7_ss500_m0.08_${w}_chr1_20.rmdup.train.re.1M.pca.eigenvec

#4. cm + rare pca
# bash scripts/gwas/gwas.sh \
# train/genotypes/genos_ukb_l1e7_ss500_m0.08_${w}_chr1_20.rmdup.train \
# phenotypes/ge/pheno_ukb_ge_s9k_${w}.train.${rep}.txt \
# train/gwas_results/fixed_effects/ge/gwas_ukb_train.ge.${rep}.${w}.cmre \
# train/genotypes/genos_ukb_l1e7_ss500_m${m}_${w}_chr1_20.rmdup.train.cmre.pca.eigenvec

echo "Zipping glm results"
gzip -f train/gwas_results/fixed_effects/ge/*.${rep}.*.glm.linear
