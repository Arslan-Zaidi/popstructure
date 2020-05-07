#!/bin/bash

rep=${1}

source /home/aazaidi/anaconda3/etc/profile.d/conda.sh
conda activate rpkgs

mkdir -p phenotypes/noge

echo "Simulating nonheritable phenotypes"
#simulate nongenetic phenotype
Rscript scripts/simphenotype/simphenotype_noge.R \
genos_complex_l1e7_ss500_m0.07.train.pop \
2 \
phenotypes/noge/pheno_complex_noge_s9k.train.${rep}.txt \
${rep}

echo "Running GWAS on non-heritable phenotypes"
#carry out GWAS with different modes of correction

#carry out GWAS for nongenetic phenotypes

#1. no correction
bash scripts/gwas/gwas.sh \
train/genotypes/genos_complex_l1e7_ss500_m0.07_chr1_20.rmdup.train \
phenotypes/noge/pheno_complex_noge_s9k.train.${rep}.txt \
train/gwas_results/fixed_effects/noge/gwas_complex_train.noge.${rep}.pcs0

#2. common pca
bash scripts/gwas/gwas.sh \
train/genotypes/genos_complex_l1e7_ss500_m0.07_chr1_20.rmdup.train \
phenotypes/noge/pheno_complex_noge_s9k.train.${rep}.txt \
train/gwas_results/fixed_effects/noge/gwas_complex_train.noge.${rep}.cm \
train/genotypes/genos_complex_l1e7_ss500_m0.07_chr1_20.rmdup.train.cmpruned.pca.eigenvec

#3. rare pca
bash scripts/gwas/gwas.sh \
train/genotypes/genos_complex_l1e7_ss500_m0.07_chr1_20.rmdup.train \
phenotypes/noge/pheno_complex_noge_s9k.train.${rep}.txt \
train/gwas_results/fixed_effects/noge/gwas_complex_train.noge.${rep}.re \
train/genotypes/genos_complex_l1e7_ss500_m0.07_chr1_20.rmdup.train.repruned.pca.eigenvec

echo "Zipping glm results"
gzip train/gwas_results/fixed_effects/noge/*.${rep}.*.glm.linear
