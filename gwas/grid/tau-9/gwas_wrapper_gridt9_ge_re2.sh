
#!/bin/bash

rep=${1}

source /home/aazaidi/anaconda3/etc/profile.d/conda.sh
conda activate rpkgs

mkdir -p phenotypes/ge

#3. rare pca
bash scripts/gwas/gwas.sh \
train/genotypes/genos_grid_d36_m0.07_s500_t9_chr1_20.rmdup.train.snps \
phenotypes/ge/pheno_gridt9_ge_s9k.train.${rep}.txt \
train/gwas_results/fixed_effects/ge/gwas_gridt9_train.ge.${rep}.re2 \
train/genotypes/genos_grid_d36_m0.07_s500_t9_chr1_20.rmdup.train.snps.re2.pca.eigenvec
