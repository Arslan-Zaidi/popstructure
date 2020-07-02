
#!/bin/bash

rep=${1}
source /home/aazaidi/anaconda3/etc/profile.d/conda.sh
#
# plink2 --pfile train/genotypes/genos_grid_d36_m0.07_s500_t9_chr1_20.rmdup.train.snps \
# --read-freq train/genotypes/genos_grid_d36_m0.07_s500_t9_chr1_20.rmdup.train.snps.frq.afreq \
# --glm hide-covar \
# --pheno phenotypes/ge/pheno_gridt9_ge_s9k.train.${rep}.txt \
# --pheno-name smooth \
# --covar genos_grid_d36_m0.07_s500_t9.train.pop \
# --covar-col-nums 4,5 \
# --out train/gwas_results/fixed_effects/ge/gwas_gridt9_train.ge.${rep}.geo

gzip train/gwas_results/fixed_effects/ge/gwas_gridt9_train.ge.${rep}.geo.smooth.glm.linear
