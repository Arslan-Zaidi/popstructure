#!/bin/bash

source /home/aazaidi/anaconda3/etc/profile.d/conda.sh
conda activate py3_clone

chr=${1}
effects_rep=${2}
npairs=${3}

#working directory
#gwas_bias2/gwas/grid/genotypes/tau100/ss500/sibs/

mkdir -p gwas_results/

echo "GWAS for assortative mating"

python sib_gwas.py \
--phenotypes phenotypes/genos_gridt100_l1e7_ss750_m0.05_chr20_e${effects_rep}_sibs_assort.pheno \
--haplotypes genotypes/genos_gridt100_l1e7_ss750_m0.05_chr${chr}_sibs_assort.haps.npz \
--legend genotypes/genos_gridt100_l1e7_ss750_m0.05_chr${chr}_sibs_assort.legend \
--outpre gwas_results/gwas_gridt100_l1e7_ss750_m0.05_chr${chr}_e${effects_rep}_n${npairs}_sibs_assort \
--sample mates/genos_gridt100_assort_1_4.sample \
--npairs ${npairs}

gzip gwas_results/gwas_gridt100_l1e7_ss750_m0.05_chr${chr}_e${effects_rep}_n${npairs}_sibs_assort.*.linear
