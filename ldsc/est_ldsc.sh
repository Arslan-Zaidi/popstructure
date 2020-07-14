#!/bin/bash

source /home/aazaidi/anaconda3/etc/profile.d/conda.sh
conda activate ldsc

geno_prefix=${1}
chrom=${2}
mkdir -p train/genotypes/ldsc
mkdir -p train/genotypes/beds

plink2 --pfile train/genotypes/${geno_prefix} \
--mac 2 \
--chr ${chrom} \
--thin-count 150000 \
--make-bed --out train/genotypes/beds/${geno_prefix}.${chrom}.b

ldsc.py --bfile train/genotypes/beds/${geno_prefix}.${chrom}.b \
--l2 \
--ld-wind-kb 100 \
--out train/genotypes/ldsc/${geno_prefix}.ldsc.chr${chrom}
