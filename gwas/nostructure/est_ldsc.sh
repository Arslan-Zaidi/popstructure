#!/bin/bash

conda activate ldsc

genotype=${1}
chrom=${2}

mkdir -p train/genotypes/chrs
mkdir -p train/ldsc

plink2 --pfile \
train/genotypes/genos_nostr_s18k_l1e7_chr1_20.rmdup.train \
--mac 1 \
--make-bed \
--out train/genotypes/chrs/chr${chrom}

ldsc.py --bfile train/genotypes/chrs/chr${chrom} \
--l2 \
--ld-wind-kb 1000 \
--out train/ldsc/chr${chrom}

rm train/genotypes/chrs/chr${chrom}.*
