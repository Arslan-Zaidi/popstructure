#!/bin/bash

echo "separating rare and common variants and variant thinning"

plink2 \
--mac 2 \
--max-mac 4 \
--out train/genotypes/genos_grid_d36_m0.07_s500_t9_chr1_20.rmdup.train.snps.re \
--pfile train/genotypes/genos_grid_d36_m0.07_s500_t9_chr1_20.rmdup.train \
--write-snplist

plink2 \
--maf 0.05 \
--out train/genotypes/genos_grid_d36_m0.07_s500_t9_chr1_20.rmdup.train.snps.cm \
--pfile train/genotypes/genos_grid_d36_m0.07_s500_t9_chr1_20.rmdup.train \
--write-snplist

echo "carrying out common PCA"

plink2 \
--pfile train/genotypes/genos_grid_d36_m0.07_s500_t9_chr1_20.rmdup.train \
--extract train/genotypes/genos_grid_d36_m0.07_s500_t9_chr1_20.rmdup.train.snps.cm.snplist \
--thin-count 200000 \
--pca 100 \
--out train/genotypes/genos_grid_d36_m0.07_s500_t9_chr1_20.rmdup.train.snps.re.200k.pca

echo "carrying out rare PCA"

plink2 \
--pfile	train/genotypes/genos_grid_d36_m0.07_s500_t9_chr1_20.rmdup.train \
--extract train/genotypes/genos_grid_d36_m0.07_s500_t9_chr1_20.rmdup.train.snps.re.snplist \
--thin-count 1000000 \
--pca 100 \
--out train/genotypes/genos_grid_d36_m0.07_s500_t9_chr1_20.rmdup.train.snps.re.1M.pca
