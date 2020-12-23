#!/bin/bash

echo "separating rare and common variants and variant thinning"

m=0.05


plink2 \
--mac 2 \
--max-mac 4 \
--out train/genotypes/genos_gridt100_l1e7_ss750_m${m}_chr1_20.rmdup.train.re.all \
--pfile train/genotypes/genos_gridt100_l1e7_ss750_m${m}_chr1_20.rmdup.train \
--write-snplist

plink2 \
--maf 0.05 \
--out train/genotypes/genos_gridt100_l1e7_ss750_m${m}_chr1_20.rmdup.train.cm \
--pfile train/genotypes/genos_gridt100_l1e7_ss750_m${m}_chr1_20.rmdup.train \
--write-snplist

echo "carrying out common PCA"

plink2 \
--pfile train/genotypes/genos_gridt100_l1e7_ss750_m${m}_chr1_20.rmdup.train \
--extract train/genotypes/genos_gridt100_l1e7_ss750_m${m}_chr1_20.rmdup.train.cm.snplist \
--thin-count 200000 \
--pca 100 \
--out train/genotypes/genos_gridt100_l1e7_ss750_m${m}_chr1_20.rmdup.train.cm.200k

echo "carrying out rare PCA"

plink2 \
--pfile	train/genotypes/genos_gridt100_l1e7_ss750_m${m}_chr1_20.rmdup.train \
--extract train/genotypes/genos_gridt100_l1e7_ss750_m${m}_chr1_20.rmdup.train.re.all.snplist \
--thin-count 1000000 \
--pca 100 \
--out train/genotypes/genos_gridt100_l1e7_ss750_m${m}_chr1_20.rmdup.train.re.1M

echo "calculating allele frequency"

plink2 \
--pfile	train/genotypes/genos_gridt100_l1e7_ss750_m${m}_chr1_20.rmdup.train \
--mac 1 \
--write-snplist \
--out train/genotypes/genos_gridt100_l1e7_ss750_m${m}_chr1_20.rmdup.train.snps

plink2 \
--pfile train/genotypes/genos_gridt100_l1e7_ss750_m${m}_chr1_20.rmdup.train \
--extract train/genotypes/genos_gridt100_l1e7_ss750_m${m}_chr1_20.rmdup.train.snps \
--freq \
--out train/genotypes/genos_gridt100_l1e7_ss750_m${m}_chr1_20.rmdup.train.snps.frq
