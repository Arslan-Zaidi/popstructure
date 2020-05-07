#!/bin/bash

echo "separating rare and common variants and pruning for LD"

m=${1}

plink2 \
--indep-pairwise 100 10 0.1 \
--mac 2 \
--max-mac 4 \
--out train/genotypes/genos_complex_l1e7_ss500_m${m}_chr1_20.rmdup.train.repruned \
--pfile train/genotypes/genos_complex_l1e7_ss500_m${m}_chr1_20.rmdup.train \
--write-snplist

plink2 \
--indep-pairwise 100 10 0.1 \
--maf 0.05 \
--out train/genotypes/genos_complex_l1e7_ss500_m${m}_chr1_20.rmdup.train.cmpruned \
--pfile train/genotypes/genos_complex_l1e7_ss500_m${m}_chr1_20.rmdup.train \
--write-snplist


echo "carrying out common PCA"

plink2 \
--pfile train/genotypes/genos_complex_l1e7_ss500_m${m}_chr1_20.rmdup.train \
--extract train/genotypes/genos_complex_l1e7_ss500_m${m}_chr1_20.rmdup.train.cmpruned.prune.in \
--pca 100 \
--out train/genotypes/genos_complex_l1e7_ss500_m${m}_chr1_20.rmdup.train.cmpruned.pca

echo "carrying out rare PCA"

plink2 \
--pfile	train/genotypes/genos_complex_l1e7_ss500_m${m}_chr1_20.rmdup.train \
--extract train/genotypes/genos_complex_l1e7_ss500_m${m}_chr1_20.rmdup.train.repruned.prune.in \
--pca 100 \
--out train/genotypes/genos_complex_l1e7_ss500_m${m}_chr1_20.rmdup.train.repruned.pca

echo "calculating allele frequency"

plink2 \
--pfile	train/genotypes/genos_complex_l1e7_ss500_m${m}_chr1_20.rmdup.train \
--mac 1 \
--make-pgen --out train/genotypes/genos_complex_l1e7_ss500_m${m}_chr1_20.rmdup.train.snps
