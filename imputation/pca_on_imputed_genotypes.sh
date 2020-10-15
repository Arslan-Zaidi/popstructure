#!/bin/bash


# bcftools concat train/imputed_genotypes/imputed_chr*.beagle.vcf.gz \
# -o train/imputed_genotypes/imputed_chr1_20.beagle.vcf.gz \
# -O z
#
#
# plink2 \
# --double-id \
# --make-pgen \
# --out train/imputed_genotypes/imputed_chr1_20.plink \
# --vcf train/imputed_genotypes/imputed_chr1_20.beagle.vcf.gz
#
# cp train/genotypes/genos_gridt100_l1e7_ss750_m0.05_chr1_20.rmdup.train.re.all.snplist \
# train/imputed_genotypes/
#
# cp train/genotypes/genos_gridt100_l1e7_ss750_m0.05_chr1_20.rmdup.train.cm.200k.snplist \
# train/imputed_genotypes/

echo "carrying out pca on common variants"
plink2 \
--pfile train/imputed_genotypes/imputed_chr1_20.plink \
--extract train/imputed_genotypes/genos_gridt100_l1e7_ss750_m0.05_chr1_20.rmdup.train.cm.200k.snplist \
--pca 100 \
--out train/imputed_genotypes/imputed_chr1_20.cm.200k.pca

echo "carrying out pca on rare variants"
plink2 \
--pfile train/imputed_genotypes/imputed_chr1_20.plink \
--extract train/imputed_genotypes/genos_gridt100_l1e7_ss750_m0.05_chr1_20.rmdup.train.re.all.snplist \
--pca 100 \
--out train/imputed_genotypes/imputed_chr1_20.re.all.pca
