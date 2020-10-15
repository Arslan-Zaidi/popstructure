#!/bin/bash

mkdir -p train/ldmats

chrom=${1}
start=${2}
stop=${3}

plink2 \
--pfile train/genotypes/genos_gridt100_l1e7_ss750_m0.05_chr1_20.rmdup.train.snps \
--chr ${chrom} --from-bp ${start} --to-bp ${stop} \
--make-bed --out train/ldmats/genos_chr${chrom}_${start}_${stop}

plink \
--bfile train/ldmats/genos_chr${chrom}_${start}_${stop} \
--r square gz \
--out train/ldmats/ld_chr${chrom}_${start}_${stop}

rm train/ldmats/genos_chr${chrom}_${start}_${stop}.*
