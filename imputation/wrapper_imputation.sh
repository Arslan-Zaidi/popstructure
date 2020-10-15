#!/bin/bash

#script used to impute variants from a reference panel

chrom=${1}
mkdir -p train/haplotypes
mkdir -p test/haplotypes
mkdir -p train/imputed_genotypes

#use training set as data to be imputed
#exclude rare variants to be imputed
plink2 --pfile train/genotypes/genos_gridt100_l1e7_ss750_m0.05_chr1_20.rmdup.train.snps \
--chr ${chrom} \
--maf 0.001 \
--export haps bgz \
--out train/haplotypes/haps_chr${chrom}.maf0.1


#use test set as reference dataset
plink2 --pfile test/genotypes/genos_gridt100_l1e7_ss750_m0.05_chr1_20.rmdup.test \
--mac 1 \
--chr ${chrom} \
--export hapslegend bgz \
--out test/haplotypes/haps_chr${chrom}.mac1

#generate recombination map file
gen_map.R \
gwas/grid/genotypes/tau100/ss500/test/haplotypes/haps_chr${chrom}.mac1.legend \
gwas/grid/genotypes/tau100/ss500/test/haplotypes/haps_chr${chrom}.mac1.map



impute2 \
-use_prephased_g \
-m test/haplotypes/haps_chr${chrom}.mac1.map \
-h test/haplotypes/haps_chr${chrom}.mac1.haps.gz \
-l test/haplotypes/haps_chr${chrom}.mac1.legend \
-known_haps_g train/haplotypes/haps_chr${chrom}.maf0.1.haps.gz \
-int 1 100 \
-allow_large_regions \
-o train/imputed_genotypes/imputed_chr${chrom}.impute2
