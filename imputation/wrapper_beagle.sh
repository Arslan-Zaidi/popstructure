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
--export vcf bgz id-paste=iid \
--out train/haplotypes/haps_chr${chrom}.maf0.1


#use test set as reference dataset
plink2 --pfile test/genotypes/genos_gridt100_l1e7_ss750_m0.05_chr1_20.rmdup.test \
--mac 1 \
--chr ${chrom} \
--export vcf bgz id-paste=iid \
--out test/haplotypes/haps_chr${chrom}.mac1


java -jar -Xmx16g ~/mathilab/aazaidi/bin/beagle.18May20.d20.jar \
ref=test/haplotypes/haps_chr${chrom}.mac1.vcf.gz \
gt=train/haplotypes/haps_chr${chrom}.maf0.1.vcf.gz \
chrom=${chrom} \
out=train/imputed_genotypes/imputed_chr${chrom}.beagle

conda activate fastx
bcftools index train/imputed_genotypes/imputed_chr${chrom}.beagle.vcf.gz

#output imputation accuracy scores
zcat train/imputed_genotypes/imputed_chr${chrom}.beagle.vcf.gz | \
grep -v '^#' | cut -f8 > imputed_chr${chrom}.beagle.dr2
