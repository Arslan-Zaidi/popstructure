#!/bin/bash

chrom=${1} #chromosome number

#step 1. save header to a separate file
zcat genotypes/genos_complex_l1e7_ss500_m0.07_chr${chrom}.vcf.gz | head -n6 > genotypes/header_${chrom}.txt

#step 2: use awk to replace the ref/alt alleles with A/T, then concatenate with header and then bgzip
#also, add ID column to each variant
cat genotypes/header_${chrom}.txt <(zcat genotypes/genos_complex_l1e7_ss500_m0.07_chr${chrom}.vcf.gz | awk -v OFS="\t" 'NR>6 {$3=$1"_"$2"_A_T"; 
$4="A"; $5="T"; print ;}') | bgzip > genotypes/genos_complex_l1e7_ss500_m0.07_chr${chrom}.ids.vcf.gz

bcftools index genotypes/genos_complex_l1e7_ss500_m0.07_chr${chrom}.ids.vcf.gz
