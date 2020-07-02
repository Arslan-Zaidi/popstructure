#!/bin/bash

source /home/aazaidi/anaconda3/etc/profile.d/conda.sh
conda activate py3_clone

chr=${1}
m=${2}

python complex_dem_2.py \
 --sample_size 500 \
 --chr ${chr} \
 --outpre genotypes/genos_complex_l1e7_ss500_m${m} \
 --migrate ${m}

#gzip genotypes/genos_complex_l1e7_ss500_m0.20_chr${chr}.vcf

conda activate fastx

#step 1. save header to a separate file
cat genotypes/genos_complex_l1e7_ss500_m${m}_chr${chr}.vcf | head -n6 > genotypes/header_${chr}.txt

#step 2: use awk to replace the ref/alt alleles with A/T, then concatenate with header and then bgzip
#also, add ID column to each variant
cat genotypes/header_${chr}.txt <(cat genotypes/genos_complex_l1e7_ss500_m${m}_chr${chr}.vcf | awk -v OFS="\t" 'NR>6 {$3=$1"_"$2"_A_T";
$4="A"; $5="T"; print ;}') | bgzip > genotypes/genos_complex_l1e7_ss500_m${m}_chr${chr}.ids.vcf.gz

bcftools index genotypes/genos_complex_l1e7_ss500_m${m}_chr${chr}.ids.vcf.gz
