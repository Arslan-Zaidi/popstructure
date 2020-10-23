#!/bin/bash
conda activate fastx

input=$1
chrom=$2

mkdir -p train/germline

echo "removing rare variants"

#bcftools view -i 'MAF>0.05' ${input} | gzip -> ${input}.cm
plink2 --pfile ${input} \
--chr ${chrom} \
--maf 0.05 \
--export haps \
--out train/germline/genos_ukb_chr${chrom}.cm

echo "converting hap/sample format to ped/map format"
#convert hap sample format to ped map format
impute_to_ped train/germline/genos_ukb_chr${chrom}.cm.haps \
train/germline/genos_ukb_chr${chrom}.cm.sample \
train/germline/germ_ukb_chr${chrom}.cm

echo "running germline"
#run germline
germline \
-input train/germline/germ_ukb_chr${chrom}.cm.ped \
train/germline/germ_ukb_chr${chrom}.cm.map \
-output germ_ukb_ibd_chr${chrom}.cm
