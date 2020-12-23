#!/bin/bash


source /home/aazaidi/anaconda3/etc/profile.d/conda.sh
conda activate fastx

echo "merging vcf files"
m=${1}
w=${2}

bcftools concat genotypes/*.ids.vcf.gz \
-o genotypes/genos_ukb_l1e7_ss500_m${m}_${w}_chr1_20.ids.vcf.gz \
-O z

echo "converting to pgen format"
plink2 \
--double-id \
--make-pgen \
--out genotypes/genos_ukb_l1e7_ss500_m${m}_${w}_chr1_20 \
--vcf genotypes/genos_ukb_l1e7_ss500_m${m}_${w}_chr1_20.ids.vcf.gz

echo "splitting vcf into test and training set"
plink2 \
--keep popfiles/ukb_ss500_d35.${w}.train_iid.txt \
--make-pgen \
--out train/genotypes/genos_ukb_l1e7_ss500_m${m}_${w}_chr1_20.rmdup.train \
--pfile genotypes/genos_ukb_l1e7_ss500_m${m}_${w}_chr1_20 \
--rm-dup exclude-all

plink2 \
--keep popfiles/ukb_ss500_d35.${w}.train_iid.txt \
--make-pgen \
--out test/genotypes/genos_ukb_l1e7_ss500_m${m}_${w}_chr1_20.rmdup.test \
--pfile genotypes/genos_ukb_l1e7_ss500_m${m}_${w}_chr1_20 \
--rm-dup exclude-all
