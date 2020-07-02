#!/bin/bash


source /home/aazaidi/anaconda3/etc/profile.d/conda.sh
conda activate fastx

echo "merging vcf files"
m=${1}
bcftools concat genotypes/genos_gridt100_l1e7_ss750_m0.05_chr*.ids.vcf.gz \
-o genotypes/genos_gridt100_l1e7_ss750_m${m}_chr1_20.ids.vcf.gz \
-O z

echo "converting to pgen format"
plink2 \
--double-id \
--make-pgen \
--out genotypes/genos_gridt100_l1e7_ss750_m${m}_chr1_20 \
--vcf genotypes/genos_gridt100_l1e7_ss750_m${m}_chr1_20.ids.vcf.gz

echo "splitting vcf into test and training set"
plink2 \
--keep iid_train.txt \
--make-pgen \
--out train/genotypes/genos_gridt100_l1e7_ss750_m${m}_chr1_20.rmdup.train \
--pfile genotypes/genos_gridt100_l1e7_ss750_m${m}_chr1_20 \
--rm-dup exclude-all

plink2 \
--keep iid_test.txt \
--make-pgen \
--out test/genotypes/genos_gridt100_l1e7_ss750_m${m}_chr1_20.rmdup.test \
--pfile genotypes/genos_gridt100_l1e7_ss750_m${m}_chr1_20 \
--rm-dup exclude-all

echo "splitting vcf into test and training set"
plink2 \
--keep iid_sib.txt \
--make-pgen \
--out sibs/genotypes/genos_gridt100_l1e7_ss750_m${m}_chr1_20.rmdup.sib \
--pfile genotypes/genos_gridt100_l1e7_ss750_m${m}_chr1_20 \
--rm-dup exclude-all
