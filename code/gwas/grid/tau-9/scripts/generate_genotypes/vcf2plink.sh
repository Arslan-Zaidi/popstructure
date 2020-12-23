#!/bin/bash

echo "merging vcf files"

bcftools concat genotypes/*.ids.vcf.gz \
-o genotypes/genos_nostr_s18k_l1e7_chr1_20.vcf.gz \
-O z


echo "converting to pgen format"
plink2 \
--double-id \
--make-pgen \
--out genotypes/genos_nostr_s18k_l1e7_chr1_20 \
--vcf genotypes/genos_nostr_s18k_l1e7_chr1_20.vcf.gz

echo "splitting vcf into test and training set"
plink2 \
--keep iid1_test.txt \
--make-pgen \
--out train/genotypes/genos_nostr_s18k_l1e7_chr1_20.rmdup.train \
--pfile genotypes/genos_nostr_s18k_l1e7_chr1_20 \
--rm-dup exclude-all

plink2 \
--keep iid2_test.txt \
--make-pgen \
--out test/genotypes/genos_nostr_s18k_l1e7_chr1_20.rmdup.test \
--pfile genotypes/genos_nostr_s18k_l1e7_chr1_20 \
--rm-dup exclude-all

