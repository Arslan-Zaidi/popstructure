#!/bin/bash

chrom=${1}

source /home/aazaidi/anaconda3/etc/profile.d/conda.sh

# conda activate py3_clone
#
# python generate_genos_nostructure.py \
# -s 18000 \
# -c ${chromosome} \
# -N 1035000 \
# -L 10000000 \
# -r 1e-08 \
# -u 1e-08 \
# -o genotypes/genos_nostr_s18k_l1e7

conda activate fastx

head -n6 genotypes/genos_nostr_s18k_l1e7_chr${chrom}.vcf > genotypes/header_${chrom}.txt

cat genotypes/genos_nostr_s18k_l1e7_chr${chrom}.vcf | awk -v OFS="\t" 'NR>6 {$3=$1"_"$2"_A_T"; $4="A"; $5="T"; print ;}'> genotypes/tmp_${chrom}

cat genotypes/header_${chrom}.txt genotypes/tmp_${chrom} | bgzip > genotypes/genos_nostr_s18k_l1e7_chr${chrom}.ids.vcf.gz

rm tmp_${chrom}
rm header_${chrom}.txt
