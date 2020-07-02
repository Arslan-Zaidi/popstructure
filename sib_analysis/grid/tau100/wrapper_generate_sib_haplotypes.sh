#!/bin/bash

source /home/aazaidi/anaconda3/etc/profile.d/conda.sh
conda activate py3_clone

chr=${1}
#effects_rep=${2}

mkdir -p mates

# first generate vcf file from plink format
# plink2 \
# --pfile genotypes/genos_gridt100_l1e7_ss750_m0.05_chr1_20.rmdup.sib \
# --chr ${chr} \
# --export vcf-4.2 bgz id-paste=iid id-delim="," \
# --out genotypes/genos_gridt100_l1e7_ss750_m0.05_chr${chr}.rmdup.sib

echo "carrying out matings"

## carry out mating using one of the vcf files - not for each file separately
# for i in {1..4}; do
# python mate4sibs.py \
# -v genotypes/genos_gridt100_l1e7_ss750_m0.05_chr1.rmdup.sib.vcf.gz \
# -i ${i} \
# -o mates/genos_gridt100; done

echo "generating sib haplotypes assortatively"

#construct sibling haplotypes based on the matings
python make_sib_haplotypes.py \
-v genotypes/genos_gridt100_l1e7_ss750_m0.05_chr${chr}.rmdup.sib.vcf.gz \
-m mates/genos_gridt100_assort_1_4.sample \
-o genotypes/genos_gridt100_l1e7_ss750_m0.05_chr${chr}_sibs_assort
