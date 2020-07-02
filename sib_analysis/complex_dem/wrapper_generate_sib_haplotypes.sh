#!/bin/bash

source /home/aazaidi/anaconda3/etc/profile.d/conda.sh
conda activate py3_clone

chr=${1}
effects_rep=${2}

mkdir -p test/sibs/genotypes

## first generate vcf file from plink format
# plink2 \
# --pfile test/genotypes/genos_complex_l1e7_ss500_m0.08_chr1_20.rmdup.test \
# --chr ${chr} \
# --export vcf-4.2 bgz id-paste=iid id-delim="," \
# --out test/genotypes/genos_complex_l1e7_ss500_m0.08_chr${chr}.rmdup.test

echo "carrying out matings"

## carry out mating using one of the vcf files - not for each file separately
# for i in {1..4}; do
# python mate4sibs.py \
# -v test/genotypes/genos_complex_l1e7_ss500_m0.08_chr1.rmdup.test.vcf.gz \
# -i ${i} \
# -o test/sibs/genotypes/genos_complex; done

echo "generating sib haplotypes assortatively"

#construct sibling haplotypes based on the matings
python make_sib_haplotypes.py \
-v test/genotypes/genos_complex_l1e7_ss500_m0.08_chr${chr}.rmdup.test.vcf.gz \
-m test/sibs/genotypes/genos_complex_assort_1_4.sample \
-e train/genotypes/genos_complex_l1e7_ss500_m0.08_chr1_20.rmdup.train.${effects_rep}.thinned_100kb.effects \
-o test/sibs/genotypes/genos_complex_l1e7_ss500_m0.08_chr${chr}_e${effects_rep}_sibs_assort

echo "generating sib haplotypes randomly"

python make_sib_haplotypes.py \
-v test/genotypes/genos_complex_l1e7_ss500_m0.08_chr${chr}.rmdup.test.vcf.gz \
-m test/sibs/genotypes/genos_complex_random_1_4.sample \
-e train/genotypes/genos_complex_l1e7_ss500_m0.08_chr1_20.rmdup.train.${effects_rep}.thinned_100kb.effects \
-o test/sibs/genotypes/genos_complex_l1e7_ss500_m0.08_chr${chr}_e${effects_rep}_sibs_random
