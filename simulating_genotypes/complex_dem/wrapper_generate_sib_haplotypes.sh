#!/bin/bash

source /home/aazaidi/anaconda3/etc/profile.d/conda.sh
conda activate py3_clone

chr=${1}
iteration=${2}

mkdir -p test/sibs/genotypes

# first generate vcf file from plink format
plink2 \
--pfile test/genotypes/genos_complex_l1e7_ss500_m0.12_chr1_20.rmdup.test \
--chr ${chr} \
--export vcf-4.2 bgz id-paste=iid id-delim="," \
--out test/genotypes/genos_complex_l1e7_ss500_m0.12_chr${chr}.rmdup.test

#carry out mating using one of the vcf files - not for each file separately
# python mate4sibs.py \
# -v test/genotypes/genos_complex_l1e7_ss500_m0.12_chr1.rmdup.test.vcf.gz \
# -i ${iteration}

#construct sibling haplotypes based on the matings
python make_sib_haplotypes.py \
-v test/genotypes/genos_complex_l1e7_ss500_m0.12_chr${chr}.rmdup.test.vcf.gz \
-m test/sibs/genotypes/geno_complex_assort_${iteration}.sample \
-o test/sibs/genotypes/genos_complex_l1e7_ss500_m0.12_chr${chr}_i${iteration}_sibs_assort

python make_sib_haplotypes.py \
-v test/genotypes/genos_complex_l1e7_ss500_m0.12_chr${chr}.rmdup.test.vcf.gz \
-m test/sibs/genotypes/geno_complex_random_${iteration}.sample \
-o test/sibs/genotypes/genos_complex_l1e7_ss500_m0.12_chr${chr}_i${iteration}_sibs_random

conda activate fastx

bcftools convert \
--haplegendsample2vcf \
test/sibs/genotypes/genos_complex_l1e7_ss500_m0.12_chr${chr}_i${iteration}_sibs_assort.hap,test/sibs/genotypes/genos_complex_l1e7_ss500_m0.12_chr${chr}_i${iteration}_sibs_assort.legend,test/sibs/genotypes/geno_complex_assort_${iteration}.sample \
-Oz -o test/sibs/genotypes/genos_complex_l1e7_ss500_m0.12_chr${chr}_i${iteration}_sibs_assort.vcf.gz

bcftools convert \
--haplegendsample2vcf \
test/sibs/genotypes/genos_complex_l1e7_ss500_m0.12_chr${chr}_i${iteration}_sibs_random.hap,test/sibs/genotypes/genos_complex_l1e7_ss500_m0.12_chr${chr}_i${iteration}_sibs_random.legend,test/sibs/genotypes/geno_complex_random_${iteration}.sample \
-Oz -o test/sibs/genotypes/genos_complex_l1e7_ss500_m0.12_chr${chr}_i${iteration}_sibs_random.vcf.gz

# annoyingly, bcftools is not retaining the IDs provided in the legend file during conversion.
# IDs will have to be manually generated and added to the vcf file
bcftools annotate \
--set-id +'%CHROM\_%POS\_%REF\_%FIRST_ALT' \
genos_complex_l1e7_ss500_m0.12_chr${chr}_i${iteration}_sibs_random.vcf.gz \
-Oz -o tmp_${chr}_${iteration}.gz

mv tmp tmp_${chr}_${iteration}.gz > genos_complex_l1e7_ss500_m0.12_chr${chr}_i${iteration}_sibs_random.vcf.gz

bcftools annotate \
--set-id +'%CHROM\_%POS\_%REF\_%FIRST_ALT' \
genos_complex_l1e7_ss500_m0.12_chr${chr}_i${iteration}_sibs_assort.vcf.gz \
-Oz -o tmp_${chr}_${iteration}.gz

mv tmp_${chr}_${iteration}.gz > genos_complex_l1e7_ss500_m0.12_chr${chr}_i${iteration}_sibs_assort.vcf.gz

# merge the vcf files together into one
bcftools concat test/sibs/genotypes/*_i${iteration}_sibs_random.vcf.gz \
-Oz -o test/sibs/genotypes/genos_complex_l1e7_ss500_m0.12_chr1_20_i${iteration}_sibs_random.vcf.gz

bcftools concat test/sibs/genotypes/*_i${iteration}_sibs_assort.vcf.gz \
-Oz -o test/sibs/genotypes/genos_complex_l1e7_ss500_m0.12_chr1_20_i${iteration}_sibs_assort.vcf.gz

#convert to plink format - needed for GWAS
plink2 \
--vcf genos_complex_l1e7_ss500_m0.12_chr1_20_i${iteration}_sibs_assort.vcf.gz \
--make-pgen --out genos_complex_l1e7_ss500_m0.12_chr1_20_i${iteration}_sibs_assort

#remove intermediate files that are taking space
