#!/bin/bash

source /home/aazaidi/anaconda3/etc/profile.d/conda.sh
conda activate py3_clone

chr=${1}
effects_rep=${2}

mkdir -p test/sibs/gwas_results/

# echo "GWAS for random mating"
#
# python sib_gwas.py \
# --phenotypes test/sibs/phenotypes/genos_complex_l1e7_ss500_m0.08_chr1_20_e${effects_rep}_sibs_random.ge.pheno \
# --haplotypes test/sibs/genotypes/genos_complex_l1e7_ss500_m0.08_chr${chr}_e${effects_rep}_sibs_random.haps.npz \
# --legend test/sibs/genotypes/genos_complex_l1e7_ss500_m0.08_chr${chr}_e${effects_rep}_sibs_random.legend \
# --outpre test/sibs/gwas_results/gwas_complex_l1e7_ss500_m0.08_chr${chr}_e${effects_rep}_sibs_random

echo "GWAS for assortative mating"

python sib_gwas.py \
--phenotypes test/sibs/phenotypes/genos_complex_l1e7_ss500_m0.08_chr1_20_e${effects_rep}_sibs_assort.ge.pheno \
--haplotypes test/sibs/genotypes/genos_complex_l1e7_ss500_m0.08_chr${chr}_e1_sibs_assort.haps.npz \
--legend test/sibs/genotypes/genos_complex_l1e7_ss500_m0.08_chr${chr}_e1_sibs_assort.legend \
--outpre test/sibs/gwas_results/gwas_complex_l1e7_ss500_m0.08_chr${chr}_e${effects_rep}_sibs_assort


gzip test/sibs/gwas_results/gwas_complex_l1e7_ss500_m0.08_chr${chr}_e${effects_rep}_sibs_*.linear
