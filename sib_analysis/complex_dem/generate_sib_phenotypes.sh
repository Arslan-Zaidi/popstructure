#!/bin/bash

source /home/aazaidi/anaconda3/etc/profile.d/conda.sh
conda activate rpkgs

# script generates phenotypes for the sibling GWAS

mkdir -p test/sibs/gvalues
mkdir -p test/sibs/phenotypes

effects_rep=${1} # iteration of effect sizes

###### generate phenotypes for the sibs

echo "concaetnating genetic values"
### 1. concatenate genetic values across chromosomes
for j in {1..20}; do \
tail -n+2 -q test/sibs/gvalues/genos_complex_l1e7_ss500_m0.08_chr${j}_e${effects_rep}_sibs_assort.sscore | \
awk -v i=${effects_rep} '{print i,$0}'; done > test/sibs/gvalues/genos_complex_l1e7_ss500_m0.08_chr1_20_e${effects_rep}_sibs_assort.sscore

# tail -n+2 -q test/sibs/genotypes/genos_complex_l1e7_ss500_m0.08_chr*_e${effects_rep}_sibs_random.sscore | \
# awk -v i=${effects_rep} '{print i,$0}' > test/sibs/gvalues/genos_complex_l1e7_ss500_m0.08_chr1_20_e${effects_rep}_sibs_random.sscore

echo "generating phenotypes"
### 2. Add environmental effects to these valuese

Rscript simphenotype_sibs_ge.R \
test/sibs/gvalues/genos_complex_l1e7_ss500_m0.08_chr1_20_e${effects_rep}_sibs_assort.sscore \
test/sibs/genotypes/genos_complex_assort_1_4.sample \
test/sibs/phenotypes/genos_complex_l1e7_ss500_m0.08_chr1_20_e${effects_rep}_sibs_assort.ge.pheno \
${effects_rep}

# Rscript simphenotype_sibs_ge.R \
# test/sibs/gvalues/genos_complex_l1e7_ss500_m0.08_chr1_20_e${effects_rep}_sibs_random.sscore \
# test/sibs/genotypes/genos_complex_random_1_4.sample \
# test/sibs/phenotypes/genos_complex_l1e7_ss500_m0.08_chr1_20_e${effects_rep}_sibs_random.ge.pheno \
# ${effects_rep}
