#!/bin/bash

source /home/aazaidi/anaconda3/etc/profile.d/conda.sh
conda activate rpkgs

# script generates phenotypes for the sibling GWAS

mkdir -p phenotypes

effects_rep=${1} # iteration of effect sizes

###### generate phenotypes for the sibs

echo "concaetnating genetic values"
### 1. concatenate genetic values across chromosomes
for j in {1..20}; do \
tail -n+2 -q gvalues/genos_gridt100_l1e7_ss750_m0.05_chr${j}_e${effects_rep}_sibs_assort.sscore | \
awk -v i=${effects_rep} '{print i,$0}'; done > gvalues/genos_gridt100_l1e7_ss750_m0.05_chr1_20_e${effects_rep}_sibs_assort.sscore

echo "generating phenotypes"
### 2. Add environmental effects to these valuese

Rscript simphenotype_sibs_ge.R \
gvalues/genos_gridt100_l1e7_ss750_m0.05_chr1_20_e${effects_rep}_sibs_assort.sscore \
mates/genos_gridt100_assort_1_4.sample \
phenotypes/genos_gridt100_l1e7_ss750_m0.05_chr1_20_e${effects_rep}_sibs_assort.pheno \
${effects_rep}
