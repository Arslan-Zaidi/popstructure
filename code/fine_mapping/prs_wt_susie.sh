#!/bin/bash


corr=${1}

for rep in {1..20}; do \
plink2 \
--pfile test/genotypes/genos_gridt100_l1e7_ss750_m0.05_chr1_20.rmdup.test \
--score train/susie_output/susie_${rep}.${corr}.betas cols=scoresums \
--score-col-nums 3,6 \
--out test/susie_prs/prs_susie_smooth.${rep}.${corr}; done

for rep in {1..20}; do \
tail -n+2 test/susie_prs/prs_susie_smooth.${rep}.${corr}.sscore | \
awk -v OFS="\t" -v i=$rep '{print i, $0}'; done >test/susie_prs/prs_susie_smooth.all.${corr}.sscore


for rep in {1..20}; do \
tail -n+2 train/susie_output/susie_${rep}.${corr}.betas | \
awk -v OFS="\t" -v i=$rep '{print i, $0}'; done >train/susie_output/susie_all.${corr}.betas
