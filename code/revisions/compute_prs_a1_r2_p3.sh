#!/bin/bash

#script used to calculate PRS given genotype file and variant file with effect size estimates

#script used to calculate PRS given genotype file and variant file with effect size estimates
#specifically for when:
# a1: variants are Ascertained in training set
# r2: effects Re-estimated in the test set
# p3: polygenic Predictions made in the third set (used to generate sib haplotypes)

#paths to R scripts are relative to the root of the project
#(i.e.) $HOME/gwas_bias2

#rep=${1}
corr=${1}

source /home/aazaidi/anaconda3/etc/profile.d/conda.sh
conda activate rpkgs

mkdir -p betas/a1_r2
mkdir -p prs/a1_r2_p3

for rep in {1..20}; do \
re_estimate_effects.R \
gwas/grid/genotypes/tau100/ss500/train/betas/est_effects.${rep}.smooth.nc.betas \
gwas/grid/genotypes/tau100/ss500/test/gwas_results/fixed_effects/ge/gwas_gridt100_test.ge.${rep}.${corr}.smooth.glm.linear.gz \
gwas/grid/genotypes/tau100/ss500/revisions/prs_prediction/betas/a1_r2/effects.smooth.a1_r2.${corr}.${rep}.betas \
plink; \
done


for rep in {1..20}; do \
plink2 --pfile \
../../sibs/genotypes/genos_gridt100_l1e7_ss750_m0.05_chr1_20.rmdup.sib \
--score betas/a1_r2/effects.smooth.a1_r2.${corr}.${rep}.betas cols=scoresums \
--out prs/a1_r2_p3/prs.smooth.a1_r2_p3.${corr}.${rep} \
--score-col-nums 4,5,6,7;
done

#concatenate files

for rep in {1..20}; do \
tail -n+2 prs/a1_r2_p3/prs.smooth.a1_r2_p3.${corr}.${rep}.sscore | awk -v OFS="\t" -v i=${rep} '{print i,$0}'; \
done > prs/a1_r2_p3.smooth.${corr}.all.sscore

for rep in {1..20}; do \
tail -n+2 betas/a1_r2/effects.smooth.a1_r2.${corr}.${rep}.betas | awk -v OFS="\t" -v i=${rep} '{print i,$0}'; \
done > betas/effects.smooth.a1_r2.${corr}.all.betas
