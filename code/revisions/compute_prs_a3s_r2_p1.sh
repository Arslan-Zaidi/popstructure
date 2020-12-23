#!/bin/bash

#script used to calculate PRS given genotype file and variant file with effect size estimates
#specifically for when:
# a3s: variants are Ascertained in Siblings from the 3rd est_effect
# r2: effects Re-estimated in the second (test) set
# p1: polygenic Predictions made in the first (training) set

#paths to R scripts are relative to the root of the project
#(i.e.) $HOME/gwas_bias2

#rep=${1}
corr=${1}

source /home/aazaidi/anaconda3/etc/profile.d/conda.sh
conda activate rpkgs

mkdir -p betas/a3s_r2
mkdir -p prs/a3s_r2_p1

for rep in {1..20}; do \
re_estimate_effects.R \
gwas/grid/genotypes/tau100/ss500/sibs/betas/est_effects.psmooth.e${rep}.all.betas \
gwas/grid/genotypes/tau100/ss500/test/gwas_results/fixed_effects/ge/gwas_gridt100_test.ge.${rep}.${corr}.smooth.glm.linear.gz \
gwas/grid/genotypes/tau100/ss500/revisions/prs_prediction/betas/a3s_r2/effects.smooth.a3s_r2.${corr}.${rep}.betas \
plink; \
done

for rep in {1..20}; do \
plink2 --pfile \
../../train/genotypes/genos_gridt100_l1e7_ss750_m0.05_chr1_20.rmdup.train \
--score betas/a3s_r2/effects.smooth.a3s_r2.${corr}.${rep}.betas cols=scoresums \
--out prs/a3s_r2_p1/prs.smooth.a3s_r2_p1.${corr}.${rep} \
--score-col-nums 4,5,6;
done

#concatenate files

for rep in {1..20}; do \
tail -n+2 prs/a3s_r2_p1/prs.smooth.a3s_r2_p1.${corr}.${rep}.sscore | awk -v OFS="\t" -v i=${rep} '{print i,$0}'; \
done > prs/a3s_r2_p1.smooth.${corr}.all.sscore

for rep in {1..20}; do \
cat betas/a3s_r2/effects.smooth.a3s_r2.${corr}.${rep}.betas | awk -v OFS="\t" -v i=${rep} '{print i,$0}'; \
done > betas/effects.smooth.a3s_r2.${corr}.all.betas
