#!/bin/bash

#script used to calculate PRS given genotype file and variant file with effect size estimates
#specifically for when:
# a3s: variants are Ascertained in Siblings from the 3rd est_effect
# r1: effects Re-estimated in the first (training) set
# p2: polygenic Predictions made in the second (test) set

#paths to R scripts are relative to the root of the project
#(i.e.) $HOME/gwas_bias2

#rep=${1}
corr=${1}

source /home/aazaidi/anaconda3/etc/profile.d/conda.sh
conda activate rpkgs

mkdir -p betas/a3s_r1
mkdir -p prs/a3s_r1_p2

for rep in {1..20}; do \
re_estimate_effects.R \
gwas/grid/genotypes/tau100/ss500/sibs/betas/est_effects.psmooth.e${rep}.all.betas \
gwas/grid/genotypes/tau100/ss500/train/gwas_results/fixed_effects/ge/gwas_gridt100_train.ge.${rep}.${corr}.smooth.glm.linear.gz \
gwas/grid/genotypes/tau100/ss500/revisions/prs_prediction/betas/a3s_r1/effects.smooth.a3s_r1.${corr}.${rep}.betas; \
done

for rep in {1..20}; do \
plink2 --pfile \
../../test/genotypes/genos_gridt100_l1e7_ss750_m0.05_chr1_20.rmdup.test \
--score betas/a3s_r1/effects.smooth.a3s_r1.${corr}.${rep}.betas cols=scoresums \
--out prs/a3s_r1_p2/prs.smooth.a3s_r1_p2.${corr}.${rep} \
--score-col-nums 5;
done

#concatenate files

for rep in {1..20}; do \
tail -n+2 prs/a3s_r1_p2/prs.smooth.a3s_r1_p2.${corr}.${rep}.sscore | awk -v OFS="\t" -v i=${rep} '{print i,$0}'; \
done > prs/a3s_r1_p2.smooth.${corr}.all.sscore
