#!/bin/bash


#script used to calculate PRS given genotype file and variant file with effect size estimates
#specifically for when:
# a1: variants are Ascertained in training set
# r3s: effects Re-estimated in siblings
# p2: polygenic Predictions made in the second set (test set)

#paths to R scripts are relative to the root of the project
#(i.e.) $HOME/gwas_bias2

#corr=${1}

source /home/aazaidi/anaconda3/etc/profile.d/conda.sh
conda activate rpkgs

mkdir -p betas/a1_r3s
mkdir -p prs/a1_r3s_p2

for rep in {1..20}; do \
re_estimate_effects.R \
gwas/grid/genotypes/tau100/ss500/train/betas/est_effects.${rep}.smooth.nc.betas \
gwas/grid/genotypes/tau100/ss500/sibs/gwas_results/e${rep}/gwas_gridt100_l1e7_ss750_m0.05_chr1_20_e${rep}_n9000_sibs_assort.smooth.glm.linear.gz \
gwas/grid/genotypes/tau100/ss500/revisions/prs_prediction/betas/a1_r3s/effects.smooth.a1_r3s.ncorr.${rep}.betas \
sibling; \
done

for rep in {1..20}; do \
plink2 --pfile \
../../test/genotypes/genos_gridt100_l1e7_ss750_m0.05_chr1_20.rmdup.test \
--score betas/a1_r3s/effects.smooth.a1_r3s.ncorr.${rep}.betas cols=scoresums \
--out prs/a1_r3s_p2/prs.smooth.a1_r3s_p2.ncorr.${rep} \
--score-col-nums 4,5,6,7;
done

#concatenate files

for rep in {1..20}; do \
tail -n+2 prs/a1_r3s_p2/prs.smooth.a1_r3s_p2.ncorr.${rep}.sscore | awk -v OFS="\t" -v i=${rep} '{print i,$0}'; \
done > prs/a1_r3s_p2.smooth.ncorr.all.sscore

for rep in {1..20}; do \
tail -n+2 betas/a1_r3s/effects.smooth.a1_r3s.ncorr.${rep}.betas | awk -v OFS="\t" -v i=${rep} '{print i,$0}'; \
done > betas/effects.smooth.a1_r3s.ncorr.all.betas
