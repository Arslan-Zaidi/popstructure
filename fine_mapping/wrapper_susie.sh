#!/bin/bash

mkdir -p train/susie_output/

rep=${1}
correction=${2}

#paths are relative to the rproject root (i.e ~/gwas_bias2)
finemap.R \
gwas/grid/genotypes/tau100/ss500/train/gwas_results/fixed_effects/ge/gwas_gridt100_train.ge.${rep}.${correction}.smooth.glm.linear.gz \
5e-04  \
gwas/grid/genotypes/tau100/ss500/train/ldmats/ \
gwas/grid/genotypes/tau100/ss500/train/susie_output/susie_${rep}.${correction}.betas
