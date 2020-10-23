#!/bin/bash

for i in {5e3,1e4,5e4,1e5}; \
do tail -n+2 "gwas/ukb/train/genotypes/revisions/genos_ukb_l1e7_ss500_m0.08_uniform_chr1_20.rmdup.train.cm.${i}.pca.eigenvec" | \
awk -v OFS="\t" -v nvar=${i} '{print nvar,$0}'; done \
> gwas/ukb/train/genotypes/revisions/genos.ukb_cm.nvars.eigenvec

for i in {5e3,1e4,5e4,1e5}; \
do tail -n+2 "gwas/grid/genotypes/tau-9/ss500/train/genotypes/revisions/genos_grid_d36_m0.07_s500_t9_chr1_20.rmdup.train.snps.cm.${i}.pca.eigenvec" | \
awk -v OFS="\t" -v nvar=${i} '{print nvar,$0}'; done \
> gwas/grid/genotypes/tau-9/ss500/train/genotypes/revisions/genost9_cm.nvars.eigenvec

for i in {5e3,1e4,5e4,1e5}; \
do tail -n+2 "gwas/grid/genotypes/tau100/ss500/train/genotypes/revisions/genos_gridt100_l1e7_ss750_m0.05_chr1_20.rmdup.train.cm.${i}.pca.eigenvec" | \
awk -v OFS="\t" -v nvar=${i} '{print nvar,$0}'; done \
> gwas/grid/genotypes/tau100/ss500/train/genotypes/revisions/genost100_cm.nvars.eigenvec
