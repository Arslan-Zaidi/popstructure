#!/bin/bash

conda activate ldsc

chrom=${1}

ldsc.py --bfile train/chrs/genos_grid_d36_m0.05_s500_t100.rmdup.train.snps.chr${chrom} \
--l2 \
--ld-wind-kb 1000 \
--out train/chrs/chr${chrom}
