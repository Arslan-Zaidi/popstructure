#!/bin/bash

cut -f 3-52 train/genotypes//genos_grid_d36_m0.07_s500_t9_chr1_20.rmdup.train.cmpruned.pca.eigenvec | \
paste -d '\t' train/genotypes/genos_grid_d36_m0.07_s500_t9_chr1_20.rmdup.train.repruned.pca.eigenvec - \
> train/plinkPCA/genos_grid_d36_m0.05_s500_t100.rmdup.train.cmre.pca.eigenvec
