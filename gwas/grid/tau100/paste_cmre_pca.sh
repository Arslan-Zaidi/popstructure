#!/bin/bash

cut -f 3-52 train/plinkPCA/genos_grid_d36_m0.05_s500_t100.rmdup.train.cmpruned.pca.eigenvec | \
paste -d '\t' train/plinkPCA/genos_grid_d36_m0.05_s500_t100.rmdup.train.repruned.pca.eigenvec - \
> train/plinkPCA/genos_grid_d36_m0.05_s500_t100.rmdup.train.cmre.pca.eigenvec
