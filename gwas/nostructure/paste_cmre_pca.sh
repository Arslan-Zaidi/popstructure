#!/bin/bash

cut -f 3-52 train/genotypes/genos_nostr_s18k_l1e7_chr1_20.rmdup.train.cmpruned.pca.eigenvec | \
paste -d '\t' train/genotypes/genos_nostr_s18k_l1e7_chr1_20.rmdup.train.repruned.pca.eigenvec - \
> train/genotypes/genos_nostr_s18k_l1e7_chr1_20.rmdup.train.cmre.pca.eigenvec
