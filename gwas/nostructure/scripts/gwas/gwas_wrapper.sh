#!/bin/bash

#carry out GWAS for nongenetic phenotypes

#1. no correction
bash gwas_ge.sh train/genotypes/genos_nostr_s18k_l1e7_chr1_20.rmdup.train \
train/phenotypes/pheno_nostructure_noge_s18k.txt \
train/gwas_results/fixed_effects/noge/gwas_nostr_s18k_noge.train.pcs0

#2. common pca
bash gwas_ge.sh train/genotypes/genos_nostr_s18k_l1e7_chr1_20.rmdup.train \
train/phenotypes/pheno_nostructure_noge_s18k.txt \
train/gwas_results/fixed_effects/noge/gwas_nostr_s18k_noge.train.cm \
train/genotypes/genos_nostr_s18k_l1e7_chr1_20.rmdup.train.cmpruned.pca.eigenvec

#3. rare pca
bash gwas_ge.sh train/genotypes/genos_nostr_s18k_l1e7_chr1_20.rmdup.train \
train/phenotypes/pheno_nostructure_noge_s18k.txt \
train/gwas_results/fixed_effects/noge/gwas_nostr_s18k_noge.train.re \
train/genotypes/genos_nostr_s18k_l1e7_chr1_20.rmdup.train.repruned.pca.eigenvec

## genetic phenotypes
#1. no correction
bash gwas_ge.sh train/genotypes/genos_nostr_s18k_l1e7_chr1_20.rmdup.train \
train/phenotypes/pheno_nostructure_ge_s18k.txt \
train/gwas_results/fixed_effects/ge/gwas_nostr_s18k_ge.train.pcs0

#2. common pca
bash gwas_ge.sh train/genotypes/genos_nostr_s18k_l1e7_chr1_20.rmdup.train \
train/phenotypes/pheno_nostructure_noge_s18k.txt \
train/gwas_results/fixed_effects/ge/gwas_nostr_s18k_ge.train.cm \
train/genotypes/genos_nostr_s18k_l1e7_chr1_20.rmdup.train.cmpruned.pca.eigenvec

#2. common pca
bash gwas_ge.sh train/genotypes/genos_nostr_s18k_l1e7_chr1_20.rmdup.train \
train/phenotypes/pheno_nostructure_noge_s18k.txt \
train/gwas_results/fixed_effects/ge/gwas_nostr_s18k_ge.train.re \
train/genotypes/genos_nostr_s18k_l1e7_chr1_20.rmdup.train.repruned.pca.eigenvec
