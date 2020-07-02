#!/bin/bash

#save header to file
head -n1 train/genotypes/genos_complex_l1e7_ss500_m0.08_chr1_20.rmdup.train.cmpruned.pca.eigenvec > cm.header

#cut eigenvec files and remove header
cut -f 3-52 train/genotypes/genos_complex_l1e7_ss500_m0.08_chr1_20.rmdup.train.cmpruned.pca.eigenvec | tail -n+2 - > cm.tmp
cut -f 1-52 train/genotypes/genos_complex_l1e7_ss500_m0.08_chr1_20.rmdup.train.repruned.pca.eigenvec | tail -n+2 - > re.tmp

#paste them together and add header
paste -d '\t' re.tmp cm.tmp >  cmre.noheader

cat cm.header cmre.noheader > train/genotypes/genos_complex_l1e7_ss500_m0.08_chr1_20.rmdup.train.cmre.pca.eigenvec

rm cmre.noheader
rm cm.header
rm re.tmp
rm cm.tmp
