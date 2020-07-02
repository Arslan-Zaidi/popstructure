#!/bin/bash

#save header to file
head -n1 train/genotypes/genos_gridt100_l1e7_ss750_m0.05_chr1_20.rmdup.train.cm.200k.eigenvec > cm.header

#cut eigenvec files and remove header
cut -f 3-52 train/genotypes/genos_gridt100_l1e7_ss750_m0.05_chr1_20.rmdup.train.cm.200k.eigenvec | tail -n+2 - > cm.tmp
cut -f 1-52 train/genotypes/genos_gridt100_l1e7_ss750_m0.05_chr1_20.rmdup.train.re.all.eigenvec | tail -n+2 - > re.tmp

#paste them together and add header
paste -d '\t' re.tmp cm.tmp >  cmre.noheader

cat cm.header cmre.noheader > train/genotypes/genos_gridt100_l1e7_ss750_m0.05_chr1_20.rmdup.train.cmre.eigenvec

rm cmre.noheader
rm cm.header
rm re.tmp
rm cm.tmp
