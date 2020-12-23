#!/bin/bash

#!/bin/bash

#save header to file
head -n1 train/genotypes/genos_grid_d36_m0.07_s500_t9_chr1_20.rmdup.train.snps.cm.200k.pca.eigenvec > cm.header

#cut eigenvec files and remove header
cut -f 3-52 train/genotypes/genos_grid_d36_m0.07_s500_t9_chr1_20.rmdup.train.snps.cm.200k.pca.eigenvec | tail -n+2 - > cm.tmp
cut -f 1-52 train/genotypes/genos_grid_d36_m0.07_s500_t9_chr1_20.rmdup.train.snps.re.1M.pca.eigenvec | tail -n+2 - > re.tmp

#paste them together and add header
paste -d '\t' re.tmp cm.tmp >  cmre.noheader

cat cm.header cmre.noheader > train/genotypes/genos_grid_d36_m0.07_s500_t9_chr1_20.rmdup.train.snps.cmre.eigenvec

rm cmre.noheader
rm cm.header
rm re.tmp
rm cm.tmp
