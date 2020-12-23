
#!/usr/bin/env bash


#carrying out pca on second set of genotypes from the grid (tau=100g) model.


### common-PCA
plink2 \
  --maf 0.05 \
  --out genos_gridt100_l1e7_ss750_m0.05_chr1_20.rmdup.test.cm \
  --pfile genos_gridt100_l1e7_ss750_m0.05_chr1_20.rmdup.test \
  --write-snplist

plink2 \
  --extract genos_gridt100_l1e7_ss750_m0.05_chr1_20.rmdup.test.cm.snplist \
  --out genos_gridt100_l1e7_ss750_m0.05_chr1_20.rmdup.test.cm.200k \
  --pfile genos_gridt100_l1e7_ss750_m0.05_chr1_20.rmdup.test \
  --thin-count 200000 \
  --write-snplist

plink2 \
  --pfile genos_gridt100_l1e7_ss750_m0.05_chr1_20.rmdup.test \
  --extract genos_gridt100_l1e7_ss750_m0.05_chr1_20.rmdup.test.cm.200k.snplist \
  --pca 100 \
  --out genos_gridt100_l1e7_ss750_m0.05_chr1_20.rmdup.test.cm.200k.pca



### rare-PCA
plink2 \
  --mac 2 \
  --max-ac 4 \
  --out genos_gridt100_l1e7_ss750_m0.05_chr1_20.rmdup.test.re.all \
  --pfile genos_gridt100_l1e7_ss750_m0.05_chr1_20.rmdup.test \
  --write-snplist

plink2 \
  --pfile genos_gridt100_l1e7_ss750_m0.05_chr1_20.rmdup.test \
  --extract genos_gridt100_l1e7_ss750_m0.05_chr1_20.rmdup.test.re.all.snplist \
  --pca 100 \
  --out genos_gridt100_l1e7_ss750_m0.05_chr1_20.rmdup.test.re.all.pca
