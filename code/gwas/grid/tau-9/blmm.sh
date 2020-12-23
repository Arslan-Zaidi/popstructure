!/bin/bash

phenofile=${1}
phenotype=${2}
out=${3}

bolt --bfile=train/genotypes/beds/genos_grid_d36_m0.07_s500_t9_chr1_20.rmdup.train.snps.b \
--phenoFile=${phenofile} \
--phenoCol=${phenotype} \
--lmm \
--LDscoresFile=train/genotypes/ldsc/genos_grid_d36_m0.07_s500_t9_chr1_20.rmdup.train.snps.ldsc.chr1_20.l2.ldscore.gz \
--statsFile=${out}.cm.stats \
--modelSnps=train/genotypes/genos_grid_d36_m0.07_s500_t9_chr1_20.rmdup.train.snps.cm.200k.list.snplist \
--covarFile=train/genotypes/genos_grid_d36_m0.07_s500_t9_chr1_20.rmdup.train.snps.cm.200k.pca.eigenvec \
--qCovarCol=PC{1:100}

bolt --bfile=train/genotypes/beds/genos_grid_d36_m0.07_s500_t9_chr1_20.rmdup.train.snps.b \
--phenoFile=${phenofile} \
--phenoCol=${phenotype} \
--lmm \
--LDscoresFile=train/genotypes/ldsc/genos_grid_d36_m0.07_s500_t9_chr1_20.rmdup.train.snps.ldsc.chr1_20.l2.ldscore.gz \
--statsFile=${out}.re.stats \
--modelSnps=train/genotypes/genos_grid_d36_m0.07_s500_t9_chr1_20.rmdup.train.snps.re.1M.snplist \
--covarFile=train/genotypes/genos_grid_d36_m0.07_s500_t9_chr1_20.rmdup.train.snps.re.1M.pca.eigenvec \
--qCovarCol=PC{1:100}
