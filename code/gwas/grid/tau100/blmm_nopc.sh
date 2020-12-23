#!/bin/bash

phenofile=${1}
phenotype=${2}
out=${3}

bolt --bfile=train/genotypes/beds/genos_gridt100_l1e7_ss750_m0.05_chr1_20.rmdup.train.snps.all \
--phenoFile=${phenofile} \
--phenoCol=${phenotype} \
--lmm \
--LDscoresFile=train/genotypes/ldsc/genos_gridt100_l1e7_ss750_m0.05_chr1_20.rmdup.train.snps.ldsc.chr1_20.l2.ldscore.gz \
--statsFile=${out}.cm.nopc.stats \
--modelSnps=train/genotypes/genos_gridt100_l1e7_ss750_m0.05_chr1_20.rmdup.train.cm.200k.snplist


bolt --bfile=train/genotypes/beds/genos_gridt100_l1e7_ss750_m0.05_chr1_20.rmdup.train.snps.all \
--phenoFile=${phenofile} \
--phenoCol=${phenotype} \
--lmm \
--LDscoresFile=train/genotypes/ldsc/genos_gridt100_l1e7_ss750_m0.05_chr1_20.rmdup.train.snps.ldsc.chr1_20.l2.ldscore.gz \
--statsFile=${out}.re.nopc.stats \
--modelSnps=train/genotypes/genos_gridt100_l1e7_ss750_m0.05_chr1_20.rmdup.train.re.all.snplist
