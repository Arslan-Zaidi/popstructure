#!/bin/bash

#conda activate py3_clone


migrate=${1}
tau=${2}

#simulate genotypes and output vcf
python generate_genos_grid.py \
--sample_size 500 \
--tmove ${tau} \
--chr 1 \
--migrate ${migrate} \
--outpre genotypes/genos_grid_t${tau}_l1e7_ss500_m${migrate}

#gzip vcf file
#gzip genotypes/genos_grid_t${tau}_l1e7_ss500_m${migrate}_chr1.vcf

#add variant IDs
cat genotypes/genos_grid_t${tau}_l1e7_ss500_m${migrate}_chr1.vcf | head -n6 > genotypes/header.txt

cat genotypes/header.txt <(cat genotypes/genos_grid_t${tau}_l1e7_ss500_m${migrate}_chr1.vcf | \
awk -v OFS="\t" 'NR>6 {$3=$1"_"$2"_A_T"; $4="A"; $5="T"; print ;}') | \
gzip > genotypes/genos_grid_t${tau}_l1e7_ss500_m${migrate}_chr1.vcf.gz

rm genotypes/genos_grid_t${tau}_l1e7_ss500_m${migrate}_chr1.vcf

#convert vcf file to plink format
plink2 \
--vcf genotypes/genos_grid_t${tau}_l1e7_ss500_m${migrate}_chr1.vcf.gz \
--make-pgen \
--double-id \
--out genotypes/genos_grid_t${tau}_l1e7_ss500_m${migrate}_chr1

#carry out birthplace GWAS
plink2 \
--mac 1 \
--pfile genotypes/genos_grid_t${tau}_l1e7_ss500_m${migrate}_chr1 \
--glm hide-covar \
--pheno genos_grid_t${tau}_d36_s500.pop \
--out gwas_results/genos_grid_t${tau}_l1e7_ss500_m${migrate}_chr1

#calculate fst across demes using plink
plink2 \
--pfile genotypes/genos_grid_t${tau}_l1e7_ss500_m${migrate}_chr1 \
--maf 0.05 \
--indep-pairwise 100 10 0.1 \
--write-snplist \
--out genotypes/genos_grid_t${tau}_d36_s500_m${migrate}_chr1.cmpruned

plink2 \
--pfile genotypes/genos_grid_t${tau}_l1e7_ss500_m${migrate}_chr1 \
--extract genotypes/genos_grid_t${tau}_d36_s500_m${migrate}_chr1.cmpruned.prune.in \
--make-bed --out genotypes/genos_grid_t${tau}_l1e7_ss500_m${migrate}_chr1.cmpruned.b

plink --bfile genotypes/genos_grid_t${tau}_l1e7_ss500_m${migrate}_chr1.cmpruned.b \
--fst \
--within genos_grid_t${tau}_d36_s500.pop \
--out fst/genos_grid_t${tau}_l1e7_ss500_m${migrate}_chr1.cmpruned.b.fst

#calculate fst using scikit allel
plink2 \
--pfile genotypes/genos_grid_t${tau}_l1e7_ss500_m${migrate}_chr1 \
--extract genotypes/genos_grid_t${tau}_d36_s500_m${migrate}_chr1.cmpruned.prune.in \
--export vcf-4.2 id-delim="-" \
--out genotypes/genos_grid_t${tau}_l1e7_ss500_m${migrate}_chr1.cmpruned

gzip genotypes/genos_grid_t${tau}_l1e7_ss500_m${migrate}_chr1.cmpruned.vcf

python cal_fst.py \
-v genotypes/genos_grid_t${tau}_l1e7_ss500_m${migrate}_chr1.cmpruned.vcf.gz \
-o fst/genos_grid_t${tau}_l1e7_ss500_m${migrate}_chr1.cmpruned
