#!/bin/bash

source /home/aazaidi/anaconda3/etc/profile.d/conda.sh
conda activate py3_clone


migrate=${1}
tau=${2}
chr=${3}

#simulate genotypes and output vcf
python generate_genos_grid.py \
--sample_size 250 \
--tmove ${tau} \
--chr 1 \
--length 1000000 \
--migrate ${migrate} \
--outpre genotypes/genos_grid_t${tau}_l1e7_ss250_m${migrate}_chr${chr}

#gzip vcf file
#gzip genotypes/genos_grid_t${tau}_l1e7_ss250_m${migrate}_chr${chr}.vcf

#add variant IDs
cat genotypes/genos_grid_t${tau}_l1e7_ss250_m${migrate}_chr${chr}.vcf | head -n6 > genotypes/header.txt

cat genotypes/header.txt <(cat genotypes/genos_grid_t${tau}_l1e7_ss250_m${migrate}_chr${chr}.vcf | \
awk -v OFS="\t" 'NR>6 {$3=$1"_"$2"_A_T"; $4="A"; $5="T"; print ;}') | \
gzip > genotypes/genos_grid_t${tau}_l1e7_ss250_m${migrate}_chr${chr}.vcf.gz

rm genotypes/genos_grid_t${tau}_l1e7_ss250_m${migrate}_chr${chr}.vcf

#convert vcf file to plink format
plink2 \
--vcf genotypes/genos_grid_t${tau}_l1e7_ss250_m${migrate}_chr${chr}.vcf.gz \
--make-pgen \
--double-id \
--out genotypes/genos_grid_t${tau}_l1e7_ss250_m${migrate}_chr${chr}

#carry out birthplace GWAS
plink2 \
--mac 1 \
--pfile genotypes/genos_grid_t${tau}_l1e7_ss250_m${migrate}_chr${chr} \
--glm hide-covar \
--pheno genos_grid_t${tau}_d36_s250.pop \
--out bplace_gwas/genos_grid_t${tau}_l1e7_ss250_m${migrate}_chr${chr}

#calculate fst across demes using plink
plink2 \
--pfile genotypes/genos_grid_t${tau}_l1e7_ss250_m${migrate}_chr${chr} \
--maf 0.05 \
--indep-pairwise 100 10 0.1 \
--write-snplist \
--out genotypes/genos_grid_t${tau}_d36_s250_m${migrate}_chr${chr}.cmpruned

plink2 \
--pfile genotypes/genos_grid_t${tau}_l1e7_ss250_m${migrate}_chr${chr} \
--extract genotypes/genos_grid_t${tau}_d36_s250_m${migrate}_chr${chr}.cmpruned.prune.in \
--make-bed --out genotypes/genos_grid_t${tau}_l1e7_ss250_m${migrate}_chr${chr}.cmpruned.b

plink --bfile genotypes/genos_grid_t${tau}_l1e7_ss250_m${migrate}_chr${chr}.cmpruned.b \
--fst \
--within genos_grid_t${tau}_d36_s250.pop \
--out fst/genos_grid_t${tau}_l1e7_ss250_m${migrate}_chr${chr}.cmpruned.b.fst

#calculate fst using scikit allel
plink2 \
--pfile genotypes/genos_grid_t${tau}_l1e7_ss250_m${migrate}_chr${chr} \
--extract genotypes/genos_grid_t${tau}_d36_s250_m${migrate}_chr${chr}.cmpruned.prune.in \
--export vcf-4.2 id-delim="-" \
--out genotypes/genos_grid_t${tau}_l1e7_ss250_m${migrate}_chr${chr}.cmpruned

gzip genotypes/genos_grid_t${tau}_l1e7_ss250_m${migrate}_chr${chr}.cmpruned.vcf

python cal_fst.py \
-v genotypes/genos_grid_t${tau}_l1e7_ss250_m${migrate}_chr${chr}.cmpruned.vcf.gz \
-o fst/genos_grid_t${tau}_l1e7_ss250_m${migrate}_chr${chr}.cmpruned
