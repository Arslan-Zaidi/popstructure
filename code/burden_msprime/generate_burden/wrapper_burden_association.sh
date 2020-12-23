#!/bin/bash

source /home/aazaidi/anaconda3/etc/profile.d/conda.sh
conda activate py3_clone

rep=${1}

for i in {1..1000}; \
do bsub -M 16000 -R "rusage[mem=16000]" \
-o logs/gwas/bgwas_r0_t100_x${i}.log -e logs/gwas/bgwas_r0_t100_x${i}.err \
python burden_association.py \
-p tau100/phenotypes/pheno_gridt100_noge_s9k.train.${rep}.txt \
-b tau100/burden/burden_r0_t100_x${i}.haps.npz \
-c tau100/genos_gridt100_l1e7_ss750_m0.05_chr1_20.rmdup.train.cm.200k.eigenvec \
-r tau100/genos_gridt100_l1e7_ss750_m0.05_chr1_20.rmdup.train.re.all.eigenvec \
-o tau100/association/bgwas_r0_t100_x${i}.txt; done


for i in {1..1000}; \
do bsub -M 16000 -R "rusage[mem=16000]" \
-o logs/gwas/bgwas_r1_t100_x${i}.log -e logs/gwas/bgwas_r1_t100_x${i}.err \
python burden_association.py \
-p tau100/phenotypes/pheno_gridt100_noge_s9k.train.${rep}.txt \
-b tau100/burden/burden_r1_t100_x${i}.haps.npz \
-c tau100/genos_gridt100_l1e7_ss750_m0.05_chr1_20.rmdup.train.cm.200k.eigenvec \
-r tau100/genos_gridt100_l1e7_ss750_m0.05_chr1_20.rmdup.train.re.all.eigenvec \
-o tau100/association/bgwas_r1_t100_x${i}.txt; done


for i in {1..1000}; \
do bsub -M 16000 -R "rusage[mem=16000]" \
-o logs/gwas/bgwas_r0_t9_x${i}.log -e logs/gwas/bgwas_r0_t9_x${i}.err \
python burden_association.py \
-p tau-9/phenotypes/pheno_gridt9_noge_s9k.train.${rep}.txt \
-b tau-9/burden/burden_r0_t9_x${i}.haps.npz \
-c tau-9/genos_grid_d36_m0.07_s500_t9_chr1_20.rmdup.train.snps.cm.200k.pca.eigenvec \
-r tau-9/genos_grid_d36_m0.07_s500_t9_chr1_20.rmdup.train.snps.re.1M.pca.eigenvec \
-o tau-9/association/bgwas_r0_t9_x${i}.txt; done


for i in {1..579}; \
do bsub -M 16000 -R "rusage[mem=16000]" \
-o logs/gwas/bgwas_r1_t9_x${i}.log -e logs/gwas/bgwas_r1_t9_x${i}.err \
python burden_association.py \
-p tau-9/phenotypes/pheno_gridt9_noge_s9k.train.${rep}.txt \
-b tau-9/burden/burden_r1_t9_x${i}.haps.npz \
-c tau-9/genos_grid_d36_m0.07_s500_t9_chr1_20.rmdup.train.snps.cm.200k.pca.eigenvec \
-r tau-9/genos_grid_d36_m0.07_s500_t9_chr1_20.rmdup.train.snps.re.1M.pca.eigenvec \
-o tau-9/association/bgwas_r1_t9_x${i}.txt; done
