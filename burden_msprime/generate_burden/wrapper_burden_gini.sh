#!/bin/bash

source /home/aazaidi/anaconda3/etc/profile.d/conda.sh
conda activate py3_clone

for i in {120..1000}; \
do bsub -M 16000 -R "rusage[mem=16000]" \
python burden_gini.py \
-b tau100/burden/burden_r0_t100_x${i}.haps.npz \
-p tau100/iid_train.txt \
-o tau100/gini/bgini_r0_t100_x${i}.txt; done

for i in {1..100}; \
do python burden_gini.py \
-b tau100/burden/burden_r1_t100_x${i}.haps.npz \
-p tau100/iid_train.txt \
-o tau100/gini/bgini_r1_t100_x${i}.txt; done

for i in {1..100}; \
do python burden_gini.py \
-b tau-9/burden/burden_r0_t9_x${i}.haps.npz \
-p tau-9/genos_grid_d36_m0.07_s500_t9.train.pop \
-o tau-9/gini/bgini_r0_t9_x${i}.txt; done

for i in {1..100}; \
do python burden_gini.py \
-b tau-9/burden/burden_r1_t9_x${i}.haps.npz \
-p tau-9/genos_grid_d36_m0.07_s500_t9.train.pop \
-o tau-9/gini/bgini_r1_t9_x${i}.txt; done
