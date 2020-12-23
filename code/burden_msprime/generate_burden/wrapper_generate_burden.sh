#!/bin/bash

source /home/aazaidi/anaconda3/etc/profile.d/conda.sh
conda activate py3_clone

#for tau100, recombination

for i in {1..1000}; \
do bsub -M 16000 -R "rusage[mem=16000]" \
-o logs/generate/gent100.r1.x${i}.log -e logs/generate/gent100.r1.x${i}.err \
python generate_burden_t100.py \
-r 1e-08 \
-L 50000 \
-m 0.05 \
-x ${i} \
-o tau100/burden/burden_r1_t100_x${i}; done

#for tau100, no recombination
for i in {1..1000}; \
do bsub -M 16000 -R "rusage[mem=16000]" \
-o logs/generate/gent100.r0.x${i}.log -e logs/generate/gent100.r0.x${i}.err \
python generate_burden_t100.py \
-r 0 \
-L 1280 \
-m 0.05 \
-x ${i} \
-o tau100/burden/burden_r0_t100_x${i}; done



########### for tau = infinity
#recombination
for i in {1..1000}; \
do bsub -M 16000 -R "rusage[mem=16000]" \
-o logs/generate/gent9.r1.x${i}.log -e logs/generate/gent9.r1.x${i}.err \
python generate_burden_t9.py \
-r 1e-08 \
-L 50000 \
-s 250 \
-m 0.07 \
-x ${i} \
-o tau-9/burden/burden_r1_t9_x${i}; done

#for tau-9, no recombination
for i in {1..1000}; \
do bsub -M 16000 -R "rusage[mem=16000]" \
-o logs/generate/gent9.r0.x${i}.log -e logs/generate/gent9.r0.x${i}.err \
python generate_burden_t9.py \
-r 0 \
-L 1280 \
-s 250 \
-m 0.07 \
-x ${i} \
-o tau-9/burden/burden_r0_t9_x${i}; done
