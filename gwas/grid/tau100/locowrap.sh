#!/bin/bash
##07/16/2020
#smooth
for i in {1..10}; do \
  for j in {1..20}; do \
    bsub -M 16000 -R "rusage[mem=16000]" \
    -o mlma_logs/ge/e${i}_c${j}.sm.log \
    -e mlma_logs/ge/e${i}_c${j}.sm.err \
    bash lmmloco_wrapper_gridt100_ge.sh ${i} ${j} 1; \
  done; \
done

#sharp
for i in {1..10}; do \
  for j in {1..20}; do \
    bsub -M 16000 -R "rusage[mem=16000]" \
    -o mlma_logs/ge/e${i}_c${j}.sharp.log \
    -e mlma_logs/ge/e${i}_c${j}.sharp.err \
    bash lmmloco_wrapper_gridt100_ge.sh ${i} ${j} 3; \
  done; \
done
