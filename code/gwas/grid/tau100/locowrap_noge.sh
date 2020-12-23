#!/bin/bash


#smooth
for j in {1..20}; do \
    bsub -M 16000 -R "rusage[mem=16000]" \
    -o mlma_logs/noge/loco/c${j}.sm.log \
    -e mlma_logs/noge/loco/c${j}.sm.err \
    bash lmmloco_wrapper_gridt100_noge.sh 1 ${j} 1; \
done

#sharp
for j in {1..20}; do \
  bsub -M 16000 -R "rusage[mem=16000]" \
  -o mlma_logs/noge/loco/e${i}_c${j}.sharp.log \
  -e mlma_logs/noge/loco/e${i}_c${j}.sharp.err \
  bash lmmloco_wrapper_gridt100_noge.sh 1 ${j} 3; \
done
