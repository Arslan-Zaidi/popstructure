#!/bin/bash

##07/16/2020
#bug record: The heritable phenotype file has phenotypes in the following order:
# 1. grandom
# 2. smooth
# 3. smooth_long
# 4. sharp

#I accidentally used column 1 for smooth and 3 for sharp for gcta-loco on polygenic traits.
#I should have used 2 and 4. Fortunately, 3 is smooth_long, which is the same as smooth except in the longitudinal direction.
#To correct for this, I will re-run prs for the sharp effect by changing the column number below from 3 to 4.

#To correct for the smooth effect, I will simply rename all the prs files and gwas files to smooth_long and plot the prs in the other direction
#making a note here to avoid confusion later

#smooth (leaving unchanged)
# for i in {1..10}; do \
#   for j in {1..20}; do \
#     bsub -M 16000 -R "rusage[mem=16000]" \
#     -o mlma_logs/ge/e${i}_c${j}.sm.log \
#     -e mlma_logs/ge/e${i}_c${j}.sm.err \
#     bash lmmloco_wrapper_gridt100_ge.sh ${i} ${j} 1; \
#   done; \
# done

#sharp (changing column number from 3 to 4) and rerun
for i in {1..10}; do \
  for j in {1..20}; do \
    bsub -M 16000 -R "rusage[mem=16000]" \
    -o mlma_logs/sharp/e${i}_c${j}.sharp.log \
    -e mlma_logs/sharp/e${i}_c${j}.sharp.err \
    bash lmmloco_wrapper_gridt100_ge.sh ${i} ${j} 4; \
  done; \
done
