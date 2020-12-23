#!/bin/bash

#script computes true genetic values for all sets (for comparison)

#echo computing for test set
for rep in {1..20}; do \
plink2 --pfile \
../../sibs/genotypes/genos_gridt100_l1e7_ss750_m0.05_chr1_20.rmdup.sib \
--score ../../train/genotypes/genos_gridt100_l1e7_ss750_m0.05_chr1_20.rmdup.train.${rep}.thinned_100kb.effects cols=scoresums \
--out gvalues/p3/gvalue.p3.${rep} \
--score-col-nums 3
done

for rep in {1..20}; do \
tail -n+2 gvalues/p3/gvalue.p3.${rep}.sscore | awk -v OFS="\t" -v i=${rep} '{print i,$0}'; \
done > gvalues/gvalue.p3.all.sscore

#echo computing for training set
for rep in {1..20}; do \
plink2 --pfile \
../../test/genotypes/genos_gridt100_l1e7_ss750_m0.05_chr1_20.rmdup.test \
--score ../../train/genotypes/genos_gridt100_l1e7_ss750_m0.05_chr1_20.rmdup.train.${rep}.thinned_100kb.effects cols=scoresums \
--out gvalues/p2/gvalue.p2.${rep} \
--score-col-nums 3
done

for rep in {1..20}; do \
tail -n+2 gvalues/p2/gvalue.p2.${rep}.sscore | awk -v OFS="\t" -v i=${rep} '{print i,$0}'; \
done > gvalues/gvalue.p2.all.sscore

#echo computing for 'sib' set
for rep in {1..20}; do \
plink2 --pfile \
../../train/genotypes/genos_gridt100_l1e7_ss750_m0.05_chr1_20.rmdup.train \
--score ../../train/genotypes/genos_gridt100_l1e7_ss750_m0.05_chr1_20.rmdup.train.${rep}.thinned_100kb.effects cols=scoresums \
--out gvalues/p1/gvalue.p1.${rep} \
--score-col-nums 3
done
