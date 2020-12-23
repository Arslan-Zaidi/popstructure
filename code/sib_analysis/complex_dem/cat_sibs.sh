#!/bin/bash



for i in {1..4}; \
      do tail -n+2 -q test/sibs/gvalues/genos_complex_l1e7_ss500_m0.12_chr*_i${i}_sibs_random.sscore | \
      awk -v i=$i '{print i,$0}'; \
      done > test/sibs/gvalues/genos_complex_l1e7_ss500_m0.12_call_i1_4_sibs_random.sscore

for i in {1..4}; \
      do tail -n+2 -q test/sibs/gvalues/genos_complex_l1e7_ss500_m0.12_chr*_i${i}_sibs_assort.sscore | \
      awk -v i=$i '{print i,$0}'; \
      done > test/sibs/gvalues/genos_complex_l1e7_ss500_m0.12_call_i1_4_sibs_assort.sscore

for i in {1..4}; \
      do tail -n+2 -q test/sibs/genotypes/genos_complex_assort_${i}.sample | \
      awk -v i=$i '{print i,$0}'; \
      done > test/sibs/genotypes/genos_complex_assort_all.sample

for i in {1..4}; \
      do tail -n+2 -q test/sibs/genotypes/genos_complex_random_${i}.sample | \
      awk -v i=$i '{print i,$0}'; \
      done > test/sibs/genotypes/genos_complex_random_all.sample
