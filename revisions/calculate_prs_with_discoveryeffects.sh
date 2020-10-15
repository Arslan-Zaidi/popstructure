
#!/bin/bash

for rep in {1..20}; do \
plink2 --pfile \
../../sibs/genotypes/genos_gridt100_l1e7_ss750_m0.05_chr1_20.rmdup.sib \
--score ../../train/betas/est_effects.${rep}.smooth.nc.betas cols=scoresums \
--out prs1sample/prs.smooth.a1_p3.${rep} \
--score-col-nums 3,4,5,6;
done
