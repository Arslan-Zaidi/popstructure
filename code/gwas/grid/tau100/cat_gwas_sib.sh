
#!/bin/bash

effects_rep=${1}

for i in {1..20};
do zcat e${effects_rep}/gwas_gridt100_l1e7_ss750_m0.05_chr${i}_e${effects_rep}_n9000_sibs_assort.smooth.glm.linear.gz | \
tail -n+2 -q ; done \
> e${effects_rep}/gwas_gridt100_l1e7_ss750_m0.05_chr1_20_e${effects_rep}_n9000_sibs_assort.smooth.glm.linear


for i in {1..20};
do zcat e${effects_rep}/gwas_gridt100_l1e7_ss750_m0.05_chr${i}_e${effects_rep}_n9000_sibs_assort.sharp.glm.linear.gz | \
tail -n+2 -q ; done \
> e${effects_rep}/gwas_gridt100_l1e7_ss750_m0.05_chr1_20_e${effects_rep}_n9000_sibs_assort.sharp.glm.linear

gzip e${effects_rep}/*.linear
