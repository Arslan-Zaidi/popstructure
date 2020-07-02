#!/bin/bash

source /home/aazaidi/anaconda3/etc/profile.d/conda.sh
conda activate rpkgs

rep=${1}

bash scripts/prs/cal_prs.sh smooth ${rep}

#bash scripts/prs/cal_prs.sh smooth_long ${rep}

bash scripts/prs/cal_prs.sh sharp ${rep}

#bash scripts/prs/cal_prs.sh grandom ${rep}

plink2 --pfile test/genotypes/genos_gridt100_l1e7_ss750_m0.05_chr1_20.rmdup.test \
--score train/genotypes/genos_gridt100_l1e7_ss750_m0.05_chr1_20.rmdup.train.${rep}.thinned_100kb.effects cols=dosagesum,scoresums \
--out test/gvalue/genos_grid_d36_m0.05_s500_t100.rmdup.test.${rep}.gvalue
