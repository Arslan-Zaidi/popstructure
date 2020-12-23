#!/bin/bash

source /home/aazaidi/anaconda3/etc/profile.d/conda.sh
conda activate rpkgs

rep=${1}

bash scripts/prs/cal_prs.sh smooth ${rep} weighted

#bash scripts/prs/cal_prs.sh smooth_long ${rep}

bash scripts/prs/cal_prs.sh sharp ${rep} weighted

#bash scripts/prs/cal_prs.sh grandom ${rep}

plink2 --pfile test/genotypes/genos_ukb_l1e7_ss500_m0.08_weighted_chr1_20.rmdup.test \
--score train/genotypes/genos_ukb_l1e7_ss500_m0.08_weighted_chr1_20.rmdup.train.${rep}.thinned_100kb.effects cols=dosagesum,scoresums \
--out test/gvalue/genos_ukb_l1e7_ss500_m0.08_weighted_chr1_20.rmdup.train.${rep}.gvalue
